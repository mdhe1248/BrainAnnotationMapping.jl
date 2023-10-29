struct BlobPos
  location::Tuple{Float64, Float64}
  σ
  amplitude::Float64
end

## visualize moving images
function pad_images(vector_of_images; h = :auto, w = :auto)
  if h == :auto
    h = max([size(vector_of_images[i],1) for i in eachindex(vector_of_images)]...)
  end
  if w == :auto
    w = max([size(vector_of_images[i],2) for i in eachindex(vector_of_images)]...)
  end
  paddedimages= [PaddedView(0, vector_of_images[i], (h, w), padOrigin((h, w), vector_of_images[i])) for i in eachindex(vector_of_images)]
  return(paddedimages)
end

#### Warp the fixed and annotation images
function warp_reference(outdir, fixed, annotationfn, slices0, fx_pxspacing, tfm)
  ## load annotation image
  annotationimg = load(annotationfn)
  annotationimg = setAxis(parent(annotationimg), fx_pxspacing)

  fixedw = warp(fixed, tfm, 0) #fill value is 0. If NaNs exist, antsRegistration does not work.
  annotationw = warp(annotationimg, tfm, fillvalue = 0, method = BSpline(Constant()))
  slices = slices0.-Base.axes(fixedw)[3].offset #offset to make "1" as the first frame.
  fixed2ds = [parent(fixedw)[:,:,slice] for slice in slices]
  
  #### select the slices matching to the moving images
  fixed2ds = map(x -> setAxis(parent(x), fx_pxspacing), fixed2ds) #Assign axes
  annotation2ds = [parent(annotationw)[:,:,slice] for slice in slices]
  annotation2ds = map(x->setAxis(parent(x), fx_pxspacing), annotation2ds) #Assign axes
  
  #### Save
  ## save file names
  fixed_tform_savefn = outdir*"fixed_tform.jld2"
  fixed2d_savefns = string.(outdir, "fixed2d_", slices ,".nrrd") #save filename
  annotation2d_savefns = string.(outdir, "annotation2d_", slices ,".nrrd") #save filename
  save.(fixed2d_savefns, fixed2ds)
  save.(annotation2d_savefns, annotation2ds)
  save(fixed_tform_savefn, Dict("x_rot"=>x_rot, "y_rot"=>y_rot, "z_rot"=>z_rot, "tform" =>tfm, "slices"=>slices))
end


"""
#### Bouondary detection
"""
function annotation_boundary(annotation_img)
  parent_img = copy(parent(annotation_img))
  output_img = zeros(size(parent_img))
  for i in Base.axes(output_img,3)
    im1 = diff(parent_img[:,:,i], dims = 1)
    output_img[2:end, :, i] .+= im1
    im2 = diff(parent_img[:,:,i], dims = 2)
    output_img[:, 2:end, i] .+= im2
  end
  output_img[output_img .!= 0] .= 1
  return(output_img)
end

"""
append a column.

`lbl` should be the id of brain structure.
`dat` is corresponding value.
empty values will have `missing`.
"""
function appendCol!(structures, lbls, dat, newcolname::Symbol)
  if length(lbls) .!= length(dat)
    warn("Length of `lbls` and `dat` must be the same.")
  end
  v = fill(missing, size(structures,1))
  v = convert(Vector{Union{Missing, Float64}}, v)
  for (i, l) in enumerate(lbls)
    idx = findfirst(structures.id .== l)
    if isnothing(idx)
      println(string(l, " does not exist."))
    else
      v[idx] = dat[i]
    end
  end
  structures[!, newcolname] = v
end

"""
append a column of Vector{Float64}.

`lbl` should be the id of brain structure.
`dat` is corresponding value.
empty values will have `missing`.
"""
function appendVecCol!(structures, lbls, dat, newcolname::Symbol)
  if length(lbls) .!= length(dat)
    warn("Length of `lbls` and `dat` must be the same.")
  end
  v = fill(missing, size(structures,1))
  v = convert(Vector{Union{Missing, Vector{Float64}}}, v)
  for (i, l) in enumerate(lbls)
    idx = findfirst(structures.id .== l)
    if isnothing(idx)
      println(string(l, " does not exist."))
    else
      v[idx] = dat[i]
    end
  end
  structures[!, newcolname] = v
end

""" get cFos counts"""
function getPosCounts(annotationimgw, pos, amp) 
  ## Load image
  lbls = Integer.(unique(annotationimgw))
  c = zeros(Int, length(lbls))
  a = [Vector{Float64}(undef,0) for i in eachindex(lbls)]
  n = 0 # to count the number of skipped blobs.
  for (i, p) in enumerate(pos)
    if minimum(p) < 1 #If the coordinate is smaller then 1, skip.
      n += 1
      println(string("Negative coordinate skipped. (", n, ")"))
    else 
      idx = findfirst(lbls .== annotationimgw[p...])
      c[idx] = c[idx].+1
      push!(a[idx], amp[i])
    end
  end
  return(lbls, c, a)
end

""" get mean intensities in each brain area"""
function getIntensityMean(annotationimgw, img)
  ## Load image
  lbls = Integer.(unique(annotationimgw))
  meanintensities = [mean(img[annotationimgw .== lbl]) for lbl in lbls]
  return(lbls, meanintensities)
end

function getIntensityMean(annotationimgw, img, thresh)
  ## Load image
  lbls = Integer.(unique(annotationimgw))
  meanintensities = [mean(img[annotationimgw .== lbl .&& img .>= thresh]) for lbl in lbls]
  return(lbls, meanintensities)
end

""" get number of pixels of each brain area"""
function getNPixels(annotationimgw, img)
  ## Load image
  lbls = Integer.(unique(annotationimgw))
  npixels = [length(img[annotationimgw .== lbl]) for lbl in lbls]
  return(lbls, npixels)
end

""" overlap img and annotation boundary
`c` is contrast limit
"""
function imshow_attn(annotationimgw, img, c)
  ## Overage an image and annotation image
  boundaryimg = annotation_boundary(annotationimgw)
  scalefun = scaleminmax(c[1], c[2])
  tmp = scalefun.(img)
  tmp[boundaryimg .== 1] .= 1
  guidict = ImageView.imshow(tmp)
  guidict
end

""" select specific area"""
function selectArea(attnimg, targetid)
  output = copy(parent(attnimg))
  output[output .!= targetid] .= 0
  output
end

""" blob amplitude filter"""
function blob_filter(blobs_log, thresh)
  blobs_filtered = Vector{typeof(blobs_log[1])}(undef, 0)
  for blob in blobs_log
    if blob.amplitude > thresh
      push!(blobs_filtered, blob)
    end
  end
  blobs_filtered
end

""" Find nearest """
findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]

"""
Image filtering and blob detection
"""
function detect_blobs(movingfn, ch, thresh_slope; show_threshold = false, show_blobs = false)
  ## Load images
  img = load(movingfn)
  imgc = img[:,:, ch] #cFos image
  blobs_filtered = detect_blobs(imgc, thresh_slope; show_threshold = show_threshold, show_blobs = show_blobs)
  return(blobs_filtered)
end

function detect_blobs(imgc::AbstractMatrix, thresh_slope; show_threshold = false, show_blobs = false)  
  ## background subtraction
  imgb = erode(imgc) # minimum filter
  imgb = erode(imgb) # minimum filter
  imgb = imfilter(imgb, centered(ones(5,5)./25)) #mean filter
  imgsig = imgc .- imgb 
  
  ## blob detection
  imgf = imfilter(imgsig, Kernel.gaussian(1))
  blobs_log = blob_LoG(imgf, [1]); #Dence detection
  
  ## thresholding
  inten = map(x -> x.amplitude, blobs_log)
  inten_sort = sort(inten)
  kernel = OffsetArray(fill(1/1001, 1001), -500:500) # moving average filter
  var = imfilter(inten_sort, kernel)
  idx = findnearest(diff(var), thresh_slope)
  thresh = inten_sort[idx] ## Threshold by slope

  #### threshold check
  if show_threshold
    fig = figure(); fig.add_subplot(2,1,1); plot(diff(var), eachindex(diff(var)))
    axvline(thresh_slope); title(string("Slope ", thresh_slope))
    ylabel("Counts"); xlabel("Slope")
    fig.add_subplot(2,1,2); plot(inten_sort, eachindex(inten)); 
    axvline(thresh); title(string("Amplitude"))
    ylabel("Counts"); xlabel("Amplitude"); tight_layout();
  end

  ## Filter by threshold
  blobs_filtered = blob_filter(blobs_log, thresh)

  #### Visualize
  if show_blobs
    guidict = ImageView.imshow(imgc, CLim(0, 0.03));
    for blob in blobs_filtered
      y, x = blob.location[1], blob.location[2]
      idx = annotate!(guidict, AnnotationText(x, y, "x", color=RGB(1,0,0), fontsize=10))
    end
  end
  return(blobs_filtered)
end

"""coordinate scaling in physical space"""
function scale_pos(blobs_filtered, mv_pxspacing_midres, xoffset, yoffset)
  s =  mv_pxspacing_midres #scale
  blobs_scaled = Vector{BlobPos}(undef, length(blobs_filtered))
  for (i, b) in enumerate(blobs_filtered)
    blobs_scaled[i] = BlobPos((s[1].val*b.location[1]+yoffset, s[2].val*b.location[2]+xoffset), b.σ, b.amplitude)
  end
  return(blobs_scaled)
end

"""blob position to data frame for antregistration point transformation
if z is not provided, default z is 0."""
function blobs2df(blobs_scaled)
  df_pos, df_amp = blobs2df(blobs_scaled, 0)
  df_amp = DataFrame(amplitude = map(x -> x.amplitude, blobs_scaled)) #data frame
  return(df_pos, df_amp)
end

function blobs2df(blobs_scaled, z::Int)
  df_pos = DataFrame(x = map(x -> x.location[1], blobs_scaled), y = map(x -> x.location[2], blobs_scaled), z = ones(Int, length(blobs_scaled)).*z, t = zeros(Int, length(blobs_scaled))) #data frame
  df_amp = DataFrame(amplitude = map(x -> x.amplitude, blobs_scaled)) #data frame
  return(df_pos, df_amp)
end

""" Load images with a vector of image file names `imgfns`"""
function loadimgs(imgfns)
  img = load(imgfns[1])
  imgvec = Vector{typeof(img)}(undef, length(imgfns))
  imgvec[1] = img
  for i in 2:length(imgfns)
    img = load(imgfns[i])
    imgvec[i] = img
  end
  return(imgvec)
end

"""Load coordinates at the fixed image physical scales"""
function load_pos_in_physical_space(ptsw_pos_savefn, fx_pxspacing)
  ptsw = CSV.read(ptsw_pos_savefn, DataFrame) #points warped in physical space
  df = DataFrame(:x => round.(Int, ptsw.x/fx_pxspacing[1].val), :y => round.(Int, ptsw.y/fx_pxspacing[1].val), :z =>0, :t => 0)
  return(df)
end
load_pos_in_physical_space(ptsw_pos_savefns::Vector{String}, fx_pxspacing) = [load_pos_in_physical_space(ptsw_pos_savefn, fx_pxspacing) for ptsw_pos_savefn in ptsw_pos_savefns]

""" Load multiple dataframe files """
loaddfs(fns::Vector{String}) = [CSV.read(fn, DataFrame) for fn in fns]

""" from a collection of multi channel images, select one channel and concatenate the frames to generate 3d image)"""
function imcat(imgwcat::AbstractVector, channel::Number)
  sz = size(imgwcat[1])
  imcat1 = zeros(eltype(imgwcat[1]), sz[1], sz[2], length(imgwcat));
  for i in 1:size(imcat1,3)
    imcat1[:,:,i] = imgwcat[i][:,:,channel]
  end
  imcat1
end

