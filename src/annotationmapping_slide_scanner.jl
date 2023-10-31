struct BlobPos
  location::Tuple{Float64, Float64}
  σ
  amplitude::Float64
end

mutable struct BlobVars
  outdir::String
  movingfn::String
  mv_pxspacing::NTuple{2, Number} #pixel spacing (mid resolution)
  thresh_slope::Float64 #cfos detection threshold
  cfos_channel::Int #cfos channel
  xoffset::Number #in micrometer
  yoffset::Number
  σ::Vector
  blobvars_fn::String
  pts_pos_scaled_savefn::String
  pts_amp_savefn::String
  ptsw_pos_savefn::String
end
BlobVars(outdir, movingfn, mv_pxspacing, thresh_slope, cfos_channel, xoffset, yoffset, σ) = BlobVars(outdir, movingfn, mv_pxspacing, thresh_slope, cfos_channel, xoffset, yoffset, σ
  outdir*"blobvars_"*first(splitext(last(splitdir(movingfn))))[end-1:end]*".jld2", #regvar fn
  string(outdir, first(splitext(last(splitdir(movingfn)))), "_cfos_points.csv"),
  string(outdir, first(splitext(last(splitdir(movingfn)))), "_cfos_amplitude.csv"),
  string(outdir, first(splitext(last(splitdir(movingfn)))), "_cfos_points_tform.csv"))

function save_blobvars(var)
  jldsave(var.blobvars_fn, blobvars = var)
end
save_blobvars(vars::AbstractVector) = [save_blobvars(var) for var in vars]

assign_blobvars(outdir, movingfns::Vector, mv_pxspacing_midres, thresh_slope, cfos_channel, xoffset, yoffset, σ) = [BlobVars(outdir, movingfn, mv_pxspacing_midres, thresh_slope, cfos_channel, xoffset, yoffset, σ) for movingfn in movingfns]

function save_blobs_in_physical_space(blobs::Vector, blobvars::BlobVars; show_threshold = false, show_blobs = false)
  blobs_scaled = pos_in_physical_space(blobs, blobvars)
  df_pos, df_amp = blobs2df(blobs_scaled)
  CSV.write(blobvars.pts_pos_scaled_savefn, df_pos)
  CSV.write(blobvars.pts_amp_savefn, df_amp) #save as CSV
  println("Saved.")
end

function runAntsApplyTransformsToPoints(regvar, blobvar)
  warped_nrrd_fns = regvar.warpout_fn
  tform1_fn = regvar.tform1_fn
  inverse_warp_fn = regvar.inverse_warp_fn
  pts_fn = blobvar.pts_pos_scaled_savefn
  ptsw_fn = blobvar.ptsw_pos_savefn
  antsApplyTransformsToPoints = `antsApplyTransformsToPoints -d 2 -i $pts_fn -o $ptsw_fn -t \[$tform1_fn, 1\] -t $inverse_warp_fn` # For points, transformation is inversed
  run(antsApplyTransformsToPoints) #FYI, data will be saved
end

function imshow_blobs(img_vec::Vector, blobs_vec::Vector{<:Vector{<:BlobLoG}}, clim; size = 1, scale = true)
  img = cat(img_vec..., dims = 3)
  guidict = ImageView.imshow(img, CLim(clim...));
  for (i, blobs) in enumerate(blobs_vec)
    pos = Vector{NTuple{2, Real}}()
    for blob in blobs
      push!(pos, (blob.location[2], blob.location[1]))
    end
    idx = annotate!(guidict, AnnotationPoints(pos, shape = '.', size=size, color=RGB(1,0,0), z = i, scale = scale))
  end
  return(guidict)
end
imshow_blobs(img1, blobs::Vector{<:BlobLoG}, clim; size = 1, scale = true) = imshow_blobs([img1], [blobs], clim; size = size, scale = scale)
function imshow_blobs(img_vec::Vector, blobs_vec::Vector{<:DataFrame}, clim; size = 1, scale = true)
  img = cat(img_vec..., dims = 3)
  guidict = ImageView.imshow(img, CLim(clim...));
  for (i, blobs) in enumerate(blobs_vec)
    pos = Vector{NTuple{2, Real}}()
    for blob in eachrow(blobs)
      push!(pos, (blob[2], blob[1]))
    end
    idx = annotate!(guidict, AnnotationPoints(pos, shape = '.', size=size, color=RGB(1,0,0), z = i, scale = scale))
  end
  return(guidict)
end
imshow_blobs(img1, blobs::DataFrame, clim; size = 1, scale = true) = imshow_blobs([img1], [blobs], clim; size = size, scale = scale)

function imshow_blobs(imgfns::Vector{<:String}, blob_fns::Vector{<:String}, clim, fx_pxspacing; size = 1, scale = true)
  img_vec = load.(imgfns)
  blobs_vec = load_pos_in_index_space(blob_fns, fx_pxspacing)
  guidict = imshow_blobs(blobs_vec, img_vec, clim; size = size, scale = scale)
  return(guidict)
end
imshow_blobs(imgfn::String, blob_fn::String, clim, fx_pxspacing; size = 1, scale = true) = imshow_blobs([imgfn], [blob_fn], clim, fx_pxspacing; size = size, scale = scale)


#function imshow_blobs(blobs::Vector, img1, fontsize, clim)
#  guidict = ImageView.imshow(img1, CLim(clim...));
#  for blob in blobs
#    y, x = blob.location[1], blob.location[2]
#    idx = annotate!(guidict, AnnotationText(x, y, "x", color=RGB(1,0,0), fontsize=fontsize))
#  end
#  return(guidict)
#end
#
#function imshow_blobs(blobs::DataFrame, img1, fontsize, clim)
#  guidict = ImageView.imshow(img1, CLim(clim...));
#  for blob in eachrow(blobs)
#    y, x = blob[1], blob[2]
#    idx = annotate!(guidict, AnnotationText(x, y, "x", color=RGB(1,0,0), fontsize=fontsize))
#  end
#  return(guidict)
#end

function overlay_boundary(var::Regvars, clim)
  img = load(var.warpedfn)
  attnimg = load(var.annotation2d_fn)
  boundaryimg = overlay_boundary(img, attnimg, clim)
  return(boundaryimg)
end

function overlay_boundary(img, attnimg, clim)
  tmp = annotation_boundary(attnimg)
  scalefun = scaleminmax(clim...) #contrast
  boundaryimg = scalefun.(img)
  boundaryimg[tmp .== 1] .= 1
  return(boundaryimg)
end

function overlay_boundary(vars::Vector{Regvars}, clim)
  boundaryimgs = Vector{Matrix}(undef, length(vars));
  for i in eachindex(boundaryimgs)
    boundaryimgs[i] = overlay_boundary(vars[i], clim)
  end
  return(cat(boundaryimgs..., dims = 3))
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
  save(fixed_tform_savefn, Dict("x_rot"=>rad2deg(tfm.linear.theta1), "y_rot"=>rad2deg(tfm.linear.theta2), "z_rot"=>rad2deg(tfm.linear.theta3), "tform" =>tfm, "slices"=>slices))
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
  attnimgsz = size(annotationimgw)
  lbls = Integer.(unique(annotationimgw))
  c = zeros(Int, length(lbls))
  a = [Vector{Float64}(undef,0) for i in eachindex(lbls)]
  n = 0 # to count the number of skipped blobs.
  for (i, p) in enumerate(pos)
    if minimum(p) < 1 || max(p, attnimgsz) > attnimgsz #If the coordinate is smaller then 1, skip.
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

function detect_blobs(movingfn, ch, thresh_slope, σ; show_threshold = false, show_blobs = false, ptsize = 0.5, clim = (0, 0.05))
  ## Load images
  img = load(movingfn)
  imgc = img[:,:, ch] #cFos image
  blobs_filtered = detect_blobs(imgc, thresh_slope, σ; show_threshold = show_threshold, show_blobs = show_blobs, ptsize = ptsize, clim = clim)
  println("Done.")
  return(blobs_filtered)
end

function detect_blobs(imgc::AbstractMatrix, thresh_slope, σ; show_threshold = false, show_blobs = false, ptsize = 0.5, clim = (0, 0.05))
  ## background subtraction
  imgb = erode(imgc) # minimum filter
  imgb = erode(imgb) # minimum filter
  imgb = imfilter(imgb, centered(ones(5,5)./25)) #mean filter
  imgsig = imgc .- imgb 
  
  ## blob detection
  imgf = imfilter(imgsig, Kernel.gaussian(1))
  blobs_log = blob_LoG(imgf, σ); #Dence detection
  
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
    imshow_blobs(imgc, blobs_filtered, clim; size = ptsize);
  end
  return(blobs_filtered)
end

detect_blobs(blobvars; show_threshold = false, show_blobs = false, ptsize = 0.5, clim = (0, 0.05)) = detect_blobs(blobvars.movingfn, blobvars.cfos_channel, blobvars.thresh_slope, blobvars.σ; show_threshold = show_threshold, show_blobs = show_blobs, ptsize = ptsize, clim = clim)

"""coordinate scaling in physical space"""
function pos_in_physical_space(blobs_filtered, mv_pxspacing_midres, xoffset, yoffset)
  s =  mv_pxspacing_midres #scale
  blobs_scaled = Vector{BlobPos}(undef, length(blobs_filtered))
  for (i, b) in enumerate(blobs_filtered)
    blobs_scaled[i] = BlobPos((s[1].val*b.location[1]+yoffset, s[2].val*b.location[2]+xoffset), b.σ, b.amplitude)
  end
  return(blobs_scaled)
end

pos_in_physical_space(blobs_filtered, blobvars) = pos_in_physical_space(blobs_filtered, blobvars.mv_pxspacing, blobvars.xoffset, blobvars.yoffset)

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
function load_pos_in_index_space(ptsw_pos_savefn, fx_pxspacing)
  ptsw = CSV.read(ptsw_pos_savefn, DataFrame) #points warped in physical space
  df = DataFrame(:x => round.(Int, ptsw.x/fx_pxspacing[1].val), :y => round.(Int, ptsw.y/fx_pxspacing[1].val), :z =>0, :t => 0)
  return(df)
end
load_pos_in_index_space(ptsw_pos_savefns::Vector{String}, fx_pxspacing) = [load_pos_in_index_space(ptsw_pos_savefn, fx_pxspacing) for ptsw_pos_savefn in ptsw_pos_savefns]

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

