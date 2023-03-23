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
  for (i, p) in enumerate(pos)
    idx = findfirst(lbls .== annotationimgw[p...])
    c[idx] = c[idx].+1
    push!(a[idx], amp[i])
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

