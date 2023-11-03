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

""" get number of pixels of each brain area"""
function getNPixels(annotationimgw, img)
  ## Load image
  lbls = Integer.(unique(annotationimgw))
  npixels = [length(img[annotationimgw .== lbl]) for lbl in lbls]
  return(lbls, npixels)
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

function addIntensities(structures::DataFrame, regvars::Vector{<:Regvars}, blobvars::Vector{<:BlobVars}, fx_pxspacing, pixel_intensity_thresh; normalize = true)
  newstructures = copy(structures)
  imgwcat = load.(map(x -> x.warpout_fn, regvars))
  attncat = load.(map(x -> x.annotation2d_fn, regvars))
  fixedcat = load.(map(x -> x.fixed2d_fn, regvars))
  ptsw_coords = load_pos_in_index_space(map(x -> x.ptsw_pos_savefn, blobvars), fx_pxspacing)
  ptsw_amps = CSV.read.(map(x -> x.pts_amp_savefn, blobvars), DataFrame)
  for imgidx in eachindex(imgwcat)
    if normalize
    ##Mean intensity
      imgn = imgwcat[imgidx][:,:,blobvars[imgidx].cfos_channel]./imgwcat[imgidx][:,:,regvars[imgidx].bg_channel] # Image normalization (mcherry / gfp)
    else
      imgn = img;
    end
    imgn[isnan.(imgn)] .= 0
    lbls, meanintensities = getIntensityMean(attncat[imgidx], imgn, pixel_intensity_thresh) ##### Maybe add Threshold
    colname = Symbol(string("meanInten_", imgidx))
    appendCol!(newstructures, lbls, meanintensities, colname)

    ##Number of pixels
    lbls, npixels = getNPixels(attncat[imgidx], imgn) #### Also, threshold?
    colname = Symbol(string("nPixels_", imgidx))
    appendCol!(newstructures, lbls, npixels, colname)

    ## Count fos cells
    p = [(ptsw_coords[imgidx][j,1], ptsw_coords[imgidx][j,2]) for j in 1:size(ptsw_coords[imgidx],1)]
    amp = ptsw_amps[imgidx][:,1]
    lbls, foscounts, fosamps = getPosCounts(attncat[imgidx], p, amp)
    colname = Symbol(string("FosCounts_", imgidx))
    appendCol!(newstructures, lbls, foscounts, colname)
    colname = Symbol(string("FosAmplitudes_", imgidx))
    appendVecCol!(newstructures, lbls, fosamps, colname)
  end

  #### Append total fos count
  counts_dat = Matrix(newstructures[:, contains.(names(newstructures), "FosCounts")])
  counts_dat[ismissing.(counts_dat)] .= 0 #replace missing => 0
  counts_colsum = vec(sum(counts_dat, dims = 2))
  newstructures[!, "FosCounts_all"] = counts_colsum # append count column

  #### Append total Fos amplitudes
  amps_dat = Matrix(newstructures[:, contains.(names(newstructures), "FosAmplitudes")])
  amps_vec = [vcat(i...) for i in eachrow(amps_dat)]
  amps_vec = map(x -> x[.!ismissing.(x)], amps_vec) #amplitude per blob in each brain area
  newstructures[!, "FosAmplitudes_all"] = amps_vec  #append amplitude column
  return(newstructures)
end
