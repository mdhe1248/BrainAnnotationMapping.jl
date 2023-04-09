using JSON, DataFrames, CSV

""" warpper function"""
function brainmapping(annotationImg, annotation, pts_filtered)
  subbrain_ids = Int.(unique(annotationImg))
  subbrain_labels = map(x -> retrieve(annotation[1], "name", "id", x), subbrain_ids)
  subbrain_fos_lbl = label_points(annotationImg, pts_filtered) # label each point
  fosncells = map(x -> count_cells(subbrain_fos_lbl, x), subbrain_ids) # Count cells in each subbrain
  parent_ids = map(x -> retrieve(annotation[1],"parent_structure_id", "id", x), subbrain_ids)
  return(subbrain_ids, subbrain_labels, subbrain_fos_lbl, fosncells, parent_ids)
end

"""
Find the Dict value from `key_of_interest` in a nested Dictionary.
A nested dictionary might be from JSON file.
"""
function retrieve(dict, key_of_interest, key_input, value_input)
  for (key, value) in dict
    if key == key_input && value == value_input
      if haskey(dict, key_of_interest)
        return(dict[key_of_interest])
      else
        return(nothing)
      end
    end
    if value isa Vector
      for dict in value
        label = retrieve(dict, key_of_interest, key_input, value_input)
        if !isnothing(label)
          return(label)
        end
      end
    end
  end
end

""" check if 'input_val_downstream' is the downstream of 'input_value'"""
function isdownstream(dict, input_key, input_value, input_val_downstream)
  parent_val= retrieve(dict, "parent_structure_id", input_key, input_val_downstream)
  while !isnothing(parent_val)
    if input_value .== parent_val
      return(true)
    end
    parent_val= retrieve(dict, "parent_structure_id", input_key, parent_val)
  end
  return(false)
end

""" find the subbrain labeling of blobs within target brain region
`annotationImg` is brain with labeling. It is from ClearMap2 reference (e.g.`ABA_25um_annotation.tif`).
"""
#`blobs` are obtained from blob_LoG; possibly contains the coordinate of c-fos.
function label_points(annotationImg, points::Vector)
#  fos_coords = map(x -> x.location, blobs)
  keep_lbl = zeros(Int, length(points))
  for (i, fos_coord) in enumerate(points)
    lbl = annotationImg[fos_coord[2], fos_coord[1], fos_coord[3]]
    keep_lbl[i] = lbl
  end
  keep_lbl
end

""" filter points outside annotation image"""
function filter_points(img, pts)
  sz = size(img) 
  keepx = map(x -> x[1], pts) .> sz[2] .|| map(x -> x[1], pts) .< 1
  keepy = map(x -> x[2], pts) .> sz[1] .|| map(x -> x[2], pts) .< 1
  keepz = map(x -> x[3], pts) .> sz[3] .|| map(x -> x[3], pts) .< 1
  keep = .!(keepx .| keepy .| keepz)
  return(pts[keep], keep)
end

"""Cound cells in each brain region given `subbrain_fos_lbl`(each cell[row] with label id) and `target_id` (target brain id)"""
count_cells(subbrain_fos_lbl::Vector, target_id::Number) = sum(subbrain_fos_lbl .== target_id)  

function load_brainmap_json(jsonfn)
  d = []
  open(jsonfn) do io
    push!(d, JSON.parse(io))
  end
  return(d)
end

function saveData(results_dir, subbrain_ids, subbrain_labels, fosncells)
  dtf = DataFrame(subbrain_ids = subbrain_ids, subbrain_labels = subbrain_labels, fos_ncells = fosncells)
  dtf.subbrain_labels = replace(dtf.subbrain_labels, nothing => missing)
  isdir(results_dir) ? nothing : mkdir(results_dir)
  CSV.write(results_dir*"cfos_counts.csv", dtf, delim = ';')
end

""" voxelize viualization
`pos` a vector of tuple in 3d
`amp` amplitude at the corresponding position.
`r` voxelize radius
`sz` image size
"""
function voxelize_roi(sz, pos, amps, r; gaussian = true)
  img1 = zeros(sz)
  if gaussian == true
    k = Kernel.gaussian((r,r,r))
  else
    k = ones(r,r,r)
  end
  for i in eachindex(pos)
     p1, amp = pos[i], amps[i]
    _voxelize_roi!(img1, p1, amp, r, k)
  end
  return(img1)
end

function _voxelize_roi!(img1, p1, amp, r, k)
  fi = max(CartesianIndex(p1.-r), CartesianIndex(1,1,1))
  li = min(CartesianIndex(p1.+r), CartesianIndex(size(img1)))
  for (i, ci) in enumerate(fi:li)
    if sqrt((ci.I[1]-p1[1])^2+(ci.I[2]-p1[2])^2+(ci.I[3]-p1[3])^2) <= r
      amp = amp
      img1[ci] = img1[ci]+amp*k[i]
    end
  end
  return(img1)
end

"""
filter blobs by position
`range_limit` is a vector of UnitRange.
e.g.) [1:300, 1:300, 1:200]
"""
function blob_pos_filter(blobs, range_limit)
  blobf = Vector{eltype(blobs)}()
  for blob in blobs
    pos = blob.location.I
    if all(pos .∈  range_limit)
      push!(blobf, blob)
    end
  end
  blobf
end 

"""
Crop blobs by position. blob locations will be re-positioned according to `coords`
e.g.) coords = [101:200, 201:300, 301:400]
if the original blob.location is (111, 201, 301), the cropped position is (11, 1, 1)
"""
function crop_blobs(blobs, coords)
  blobf = blob_pos_filter(blobs, coords)
  blobf1 = similar(blobf)
  for (i, blob) in enumerate(blobf) #repositioning
    loc = CartesianIndex(blob.location[1]-first(coords[1])+1, blob.location[2]-first(coords[2])+1, blob.location[3]-first(coords[3])+1)
    amp = blob.amplitude
    σ = blob.σ
    blobf1[i] = BlobLoG(loc, σ, amp)
  end
  blobf1
end
