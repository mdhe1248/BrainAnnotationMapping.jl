module BrainAnnotationMapping
using JSON

export retrieve, isdownstream, label_points, filter_points, count_cells
# Write your package code here.
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
  keepx = map(x -> x[1], pts) .> sz[2]
  keepy = map(x -> x[2], pts) .> sz[1]
  keepz = map(x -> x[3], pts) .> sz[3]
  keep = .!(keepx .| keepy .| keepz)
  return(pts[keep])
end

"""Cound cells in each brain region given `subbrain_fos_lbl`(each cell[row] with label id) and `target_id` (target brain id)"""
count_cells(subbrain_fos_lbl::Vector, target_id::Number) = sum(subbrain_fos_lbl .== target_id)  

end
