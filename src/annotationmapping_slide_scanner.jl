"""
#### Bouondary detection
"""
function annotation_boundary(annotation_img)
  parent_img = copy(parent(annotation_img))
  output_img .= zeros(size(parent_img))
  for i in Base.axes(output_img,3)
    im1 = diff(parent_img[:,:,i], dims = 1)
    output_img[2:end, :, i] .+= im1
    im2 = diff(parent_img[:,:,i], dims = 2)
    output_img[:, 2:end, i] .+= im2
  end
  output_img[output_img .!= 0] .= 1
  return(output_img)
end
