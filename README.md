# BrainAnnotationMapping.jl
Example script. See ElastixRegistration.jl before running this.
```jl
###### 3 The last: Load annotation file and count the number of c-fos cells
## Initialize variable
jsonfn = "/home/donghoon/usr/ClearMap2/ClearMap/Resources/Atlas/ABA_annotation.json"
results_dir = "counts/"
annotation_map = load_brainmap_json(jsonfn) #Load json file (brain annotation map)
subbrain_ids, subbrain_labels, subbrain_fos_lbl, fosncells, parent_ids = brainmapping(annotationImg, annotation_map, pts_filtered) #point-to-brain mapping
saveData(results_dir, subbrain_ids, subbrain_labels, fosncells) #save counting data
```
