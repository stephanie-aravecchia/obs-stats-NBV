This git repository contains the code related to the IROS 2023 submitted paper:
"Next-Best-View selection from observation viewpoint statistics", Stephanie Aravecchia, Antoine Richard, Marianne Clausel, Cedric Pradalier, 2023.

The main packages are:

mesh_slicer: to slice a mesh into images
slice_octomap: to slice an octomap octree into images
comparator: to compare the dataset reconstruction and ground-truth obtained from the two previous slicers

obs_stat: the ROS node to compute viewpoint statistics

explo_quality_aware: the ROS node to explore with NBV policies based on viewpoint statistics

The simulated environments can be found here:
https://github.com/stephanie-aravecchia/unstructured-env-simulator


The comparaison is implemented on a specific folder hierarchy, please follow the guidelines below to run them as they are provided
in base_folder/
ground_truth_stl/stl/xp_id  the stl and off files (the slicer input is a .off file)
ground_truth_stl/img/xp_id  the sliced images from mesh_slicer are to be put here

3d_grid/ot/xp_id  the octomap xp_id.ot
3d_grid/img/xp_id slice_octomap will output the images here

3d_grid/res/obs/xp_id the comparator wil read in the two img folders, and will output the result here