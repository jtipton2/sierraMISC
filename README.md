# sierraMISC
Miscellaneous Python tools for use with Sierra FEA software

## fitData.py
Requires Numpy, PyTetGen, and SKLearn modules.  Uses Delaunay triangularization to mesh a source temperature point cloud.  For a given set of target nodal coordinates, indexes which node falls inside which tetrahedron.  Interpolates temperature using Barycentric (areal) weighting.  If any target nodes fall outside the source temperature point cloud, uses K-Nearest Neighbors algoritm ot find the 3 nearest temperatures and average using a distance weighting.

## fitTemps.py
Requires the SEACAS exodus python module.  Uses fitData to interpolate temperature point clouds onto a mesh in ExodusII format.  Saves the temperature fields to a new ExodusII file as a function of time.
