# sierraMISC
Miscellaneous Python tools for use with Sierra FEA software

## fitData.py
Requires Numpy, PyTetGen, and SKLearn modules.  Uses Delaunay triangularization to mesh a source temperature point cloud.  For a given set of target nodal coordinates, indexes which node falls inside which tetrahedron.  Interpolates temperature using Barycentric (areal) weighting.  If any target nodes fall outside the source temperature point cloud, uses K-Nearest Neighbors algoritm ot find the 3 nearest temperatures and average using a distance weighting.

## fitHeating.py & fitPulse.py
Requires the SEACAS exodus python module.  Uses fitData to interpolate neutronics point clouds onto a mesh in ExodusII format.

## makeTemps.py
Requires the SEACAS exodus python module.  Opens a heat transfer solution and opens a neutronics proton pulse temperature rise field.  Creates a new ExodusII file that contains combined temperature fields at defeined timestamps.  (This file is then used to drive a Sierra Explicit Dynamic simulation driven by temperature field-induced thermal expansion.)

## sierraExport.py
Requires the SEACAS exodus python module.  Gathers the elemental stress tensors at all available time steps and saves the results in binary Numpy format.

## sierra2ODB.py
Requires an empty Abaqus ODB file with the mesh loaded.  The suggested procedure to do this is to import the Sierra analysis mesh into CUBIT and then export as Abaqus INP format.  This creates a full Abaqus input deck including fake material definition and fake solution step.  Next, the Abaqus “datacheck” command can be used to convert the input file into an ODB results file.  The results file will be empty except for the mesh.

```bash
abaqus datacheck job=LasagnaOpt_Shroud input=LasagnaOpt_Shroud.inp interactive
```

The sierra2ODB script uses Abaqus python to write the Sierra stress tensor results into the empty ODB file.  The script first checks if a step named “Sierra” exists in the ODB file.  If it does not, the step is created.  Then a count-controlled loop writes time frames and stress results using the previously exported Numpy arrays.  This can be very memory intensive, depending on the analysis size, and may crash before completion.  If that is the case, the script should be rerun.  If the step “Sierra” already exists in the ODB, the code script will find the last written time frame and then continue writing data from that point.
