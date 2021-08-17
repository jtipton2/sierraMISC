import numpy as np
import pytetgen as pytet
from sklearn import neighbors

def fitData(sourceCoord,sourceTemp,targetID,targetIndex,targetCoord):
    """Fit point cloud temperatures from a neutronics mesh onto FEA mesh nodes for
    thermostructural simulations.  Interpolation uses Delaunay triangularization
    with a Barycentric interpolation.  Any nodes external to the triangularization
    (due to coordinate rounding) are then extrapolated using a distance weighting of the
    3 nearest neighbors.
    
    Parameters:
    
    sourceCoord (float): x,y,z coordinates (m) of source point cloud
    sourceTemp (float): temperatures (deg C) of source point cloud
    targetID (int): Sierra mesh node_id_map
    targetIndex (int): Sierra mesh node index (1 to num_nodes)
    targetCoord (float): Sierra mesh node x,y,z coordinates (m)
    
    Returns:
    
    Numpy structured array of target Sierra nodes with columns:
    'Index' (int): echo of targetIndex
    'ID' (int): echo of targetID
    'x' (float): echo of targetCoord[:,0]
    'y' (float): echo of targetCoord[:,1]
    'z' (float): echo of targetCoord[:,2]
    'Temp' (float): nodal temperatures (deg-C)
       
    """
    
    #
    # Build Delaunay Tet Mesh from Source Point Cloud
    #
    tri = pytet.Delaunay(sourceCoord)
    # Any targets that fall outside the Delaunay cells return a
    # simplex value of -1.  Need to filter these and later
    # determine via some sort of extrapolation.
    tetsRaw = tri.find_simplex(targetCoord)
    outMask = tetsRaw > 0
    tets = tetsRaw[outMask]
    R = targetCoord[outMask]
    
    
    #
    # Calculate Barycentric Coordinates | Local Weights for Each Target
    # Coordinate Inside Source Tetrahedron
    #
    """
    Background:  Interpolation of point clouds in Abaqus is painfully slow.  SciPy has a way to do this
    using LinearNDInterpolator, however, it is also painfully slow.  I tried using <stackoverflow.com/a/56628900>
    and discovered that it is the Delaunay triangularization that is taking a long time.  Amazingly, a set of
    points that still had not completed in over a day was handled in seconds using pyTetGen.  The issue is
    that SciPy works in python while pyTetGen interfaces directly with C libraries.  What follows is then my
    reverse engineering of barycentric interpolation
    
    Basic idea and formula in:
    https://codeplea.com/triangular-interpolation
    https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    
    Great code example:  https://codereview.stackexchange.com/a/41089
    But uses scipy.spatial.Delaunay.transform function
    I had to work from the Wikipedia example to get a workaround.
    If I use scipy.spatial.Delaunay(points), I get the same answer as
    in the link.
    
    How to piecewise invert a collection of matrices:
    https://stackoverflow.com/questions/17924411/vectorized-partial-inverse-of-an-nmm-tensor-with-numpy
    
    Using np.newaxis to add a dimension to an array:
    https://stackoverflow.com/questions/26333005/numpy-subtract-every-row-of-matrix-by-vector/26333184
    
    How to transpose individual matrices nested inside an array:
    https://jameshensman.wordpress.com/2010/06/14/multiple-matrix-multiplication-in-numpy/
    """
    T = sourceCoord[tri.simplices[tets,:3]] - sourceCoord[tri.simplices[tets,3]][:, np.newaxis]
    Ttrans = np.transpose(T,(0,2,1))
    Tinv = np.linalg.inv(Ttrans)
    
    RminusR4 = R - sourceCoord[tri.simplices[tets,3]]
    
    b = np.einsum('ijk,ik->ij', Tinv, RminusR4)
    
    bcoords = np.c_[b, 1 - b.sum(axis=1)]
    
    
    #
    # Calculate Interpolated Target Temperatures with Coordinates and Node ID
    #
    outTemp = np.multiply(sourceTemp[tri.simplices[tets]],bcoords).sum(axis=1)
    outCoord = targetCoord[outMask]
    outID = targetID[outMask]
    outIndex = targetIndex[outMask]
    
    
    #
    # Unmatched Nodes - KNN Distance Weighted Averaging
    #
    # np.savetxt('unmatched.txt',targetCoord[outMask!=True],delimiter=',')
    #
    # (Visual inspection of this node set in ParaView showed it to be the
    #  outer surface nodes.  There's a rounding problem on node coords
    #  between Abaqus, Atilla, and Sierra.  Use a distance weighted average
    #  of the 3 closest neighbors for these elements.)
    #
    # https://stackoverflow.com/questions/48312205/find-the-k-nearest-neighbours-of-a-point-in-3d-space-with-python-numpy
    #
    outerCoord = targetCoord[outMask!=True]
    outerID = targetID[outMask!=True]
    outerIndex = targetIndex[outMask!=True]
    knn = neighbors.KNeighborsRegressor(3, weights='distance')
    outerTemp = knn.fit(sourceCoord,sourceTemp).predict(outerCoord)
    
    
    #
    # Join the results from the interpolation and extrapolations
    # Join everything into 1 array
    # Sort the array by Sierra node index
    # 
    finalIndex = np.hstack((outIndex,outerIndex))
    finalID = np.hstack((outID,outerID))
    finalCoord = np.vstack((outCoord,outerCoord))
    finalTemp = np.hstack((outTemp,outerTemp))
    
    #final = np.column_stack([finalIndex, finalID, finalCoord, finalTemp])
    #final[final[:,0].argsort()]
    
    final = np.empty(len(finalIndex), dtype=([('Index', int), 
                                              ('ID', int), 
                                              ('x', float), 
                                              ('y', float), 
                                              ('z', float), 
                                              ('Temp', float)]))
    
    final['Index'] = finalIndex
    final['ID'] = finalID
    final['x'] = finalCoord[:,0]
    final['y'] = finalCoord[:,1]
    final['z'] = finalCoord[:,2]
    final['Temp'] = finalTemp
    
    final.sort(order='Index')
    
    return final

