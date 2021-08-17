from fitData import *
from exodus import exodus, copyTransfer

#======================================================================
# FIT POINT CLOUD TEMPERATURE SETS TO SIERRA MESH
#======================================================================

#
# Load Target Mesh
#
mesh = exodus('Mesh_Target.g',mode='r',array_type='numpy')

targetID = np.array(mesh.get_node_id_map())
targetIndex = np.arange(1,len(targetID)+1)
targetCoord = np.transpose(mesh.get_coords())

mesh.get_elem_blk_names()
mesh.get_node_set_ids()
mesh.get_node_set_names()

mesh.get_node_set_name(3)
tantalumIndex = mesh.get_node_set_nodes(3)
tantalumID = tantalumIndex - 1
tantalumCoord = targetCoord[tantalumID]

mesh.get_node_set_name(6)
tungstenIndex = mesh.get_node_set_nodes(6)
tungstenID = tungstenIndex - 1
tungstenCoord = targetCoord[tungstenID]

mesh.close()

#
# Load Source Mesh and Temperatures and Fit the Data
#

source = np.loadtxt(open('DTTa','rb'), delimiter=',', skiprows=0)
sourceCoord = source[:,:3]
sourceTemp = source[:,3]
TaPulseA= fitData(sourceCoord,sourceTemp,tantalumID,tantalumIndex,tantalumCoord)

source = np.loadtxt(open('DTTa121','rb'), delimiter=',', skiprows=0)
sourceCoord = source[:,:3]
sourceTemp = source[:,3]
TaPulseB = fitData(sourceCoord,sourceTemp,tantalumID,tantalumIndex,tantalumCoord)

source = np.loadtxt(open('Tantalum_Temperature_deg_C.csv','rb'), delimiter=',', skiprows=1)
sourceCoord = source[:,:3]
sourceTemp = source[:,3]
TaSS = fitData(sourceCoord,sourceTemp,tantalumID,tantalumIndex,tantalumCoord)

source = np.loadtxt(open('DTW','rb'), delimiter=',', skiprows=0)
sourceCoord = source[:,:3]
sourceTemp = source[:,3]
WPulseA= fitData(sourceCoord,sourceTemp,tungstenID,tungstenIndex,tungstenCoord)

source = np.loadtxt(open('DTW121','rb'), delimiter=',', skiprows=0)
sourceCoord = source[:,:3]
sourceTemp = source[:,3]
WPulseB = fitData(sourceCoord,sourceTemp,tungstenID,tungstenIndex,tungstenCoord)

source = np.loadtxt(open('Tungsten_Temperature_deg_C.csv','rb'), delimiter=',', skiprows=1)
sourceCoord = source[:,:3]
sourceTemp = source[:,3]
WSS = fitData(sourceCoord,sourceTemp,tungstenID,tungstenIndex,tungstenCoord)

finalPulseA = np.hstack((TaPulseA,WPulseA))
finalPulseB = np.hstack((TaPulseB,WPulseB))
finalSS = np.hstack((TaSS,WSS))


#
# Save Temperature/Time Fields to New Exodus File
#
addGlobalVariables = []
addNodeVariables = ["NodalTempField"] 
addElementVariables = []
# copy exodus file and add variables to all nodes
exo = copyTransfer('Mesh_Target.g','Temps_Target.g','ctype',addGlobalVariables,addNodeVariables, addElementVariables)

times = np.array([0.0, 
                  0.0024, 
                  0.0048, 
                  0.0048007, 
                  0.0088007, 
                  0.0088014, 
                  0.0128014])

temps = np.zeros([len(targetID),len(times)])

temps[:,0] = 450.0
temps[:,1] = 30.0
temps[:,2] = finalSS['Temp']
temps[:,3] = finalSS['Temp'] + finalPulseA['Temp']
temps[:,4] = finalSS['Temp'] + finalPulseA['Temp']
temps[:,5] = finalSS['Temp'] + finalPulseB['Temp']
temps[:,6] = finalSS['Temp'] + finalPulseB['Temp']

for ii in range(len(times)):
  exo.put_node_variable_values('NodalTempField',ii+1,temps[:,ii])
  exo.put_time(ii+1,float(times[ii]))

exo.close()

