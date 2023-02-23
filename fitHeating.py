from fitData import *
from exodus import exodus, copyTransfer
import pandas as pd

#==============================================================================
# FIT POINT CLOUD NEUTRONICS SETS TO SIERRA MESH
#==============================================================================
#
# Change Log:
#
# 2023-02-03
#   --> Expanded script to handle cases with shared nodes 
#   --> Fixed node sets to order by Index instead of ID
#   --> Incorporated Pandas DataFrames
#   --> Used DataFrames to take the largest values on shared nodes and
#       merge on node Index values
#
# 2023-02-17
#   --> Created a NodeSet dictionary to search on NodeSet names instead
#       of ID numbers
#


#
# Target Parameters -----------------------------------------------------------
#
PulseFreq = 15.0 / 21.0  # [s]
from MatlProps_Pulse import WDensity, In718Density, CuDensity


#
# Load Target Mesh ------------------------------------------------------------
#
mesh = exodus('LasagnaOpt_noShell.g',mode='r',array_type='numpy')

nodesetIDs = mesh.get_node_set_ids()
nodesetNames = mesh.get_node_set_names()
nodesetDict = {nodesetNames[i]: nodesetIDs[i] for i in range(len(nodesetIDs))}

targetID = np.array(mesh.get_node_id_map())
targetIndex = np.arange(1,len(targetID)+1)  # use this at the end to get the nodes back in the right order before saving
targetCoord = np.transpose(mesh.get_coords())
target = pd.DataFrame({'ID':targetID,'Index':targetIndex,'CoordX':targetCoord[:,0],'CoordY':targetCoord[:,1],'CoordZ':targetCoord[:,2]})

cladIndex = mesh.get_node_set_nodes(nodesetDict["NS_CLAD"])
clad = pd.DataFrame({'Index':cladIndex})
clad = pd.merge(clad,target)
#clad.to_csv('NS_CLAD.csv')

blockIndex = mesh.get_node_set_nodes(nodesetDict["NS_BLOCK"])
block = pd.DataFrame({'Index':blockIndex})
block = pd.merge(block,target)
#block.to_csv('NS_BLOCK.csv')

shroudIndex = mesh.get_node_set_nodes(nodesetDict["NS_SHROUD"])
shroud = pd.DataFrame({'Index':shroudIndex})
shroud = pd.merge(shroud,target)
#shroud.to_csv('NS_SHROUD.csv')

mesh.close()


#
# Load Source Mesh and Energy Deposition and Fit the Data ---------------------
#

# CLADDING - NEUTRONICS MODEL
#sourceL = np.loadtxt(open('Inputs/target_3blk_Cu_L.csv','rb'), delimiter=',', skiprows=0)
#sourceC = np.loadtxt(open('Inputs/target_3blk_Cu_C.csv','rb'), delimiter=',', skiprows=0)
#sourceR = np.loadtxt(open('Inputs/target_3blk_Cu_R.csv','rb'), delimiter=',', skiprows=0)
#sourceCoord = sourceC[:,:3]*0.01                            # [m]
#sourceVal = (sourceL[:,3]+sourceC[:,3]+sourceR[:,3])*1.0e6*PulseFreq  # [W/m^3]
# CLADDING - MATCAD SIMPLIFIED MODEL
source = np.loadtxt(open('Inputs/60cm2_bodyround.csv','rb'), delimiter=',', skiprows=1)
sourceCoord = source[:,:3] # [m]
sourceVal = source[:,3]/WDensity*CuDensity    # [W/m^3]
sourceVal[sourceVal < 0.0] = 0.0
# CLADDING - DATA FIT 
cladSS = fitData(sourceCoord,sourceVal,clad['ID'].values,clad['Index'].values,clad[['CoordX','CoordY','CoordZ']].values)

# BLOCK - NEUTRONICS MODEL
#sourceL = np.loadtxt(open('Inputs/target_3blk_W_L.csv','rb'), delimiter=',', skiprows=0)
#sourceC = np.loadtxt(open('Inputs/target_3blk_W_C.csv','rb'), delimiter=',', skiprows=0)
#sourceR = np.loadtxt(open('Inputs/target_3blk_W_R.csv','rb'), delimiter=',', skiprows=0)
#sourceCoord = sourceC[:,:3]*0.01 # [m]      
#sourceVal = (sourceL[:,3]+sourceC[:,3]+sourceR[:,3])*1.0e6*PulseFreq  # [W/m^3]
# BLOCK - MATCAD SIMPLIFIED MODEL
source = np.loadtxt(open('Inputs/60cm2_round.csv','rb'), delimiter=',', skiprows=1)
sourceCoord = source[:,:3] # [m]
sourceVal = source[:,3]    # [W/m^3]
sourceVal[sourceVal < 0.0] = 0.0
# BLOCK - DATA FIT
blockSS = fitData(sourceCoord,sourceVal,block['ID'].values,block['Index'].values,block[['CoordX','CoordY','CoordZ']].values)

# SHROUD - NEUTRONICS MODEL
#sourceL = np.loadtxt(open('Inputs/target_3blk_Inconel_L.csv','rb'), delimiter=',', skiprows=0)
#sourceC = np.loadtxt(open('Inputs/target_3blk_Inconel_C.csv','rb'), delimiter=',', skiprows=0)
#sourceR = np.loadtxt(open('Inputs/target_3blk_Inconel_R.csv','rb'), delimiter=',', skiprows=0)
#sourceCoord = sourceC[:,:3]*0.01                         # [m]
#sourceVal = (sourceL[:,3]+sourceC[:,3]+sourceR[:,3])*1.0e6*PulseFreq  # [W/m^3]
# SHROUD - MATCAD SIMPLIFIED MODEL
source = np.loadtxt(open('Inputs/60cm2_bodyround.csv','rb'), delimiter=',', skiprows=1)
sourceCoord = source[:,:3] # [m]
sourceVal = source[:,3]/WDensity*In718Density    # [W/m^3]
sourceVal[sourceVal < 0.0] = 0.0
# SHROUD - DATA FIT
shroudSS = fitData(sourceCoord,sourceVal,shroud['ID'].values,shroud['Index'].values,shroud[['CoordX','CoordY','CoordZ']].values)


#
# Join and Order Data ---------------------------------------------------------
#
finalSS = np.hstack((blockSS,cladSS,shroudSS))
finalSS = pd.DataFrame(finalSS)
# Sort, take the max of duplicate node entries
finalSS = finalSS.groupby('Index').max().reset_index(level=0)


#
# Save Heat Generation Values to New Exodus File ------------------------------
#
addGlobalVariables = []
addNodeVariables = ["VolHeatGen"] 
addElementVariables = []

# copy exodus file and add values to all nodes
exo1 = copyTransfer('LasagnaOpt_noShell.g','VolHeatGen.g','ctype',addGlobalVariables,addNodeVariables, addElementVariables)
exo1.put_node_variable_values('VolHeatGen',1,finalSS['Val'].values)
exo1.put_time(1,1.0)
exo1.close()

