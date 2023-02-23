#======================================================================
# CREATE TEMPERATURE FIELDS FOR SIERRA DYNAMIC PULSE ANALYSIS
#======================================================================
import numpy as np
from exodus import exodus, copyTransfer


#
# Load Mesh with Beam Pulse Temperatures
#
mesh = exodus('PulseDT.g',mode='r',array_type='numpy')

numNodes = mesh.num_nodes()

TempDT = np.empty(numNodes, dtype=([('Index', int),('ID', int),('Val', float)]))

TempDT['ID'] = mesh.get_node_id_map()
TempDT['Index'] = range(0,numNodes,1)
TempDT['Val'] = mesh.get_node_variable_values('PulseDT',1)

# Sort by element ID.  Then I can use the sorted INDEX to 
# put the SS temperatures in the correct order
TempDT.sort(order='ID')
Arrangement = TempDT['Index'].copy()
TempDT.sort(order='Index')

mesh.close()


#
# Load Mesh with Steady-State Temperatures
#
mesh = exodus('LasagnaOpt_Thermal.e',mode='r',array_type='numpy')

#mesh.get_node_id_map()
# By inspection, it seems that nodal values in a results file
# are sorted sequentially.  The mesh files, however, have a
# different sorting.

numNodes = mesh.num_nodes()

TempSS = np.empty(numNodes, dtype=([('Index', int),('ID', int),('Val', float)]))

TempSS['ID'] = mesh.get_node_id_map()
TempSS['Index'] = Arrangement
TempSS['Val'] = mesh.get_node_variable_values('TEMP',1)
TempSS.sort(order='Index')

mesh.close()



#
# Save Temperature/Time Fields to New Exodus File
#
addGlobalVariables = []
addNodeVariables = ["NodalTempField"] 
addElementVariables = []
# copy exodus file and add variables to all nodes
exo = copyTransfer('LasagnaOpt_noShell.g','LasagnaOpt_Temps.g','ctype',addGlobalVariables,addNodeVariables, addElementVariables)

times = np.array([0.0, 
                  0.007, 
                  0.01, 
                  0.0100007, 
                  0.0140007]) 

temps = np.zeros([numNodes,len(times)])

temps[:,0] = 30.0  # HIP lock-in temperature
temps[:,1] = 30.0
temps[:,2] = TempSS['Val']
temps[:,3] = TempSS['Val']+TempDT['Val']
temps[:,4] = TempSS['Val']+TempDT['Val']

for ii in range(len(times)):
  exo.put_node_variable_values('NodalTempField',ii+1,temps[:,ii])
  exo.put_time(ii+1,float(times[ii]))

exo.close()
