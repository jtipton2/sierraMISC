"""
Post-processes Sierra output file to save elemental stress tensors
at all available time steps.  Saves results in binary numpy format.

Built by Joseph B. Tipton, Jr. from original concept by Lianshan Lin (ORNL).

"""

import numpy as np
from exodus import exodus, copyTransfer
import exodusCalcs


#read stress data from sierra output .e file
filename = 'Sierra_STS_ST_Hex_Ta_Pulses'
ei = exodus(filename+'.e',mode='r',array_type='numpy')


# exodus block dictionary
blkIDs = np.array(ei.get_elem_blk_ids())   # block integer ids
blkNames = np.array(ei.get_elem_blk_names()) # block names
blkDict = {blkNames[i]: blkIDs[i] for i in range(len(blkIDs))}


#element IDs
eleIDs = np.transpose(ei.get_elem_num_map())


# Lists the times saved in seconds
# The exodus time_step index starts at 1 (not zero)
timeArray = ei.get_times()
numTimes = len(timeArray)


# Preallocate master array
numElem = ei.num_elems()
outData = np.empty((numTimes,numElem,6))


for jj in range(numTimes):
  ii = 0
  #read the stress components into array
  #                                  (elem_blk_id, evar_name, time_step)
  evar_xx = ei.get_element_variable_values(ii+1, 'Stress_xx', jj+1)
  evar_yy = ei.get_element_variable_values(ii+1, 'Stress_yy', jj+1)
  evar_zz = ei.get_element_variable_values(ii+1, 'Stress_zz', jj+1)
  evar_xy = ei.get_element_variable_values(ii+1, 'Stress_xy', jj+1)
  evar_zx = ei.get_element_variable_values(ii+1, 'Stress_zx', jj+1)
  evar_yz = ei.get_element_variable_values(ii+1, 'Stress_yz', jj+1)
  for ii in range(len(blkIDs)-1):
    evar_xx = np.hstack((evar_xx,ei.get_element_variable_values(ii+2, 'Stress_xx', jj+1)))
    evar_yy = np.hstack((evar_yy,ei.get_element_variable_values(ii+2, 'Stress_yy', jj+1)))
    evar_zz = np.hstack((evar_zz,ei.get_element_variable_values(ii+2, 'Stress_zz', jj+1)))
    evar_xy = np.hstack((evar_xy,ei.get_element_variable_values(ii+2, 'Stress_xy', jj+1)))
    evar_zx = np.hstack((evar_zx,ei.get_element_variable_values(ii+2, 'Stress_zx', jj+1)))
    evar_yz = np.hstack((evar_yz,ei.get_element_variable_values(ii+2, 'Stress_yz', jj+1)))

  # construct stress array
  evar_all=(np.array([evar_xx, evar_yy, evar_zz, evar_xy, evar_zx, evar_yz])).T
  
  # add to output array
  outData[jj] = evar_all 


ei.close()

# save data
#np.savez(filename, eleIDs=eleIDs, blkDict=blkDict, timeArray=timeArray, outData=outData)
with open(filename+'.npy', 'wb') as f:
  np.save(f, eleIDs)
  np.save(f, blkDict)
  np.save(f, timeArray)
  np.save(f, outData)


