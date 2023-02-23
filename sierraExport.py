"""
Post-processes Sierra output file to save elemental stress tensors
at all available time steps.  Saves results in binary numpy format.

Built by Joseph B. Tipton, Jr. from original concept by Lianshan Lin (ORNL).

Change Log:

2021-10-11
  --> Original issue

2023-02-20
  --> Added support for elements with four integration points
  --> This amount of additional data then exceeded memory limits, fix
      was to then write the stress array to a npy file at each time
      increment

"""

import numpy as np
from exodus import exodus, copyTransfer
import exodusCalcs


# set element type
numIP = 4


# read stress data from sierra output .e file
filename = 'LasagnaOpt_Dynamic_Shroud_Pulses'
ei = exodus(filename+'.e',mode='r',array_type='numpy')


# exodus block dictionary
blkIDs = np.array(ei.get_elem_blk_ids())   # block integer ids
blkNames = np.array(ei.get_elem_blk_names()) # block names
blkDict = {blkNames[i]: blkIDs[i] for i in range(len(blkIDs))}


# element IDs
eleIDs = np.transpose(ei.get_elem_num_map())


# Lists the times saved in seconds
# The exodus time_step index starts at 1 (not zero)
timeArray = ei.get_times()
numTimes = len(timeArray)

# Preallocate master array
numElem = ei.num_elems()

# assumes output file has only 1 block and selects that ID
ii = blkIDs[0][0]


if numIP == 1:
  # Element has 1 integration point (e.g. C3D8R or HEX)
  for jj in range(numTimes):
    # read the stress components into array
    #                                (elem_blk_id, evar_name, time_step)
    evar_xx = ei.get_element_variable_values(ii, 'Stress_xx', jj+1)
    evar_yy = ei.get_element_variable_values(ii, 'Stress_yy', jj+1)
    evar_zz = ei.get_element_variable_values(ii, 'Stress_zz', jj+1)
    evar_xy = ei.get_element_variable_values(ii, 'Stress_xy', jj+1)
    evar_zx = ei.get_element_variable_values(ii, 'Stress_zx', jj+1)
    evar_yz = ei.get_element_variable_values(ii, 'Stress_yz', jj+1)
    # construct stress array
    eleStress=(np.array([evar_xx, evar_yy, evar_zz, evar_xy, evar_zx, evar_yz])).T
    # save stress array at given time step increment
    np.save(filename+'_'+str(jj).zfill(3)+'.npy', eleStress, allow_pickle=True)
elif numIP == 4:
  # Element has 4 integration points (e.g. C3D10M or TETRA10)
  for jj in range(numTimes):
    # read the stress components into array
    #                                (elem_blk_id, evar_name, time_step)
    evar_xx_1 = ei.get_element_variable_values(ii, 'Stress_xx_1', jj+1)
    evar_yy_1 = ei.get_element_variable_values(ii, 'Stress_yy_1', jj+1)
    evar_zz_1 = ei.get_element_variable_values(ii, 'Stress_zz_1', jj+1)
    evar_xy_1 = ei.get_element_variable_values(ii, 'Stress_xy_1', jj+1)
    evar_zx_1 = ei.get_element_variable_values(ii, 'Stress_zx_1', jj+1)
    evar_yz_1 = ei.get_element_variable_values(ii, 'Stress_yz_1', jj+1)
    evar_xx_2 = ei.get_element_variable_values(ii, 'Stress_xx_2', jj+1)
    evar_yy_2 = ei.get_element_variable_values(ii, 'Stress_yy_2', jj+1)
    evar_zz_2 = ei.get_element_variable_values(ii, 'Stress_zz_2', jj+1)
    evar_xy_2 = ei.get_element_variable_values(ii, 'Stress_xy_2', jj+1)
    evar_zx_2 = ei.get_element_variable_values(ii, 'Stress_zx_2', jj+1)
    evar_yz_2 = ei.get_element_variable_values(ii, 'Stress_yz_2', jj+1)
    evar_xx_3 = ei.get_element_variable_values(ii, 'Stress_xx_3', jj+1)
    evar_yy_3 = ei.get_element_variable_values(ii, 'Stress_yy_3', jj+1)
    evar_zz_3 = ei.get_element_variable_values(ii, 'Stress_zz_3', jj+1)
    evar_xy_3 = ei.get_element_variable_values(ii, 'Stress_xy_3', jj+1)
    evar_zx_3 = ei.get_element_variable_values(ii, 'Stress_zx_3', jj+1)
    evar_yz_3 = ei.get_element_variable_values(ii, 'Stress_yz_3', jj+1)
    evar_xx_4 = ei.get_element_variable_values(ii, 'Stress_xx_4', jj+1)
    evar_yy_4 = ei.get_element_variable_values(ii, 'Stress_yy_4', jj+1)
    evar_zz_4 = ei.get_element_variable_values(ii, 'Stress_zz_4', jj+1)
    evar_xy_4 = ei.get_element_variable_values(ii, 'Stress_xy_4', jj+1)
    evar_zx_4 = ei.get_element_variable_values(ii, 'Stress_zx_4', jj+1)
    evar_yz_4 = ei.get_element_variable_values(ii, 'Stress_yz_4', jj+1)
    # construct stress array
    SIP1 = (np.array([evar_xx_1, evar_yy_1, evar_zz_1, evar_xy_1, evar_zx_1, evar_yz_1])).T
    SIP2 = (np.array([evar_xx_2, evar_yy_2, evar_zz_2, evar_xy_2, evar_zx_2, evar_yz_2])).T                                                  
    SIP3 = (np.array([evar_xx_3, evar_yy_3, evar_zz_3, evar_xy_3, evar_zx_3, evar_yz_3])).T                         
    SIP4 = (np.array([evar_xx_4, evar_yy_4, evar_zz_4, evar_xy_4, evar_zx_4, evar_yz_4])).T                         
    # stackoverflow.com/questions/5347065/interweaving-two-numpy-arrays
    eleStress = np.empty((4*SIP1.shape[0],SIP1.shape[1]), dtype=SIP1.dtype)
    eleStress[0::4] = SIP1
    eleStress[1::4] = SIP2
    eleStress[2::4] = SIP3
    eleStress[3::4] = SIP4
    # save stress array at given time step increment
    np.save(filename+'_'+str(jj).zfill(3)+'.npy', eleStress, allow_pickle=True)

ei.close()

# save data
#np.savez(filename, eleIDs=eleIDs, blkDict=blkDict, timeArray=timeArray, outData=outData)
with open(filename+'.npy', 'wb') as f:
  np.save(f, eleIDs)
  np.save(f, blkDict)
  np.save(f, timeArray)



