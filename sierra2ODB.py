"""
sierra2ODB.py

Depends upon sierraExport.py script which saves Sierra simulation stress
tensors in numpy binary format.  This program loads this data along with
an Abaqus ODB that already contains the corresponding mesh.  It then writes
the data along with time steps to the ODB.

Built by Joseph B. Tipton, Jr. from original concept by Lianshan Lin (ORNL).

"""

from odbAccess import *
from abaqusConstants import *
import numpy as np


#
# LOAD DATA FROM SIERRA
#
filename = 'Sierra_STS_ST_Hex_Ta_Pulses'

with open(filename+'.npy', 'rb') as f:
  eleIDs = np.load(f, allow_pickle=True).tolist()
  blkDict = np.load(f, allow_pickle=True)
  timeArray = np.load(f, allow_pickle=True)
  eleStress = np.load(f, allow_pickle=True)


#
# ACCESS ODB and ADD STEP with FRAMES
#
odb = openOdb(filename+'.odb',readOnly=False)
allInstances = odb.rootAssembly.instances.keys()
odbInstance = odb.rootAssembly.instances[allInstances[-1]]

if odb.steps.keys()[-1] != 'Sierra':
  #
  # The data transfer has not yet started,
  # so, make the step and start the transfer
  #
  step1 = odb.Step(name='Sierra',
                   description='Sierra pulse results',
                   domain=TIME,
                   timePeriod=timeArray[-1])
  
  for ii in range(len(timeArray)):
    print('Starting frame ',ii)
    Sframe = step1.Frame(incrementNumber=ii, 
                         frameValue=timeArray[ii])
    Sfield = Sframe.FieldOutput(name='S',
                                description='Stress', 
                                type=TENSOR_3D_FULL, 
                                componentLabels=('S11', 'S22', 'S33', 'S12', 'S13', 'S23'), 
                                validInvariants=(MISES,TRESCA,MAX_PRINCIPAL))
    Sfield.addData(position=INTEGRATION_POINT, 
                   instance=odbInstance, 
                   labels=eleIDs, 
                   data=eleStress[ii].tolist())
    # save every 1000th time step
    if (ii % 1000 == 0):
      odb.save()
      # gc.collect(generation=2) doesn't really help
else:
  #
  # The data transfer previously started but might
  # have ended prematurely due to a memory error.
  # Pickup where we left off...
  #
  step1 = odb.steps['Sierra']
  for ii in range(len(step1.frames),len(timeArray)):
    print('Starting frame ',ii)
    Sframe = step1.Frame(incrementNumber=ii, 
                         frameValue=timeArray[ii])
    Sfield = Sframe.FieldOutput(name='S',
                                description='Stress', 
                                type=TENSOR_3D_FULL, 
                                componentLabels=('S11', 'S22', 'S33', 'S12', 'S13', 'S23'), 
                                validInvariants=(MISES,TRESCA,MAX_PRINCIPAL))
    Sfield.addData(position=INTEGRATION_POINT, 
                   instance=odbInstance, 
                   labels=eleIDs, 
                   data=eleStress[ii].tolist())
    # save every 1000th time step
    if (ii % 1000 == 0):
      odb.save()
      # gc.collect(generation=2) doesn't really help
  

odb.save()
odb.close()
