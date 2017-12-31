#!/bin/bash
# this script is used for testing
# it runs both cpu version and gpu version of thermal simulation 
# and compare the results
if [ $1 == 1 ]
then
  echo load xyzs.txt and radius.txt from infiles/
else
  echo initialize particle position and radius randomly
fi
./cpu_thermal $1
./gpu_thermal $1
diff outfiles/temperatures.txt outfiles/temperatures_cuda.txt
