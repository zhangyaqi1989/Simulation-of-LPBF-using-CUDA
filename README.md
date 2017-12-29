# Simulation of Laser Powder Based Fusion (L-PBF) using CUDA
In this project, both dynamic and thermal simulation of laser powder bed fusion are implemented using CUDA
## 1. Dynamic Simulation
### 1.1 Build Instructions
On Linux

```
>> cd dynamics/
>> make -f makefile-pc clean
>> make -f makefile-pc

```
On Euler (a multi-core supercomputer cluster)

```
>> cd dynamics/
>> make -f makefile-euler clean
>> make -f makefile-euler
```
### 1.2 Run Instructions
CPU version
```
>> ./cpu_dynamic
```
GPU version
```
>>> ./gpu_dynamic
```

Both CPU and GPU version will 

* output particle radius to radius.txt

* output particle initial velocity to vxyzs0.txt

* output particle initial position to xyzs0.txt

* output particle final velocity to vxyz.txt

* output particle final position to xyzs.txt

## 1. Thermal Simulation
### 1.1 Build Instructions
On Linux

```
>> cd thermal/
>> make -f makefile-pc clean
>> make -f makefile-pc

```
On Euler (a multi-core supercomputer cluster)

```
>> cd thermal/
>> make -f makefile-euler clean
>> make -f makefile-euler
```
### 1.2 Run Instructions
CPU version
```
>> ./cpu_thermal <test_type>
```
GPU version
```
>>> ./gpu_thermal <test_type>
```
test_type = 1: load xyzs.txt and radius.txt from infiles/
test_type = 2: initialize particle position and radius randomly

CPU version will output temperatures of final time step to outfiles/temperatures.txt

GPU version will output temperatures of final time step to outfiles/temperatures_cuda.txt

when test_type = 2, both CPU and GPU version will also output outfiles/xyzs_after.txt and outfiles/radius_after.txt which are the particle positions and radiuses initialized randomly. 