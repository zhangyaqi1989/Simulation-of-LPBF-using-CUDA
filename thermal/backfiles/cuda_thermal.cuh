#ifndef CUDA_THERMAL_H
#define CUDA_THERMAL_H

#include "particle.hpp"

typedef struct
{
  float_value_t intensity;
  float_value_t cx; // center x
  float_value_t cy; // center y
  float_value_t cz; // center z
} laser_t;

/*
typedef struct
{
  float_value_t C_SL; // capacity between solid and liquid
  float_value_t C_LG; // capacity between liquid and gas
} capacity_genor_t;
*/

// initialize laser
__host__ __device__
laser_t* initLaser(float_value_t new_x, float_value_t new_y, float_value_t new_z);

// update laser position
__host__ __device__
void updateLaser(laser_t * laser_ptr, float_value_t new_x, float_value_t new_y, float_value_t new_z);

// print laser position
__host__ __device__
void printLaser(laser_t * laser_ptr);

// gaussian laser
// H_{i} = alpha * I(r, z) * A_{i}
__device__
float_value_t computeH(laser_t * laser_ptr, particles_t * particles_ptr, int i);
// I(r, z) is the laser intensive
// r and z : radial distance and depth from the beam center
__device__
float_value_t computeIntensity(laser_t * laser_ptr, float_value_t r, float_value_t z, float_value_t D);
__device__
float_value_t computeBeta(float_value_t D);

__device__
float_value_t computeDistance(particles_t * particles_ptr, int i, int j);
__device__
float_value_t computeCapacity(float_value_t T);
__device__
float_value_t computeContactArea(float_value_t d, float_value_t Ri, float_value_t Rj);

__device__
float_value_t computeRadiation(particles_t * particles_ptr, int index);
__device__
float_value_t computeConvection(particles_t * particles_ptr, int index);
__device__
float_value_t computeConduction(particles_t * particles_ptr, int i, int j);

__host__
void initParticlesAfterDynamics(float_value_t length, float_value_t width, particles_t * particles);
#endif
