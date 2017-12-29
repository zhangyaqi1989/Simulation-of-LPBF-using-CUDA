#ifndef THERMAL_HPP
#define THERMAL_HPP

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
laser_t* initLaser(float_value_t new_x, float_value_t new_y, float_value_t new_z);

// update laser position
void updateLaser(laser_t * laser_ptr, float_value_t new_x, float_value_t new_y, float_value_t new_z);

// print laser position
void printLaser(laser_t * laser_ptr);

// gaussian laser
// H_{i} = alpha * I(r, z) * A_{i}
float_value_t computeH(laser_t * laser_ptr, particles_t * particles_ptr, int i);
// I(r, z) is the laser intensive
// r and z : radial distance and depth from the beam center
float_value_t computeIntensity(laser_t * laser_ptr, float_value_t r, float_value_t z, float_value_t D);
float_value_t computeBeta(float_value_t D);

float_value_t computeDistance(particles_t * particles_ptr, int i, int j);
float_value_t computeCapacity(float_value_t T);
float_value_t computeContactArea(float_value_t d, float_value_t Ri, float_value_t Rj);

float_value_t computeRadiation(particles_t * particles_ptr, int index);
float_value_t computeConvection(particles_t * particles_ptr, int index);
float_value_t computeConduction(particles_t * particles_ptr, int i, int j);
float_value_t computeConductionToFloor(particles_t * particles_ptr, int i);

/* init particles after dynamic analysis */
void initParticlesAfterDynamics(float_value_t length, float_value_t width, particles_t * particles);



#endif
