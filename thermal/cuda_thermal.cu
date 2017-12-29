#include "cuda_thermal.cuh"
#include "constants.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <math.h>

// define some global variables
__const__ float_value_t C_SL = (C_SOLID + C_LIQUID)/2.0 + MELTING_LATENT_HEAT/DELTA_T;
__const__ float_value_t C_LG = (C_LIQUID + C_GAS)/2.0 + VAPORIZATION_LATENT_HEAT/DELTA_T;
__const__ float_value_t QUAD_ENVIORNMENT_T = ENVIORNMENT_T * ENVIORNMENT_T * ENVIORNMENT_T * ENVIORNMENT_T;
__const__ float_value_t EPSILON_SB = EMMISIVITY*STEFAN_BOLTZMANN;

// initialize laser
__host__ __device__
laser_t* initLaser(float_value_t x_zero, float_value_t y_zero, float_value_t z_zero)
{
  laser_t * laser_ptr = (laser_t *) malloc(sizeof(laser_t));
  laser_ptr->intensity = (2*LASER_POWER)/(M_PI*LASER_SPOT_SIZE*LASER_SPOT_SIZE);
  laser_ptr->cx = x_zero;
  laser_ptr->cy = y_zero;
  laser_ptr->cz = z_zero;
  return laser_ptr;
}

// update laser position
__host__ __device__
void updateLaser(laser_t * laser_ptr, float_value_t new_x, float_value_t new_y, float_value_t new_z)
{
  laser_ptr->cx = new_x;
  laser_ptr->cy = new_y;
  laser_ptr->cz = new_z;
}

// print laser position
__host__ __device__
void printLaser(laser_t * laser_ptr)
{
  printf("@(%.3f, %.3f, %.3f)\n", laser_ptr->cx, laser_ptr->cy, laser_ptr->cz);
}

// gaussian laser
// H_{i} = alpha * I(r, z) * A_{i}
__device__
float_value_t computeH(laser_t * laser_ptr, particles_t * particles_ptr, int i)
{
  float_value_t px = particles_ptr->xyzs[i*3];
  float_value_t py = particles_ptr->xyzs[i*3 + 1];
  float_value_t pz = particles_ptr->xyzs[i*3 + 2];
  float_value_t lx = laser_ptr->cx;
  float_value_t ly = laser_ptr->cy;
  float_value_t lz = laser_ptr->cz;
  // float_value_t r = sqrt((px - lx)*(px - lx) + (py - ly)*(py - ly) + (pz - lz)*(pz - lz));
  float_value_t r = sqrt((px - lx)*(px - lx) + (py - ly)*(py - ly));// + (fabs(pz - lz) + particles_ptr->rs[i])*(fabs(pz - lz) + particles_ptr->rs[i]));
  float_value_t z = fabs(pz - lz) + particles_ptr->rs[i];
  float_value_t I = computeIntensity(laser_ptr, r, z, 2*particles_ptr->rs[i]);
  return POROSITY*I*particles_ptr->areas[i];
}

// compute I(r, z)
// I(r, z) is the laser intensive
// r and z : radial distance and depth from the beam center
// I(r, z) = I_{0} * e^{-beta*z} * e^{-2r^{2}/omega^{2}}
__device__
float_value_t computeIntensity(laser_t * laser_ptr, float_value_t r, float_value_t z, float_value_t D)
{
  float_value_t beta = computeBeta(D);
  return laser_ptr->intensity * pow(M_E, -beta*z) * pow(M_E, -2*r*r/(LASER_SPOT_SIZE*LASER_SPOT_SIZE));
}

// compute beta based on particle diameter
__device__
float_value_t computeBeta(float_value_t D)
{
  return 3*(1 - POROSITY)/(2 * POROSITY * D);
}

__device__
float_value_t computeCapacity(float_value_t T)
{
  // printf("C_SL = %f\n", C_SL);
  // printf("C_LG = %f\n", C_LG);
  if (T < MELTING_T - DELTA_T/2.0)
    return C_SOLID;
  else if (T <= MELTING_T + DELTA_T/2.0)
    return C_SL;
  else if (T < BOILING_T - DELTA_T/2.0)
    return C_LIQUID;
  else if (T <= BOILING_T + DELTA_T/2.0)
    return C_LG;
  else
    return C_LG;//C_GAS;
}

__device__
float_value_t computeDistance(particles_t * particles_ptr, int i, int j)
{
  float_value_t xi = particles_ptr->xyzs[i*3];
  float_value_t yi = particles_ptr->xyzs[i*3 + 1];
  float_value_t zi = particles_ptr->xyzs[i*3 + 2];
  float_value_t xj = particles_ptr->xyzs[j*3];
  float_value_t yj = particles_ptr->xyzs[j*3 + 1];
  float_value_t zj = particles_ptr->xyzs[j*3 + 2];
  return sqrt((xi - xj)*(xi - xj) + (yi - yj)*(yi - yj) + (zi - zj)*(zi - zj));
}

__device__
float_value_t computeContactArea(float_value_t d, float_value_t Ri, float_value_t Rj)
{
  float_value_t cos_value = (Rj*Rj - Ri*Ri - d*d)/(-2*d*Ri);
  float_value_t theta = acos(cos_value);
  float_value_t h = Ri * sin(theta);
  return M_PI*h*h;
}

__device__
float_value_t computeRadiation(particles_t * particles_ptr, int index)
{
  // printf("QUAD_ENVIORNMENT_T = %.16f\n", QUAD_ENVIORNMENT_T);
  // printf("EPSILON_SB = %.16f\n", EPSILON_SB);
  float_value_t T = particles_ptr->Ts[index];
  return EPSILON_SB * (T*T*T*T - QUAD_ENVIORNMENT_T) * particles_ptr->areas[index];
}

__device__
float_value_t computeConvection(particles_t * particles_ptr, int index)
{
  return HEAT_TRANSFER_COEFFICIENT * (particles_ptr->Ts[index] - ENVIORNMENT_T) * particles_ptr->areas[index];
}

__device__
float_value_t computeConduction(particles_t * particles_ptr, int i, int j)
{
  float_value_t d = computeDistance(particles_ptr, i, j);
  float_value_t Ri = particles_ptr->rs[i];
  float_value_t Rj = particles_ptr->rs[j];
  if (d <= Ri + Rj)
  {
    if (d <= (Ri + Rj)*0.99) d = 0.99*(Ri + Rj);
    float_value_t Aij = computeContactArea(d, Ri, Rj);
    return CONDUCTIVITY*(particles_ptr->Ts[i] - particles_ptr->Ts[j])*Aij/d;
  }
  else
    return 0.0f;
}


// init particles after dynamic analysis
/* initialize particles */
__host__
void initParticlesAfterDynamics(float_value_t length, float_value_t width, particles_t * particles)
{
  memgrowParticles(particles);
  float_value_t mu = MU;
  float_value_t sigma = SIGMA;
  // float_value_t min_height = MIN_HEIGHT;
  // float_value_t max_height = MAX_HEIGHT;
  bool find_particle = false;
  float_value_t x, y, z, radius;
  srand(time(NULL));
  std::default_random_engine generator(rand());
  std::normal_distribution<float_value_t> radius_distribution(mu, sigma);
  std::uniform_real_distribution<float_value_t> x_distribution(0.0, length);
  std::uniform_real_distribution<float_value_t> y_distribution(0.0, width);
  // std::uniform_real_distribution<float_value_t> z_distribution(min_height, max_height);
  for (int i = 0; i < particles->n_particles; i++)
  {
    while(!find_particle)
    {
      // random a particle within table
      while(true)
      {
        x = x_distribution(generator);
        y = y_distribution(generator);
        // z = z_distribution(generator);
        radius = radius_distribution(generator);
        z = radius;
        if(x + radius <= length && x - radius >= 0 && y + radius <= width && y - radius >= 0) break;
      }

      find_particle = true;
      for(int j = 0; j < i; j++)
      {
        float_value_t dist = 0.0f;
        dist += (x - particles->xyzs[j*3])*(x - particles->xyzs[j*3]);
        dist += (y - particles->xyzs[j*3 + 1])*(y - particles->xyzs[j*3 + 1]);
        dist += (z - particles->xyzs[j*3 + 2])*(z - particles->xyzs[j*3 + 2]);
        dist = sqrt(dist);
        if(dist < (radius + particles->rs[j])*0.6)
        {
          find_particle = false;
          break;
        }
      }
    }

    find_particle = false;
    particles->xyzs[i*3 + 0] = x;
    particles->xyzs[i*3 + 1] = y;
    particles->xyzs[i*3 + 2] = z;
    particles->rs[i] = radius;
    particles->ms[i] = 4.0/3.0*M_PI*radius*radius*radius;
    particles->areas[i] = M_PI*radius*radius;
    particles->Ts[i] = PREHEAT_T;
  }

}

