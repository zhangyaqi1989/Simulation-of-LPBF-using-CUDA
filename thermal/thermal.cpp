/* written by Yaqi Zhang
 * Dec 2017
 */
#include "thermal.hpp"
#include "constants.hpp"
#include "utility.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <math.h>

// define some global variables
const float_value_t C_SL = (C_SOLID + C_LIQUID)/2.0 + MELTING_LATENT_HEAT/DELTA_T;
const float_value_t C_LG = (C_LIQUID + C_GAS)/2.0 + VAPORIZATION_LATENT_HEAT/DELTA_T;
const float_value_t QUAD_ENVIORNMENT_T = ENVIORNMENT_T * ENVIORNMENT_T * ENVIORNMENT_T * ENVIORNMENT_T;
const float_value_t EPSILON_SB = EMMISIVITY*STEFAN_BOLTZMANN;


// initialize laser
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
void updateLaser(laser_t * laser_ptr, float_value_t new_x, float_value_t new_y, float_value_t new_z)
{
  laser_ptr->cx = new_x;
  laser_ptr->cy = new_y;
  laser_ptr->cz = new_z;
}


// print laser position
void printLaser(laser_t * laser_ptr)
{
  printf("@(%.3f, %.3f, %.3f)\n", laser_ptr->cx, laser_ptr->cy, laser_ptr->cz);
}


/* gaussian laser
 * H_{i} = alpha * I(r, z) * A_{i}
 */
float_value_t computeH(laser_t * laser_ptr, particles_t * particles_ptr, int i)
{
  float_value_t px = particles_ptr->xyzs[i*3];
  float_value_t py = particles_ptr->xyzs[i*3 + 1];
  float_value_t pz = particles_ptr->xyzs[i*3 + 2];
  float_value_t lx = laser_ptr->cx;
  float_value_t ly = laser_ptr->cy;
  float_value_t lz = laser_ptr->cz;
  // float_value_t r = sqrt((px - lx)*(px - lx) + (py - ly)*(py - ly) + (pz - lz)*(pz - lz)); // radial distance does not include z direction
  float_value_t r = sqrt((px - lx)*(px - lx) + (py - ly)*(py - ly));
  // float_value_t z = fabs(pz - lz) + particles_ptr->rs[i];
  float_value_t z = fabs(pz - lz); // distance of laser and particle center 
  float_value_t I = computeIntensity(laser_ptr, r, z, 2*particles_ptr->rs[i]);
  // return POROSITY*I*particles_ptr->areas[i]; // this is a bug
  return ABSORPTIVITY*I*particles_ptr->areas[i];
}


/* compute I(r, z)
 * I(r, z) is the laser intensive
 * r and z : radial distance and depth from the beam center
 * I(r, z) = I_{0} * e^{-beta*z} * e^{-2r^{2}/omega^{2}}
 */
float_value_t computeIntensity(laser_t * laser_ptr, float_value_t r, float_value_t z, float_value_t D)
{
  float_value_t beta = computeBeta(D);
  return laser_ptr->intensity * pow(M_E, -beta*z) * pow(M_E, -2*r*r/(LASER_SPOT_SIZE*LASER_SPOT_SIZE));
}


// compute beta based on particle diameter
float_value_t computeBeta(float_value_t D)
{
  return 3*(1 - POROSITY)/(2 * POROSITY * D);
}

/* compute capacity based on temperature */
float_value_t computeCapacity(float_value_t T)
{
  if (T < MELTING_T - DELTA_T/2.0)
    return C_SOLID;
  else if (T <= MELTING_T + DELTA_T/2.0)
    return C_SL;
  else if (T < BOILING_T - DELTA_T/2.0)
    return C_LIQUID;
  else if (T <= BOILING_T + DELTA_T/2.0)
    return C_LG;
  else
    return C_LG; //C_GAS;
}


/* compute distance of particle i and j */
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


/* compute contact area of particle i and j */
float_value_t computeContactArea(float_value_t d, float_value_t Ri, float_value_t Rj)
{
  float_value_t cos_value = (Rj*Rj - Ri*Ri - d*d)/(-2*d*Ri);
  float_value_t theta = acos(cos_value);
  float_value_t h = Ri * sin(theta);
  return M_PI*h*h;
}


/* compute radiation */
float_value_t computeRadiation(particles_t * particles_ptr, int index)
{
  float_value_t T = particles_ptr->Ts[index];
  return EPSILON_SB * (T*T*T*T - QUAD_ENVIORNMENT_T) * particles_ptr->areas[index];
}


/* compute convection */
float_value_t computeConvection(particles_t * particles_ptr, int index)
{
  return HEAT_TRANSFER_COEFFICIENT * (particles_ptr->Ts[index] - ENVIORNMENT_T) * particles_ptr->areas[index];
}


/* compute conduction to floor (not used) */
float_value_t computeConductionToFloor(particles_t * particles_ptr, int i)
{
  float_value_t z = particles_ptr->xyzs[i*3 + 2];
  float_value_t r = particles_ptr->rs[i];
  if (fabs(z - r) > 1e-12)  return 0.0;
  float_value_t area = FLOOR_CONTACT_RATIO * particles_ptr->areas[i];
  float_value_t T = particles_ptr->Ts[i];
  return CONDUCTIVITY * (T - TABLE_T) / (particles_ptr->rs[i]) * area;
}


/* compute heat conduction between particle i and j */
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
void initParticlesAfterDynamics(float_value_t length, float_value_t width, particles_t * particles)
{
  memgrowParticles(particles);
  float_value_t mu = ELEMENT_DIAMETER_MEAN;
  float_value_t sigma = ELEMENT_DIAMETER_STD;
  // float_value_t min_height = MIN_HEIGHT;
  // float_value_t max_height = MAX_HEIGHT;
  bool find_particle = false;
  float_value_t x, y, z, diameter, radius;
  // srand(time(NULL));
  srand(1);
  std::default_random_engine generator(rand());
  std::normal_distribution<float_value_t> diameter_distribution(mu, sigma);
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
        diameter = diameter_distribution(generator);
        radius = 0.5*diameter;
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
