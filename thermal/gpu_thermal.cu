#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>

#include "particle.hpp"
#include "constants.hpp"
#include "dynamics.hpp"
#include "utility.hpp"

#define NUM_THREADS 1024

#include "cuda_thermal.cuh"
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
  float_value_t r = sqrt((px - lx)*(px - lx) + (py - ly)*(py - ly));
  float_value_t z = fabs(pz - lz); // distance of laser and particle center
  float_value_t I = computeIntensity(laser_ptr, r, z, 2*particles_ptr->rs[i]);
  return ABSORPTIVITY*I*particles_ptr->areas[i];
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
  float_value_t T = particles_ptr->Ts[index];
  return EPSILON_SB * (T*T*T*T - QUAD_ENVIORNMENT_T) * particles_ptr->areas[index];
}

  __device__
float_value_t computeConvection(particles_t * particles_ptr, int index)
{
  return HEAT_TRANSFER_COEFFICIENT * (particles_ptr->Ts[index] - ENVIORNMENT_T) * particles_ptr->areas[index];
}

__device__
float_value_t computeConductionToFloor(particles_t * particles_ptr, int i)
{
  float_value_t area = FLOOR_CONTACT_RATIO * particles_ptr->areas[i];
  float_value_t T = particles_ptr->Ts[i];
  return CONDUCTIVITY * (T - TABLE_T) / (particles_ptr->rs[i]) * area;
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
        do
        {
          diameter = diameter_distribution(generator);
        } while (diameter < ELEMENT_DIAMETER_MIN || diameter > ELEMENT_DIAMETER_MAX);
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

/*
typedef struct
{
  float_value_t intensity;
  float_value_t cx; // center x
  float_value_t cy; // center y
  float_value_t cz; // center z
} laser_t;

// initialize laser
extern __host__ __device__
laser_t* initLaser(float_value_t new_x, float_value_t new_y, float_value_t new_z);

// update laser position
extern __host__ __device__
void updateLaser(laser_t * laser_ptr, float_value_t new_x, float_value_t new_y, float_value_t new_z);

// print laser position
extern __host__ __device__
void printLaser(laser_t * laser_ptr);

extern __device__
float_value_t computeH(laser_t * laser_ptr, particles_t * particles_ptr, int i);

extern __device__
float_value_t computeIntensity(laser_t * laser_ptr, float_value_t r, float_value_t z, float_value_t D);

extern __device__
float_value_t computeBeta(float_value_t D);

extern __device__
float_value_t computeDistance(particles_t * particles_ptr, int i, int j);

extern __device__
float_value_t computeCapacity(float_value_t T);

extern __device__
float_value_t computeContactArea(float_value_t d, float_value_t Ri, float_value_t Rj);

extern __device__
float_value_t computeRadiation(particles_t * particles_ptr, int index);

extern __device__
float_value_t computeConvection(particles_t * particles_ptr, int index);

extern __device__
float_value_t computeConduction(particles_t * particles_ptr, int i, int j);

extern __host__
void initParticlesAfterDynamics(float_value_t length, float_value_t width, particles_t * particles);
*/

__global__
void thermalkernel(float_value_t x, float_value_t y, float_value_t z, particles_t* d_particles_ptr, laser_t* d_laser_ptr)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  // use the first thread to update laser
  if(index == 0)
  {
    updateLaser(d_laser_ptr, x, y, z);
    // printf("x = %0.12f\n", x);
  }
  __syncthreads();
  if (index < N_PARTICLES)
  {
    float_value_t Ti = d_particles_ptr->Ts[index];
    float_value_t mi = d_particles_ptr->ms[index];
    float_value_t ci = computeCapacity(Ti);
    float_value_t Hi = computeH(d_laser_ptr, d_particles_ptr, index);
    float_value_t Qi_conv = computeConvection(d_particles_ptr, index);
    float_value_t Qi_rad = computeRadiation(d_particles_ptr, index);

    float_value_t Qi_floor = computeConductionToFloor(d_particles_ptr, index); // add floor
    Qi_floor = 0.0;
    float_value_t Qi_cond = 0.0f;
    for (int j = 0; j < N_PARTICLES; j++)
    {
      if (j == index) continue;
      Qi_cond += computeConduction(d_particles_ptr, index, j);
    }
    // float_value_t deltTi = THERMAL_TIME_STEP*(Hi - 4*Qi_conv - 4*Qi_rad - Qi_cond)/(mi*ci);
    float_value_t deltTi = THERMAL_TIME_STEP*(Hi - Qi_conv - Qi_rad - Qi_cond - Qi_floor)/(mi*ci);
    __syncthreads();
    d_particles_ptr->Ts[index] += deltTi;
    /*
    if(index == 7)
    {
      printf("Hi = %.18f, Qi_conv = %.18f, Qi_rad = %.18f, Qi_cond = %.18f\n", Hi, Qi_conv, Qi_rad, Qi_cond);
      printf("T = %0.6f\n", d_particles_ptr->Ts[index]);
    }
    */
  } // end if
}

int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    fprintf(stderr, "Usage: >> ./gpu_thermal test_type\n");
    fprintf(stderr, "test type = 1: load xyzs.txt and radius.txt from infiles/ \n");
    fprintf(stderr, "test type = 2: initialize particle position and radius randomly \n");
    exit(1);
  }

  float_value_t length = TABLE_LENGTH;
  float_value_t width = TABLE_WIDTH;

  cudaEvent_t start_event, end_event;
  cudaEventCreate(&start_event);
  cudaEventCreate(&end_event);

  particles_t* h_particles_ptr = (particles_t *) malloc(sizeof(particles_t));

  initParticlesAfterDynamics(length, width, h_particles_ptr);

  int test_type = atoi(argv[1]);


  if (test_type == 2) // initialize particle randomly, write positions and radius to files
  {
    std::string coord_file = "outfiles/xyzs_after.txt";
    std::string radius_file = "outfiles/radius_after.txt";
    writeFile(h_particles_ptr->xyzs, h_particles_ptr->n_particles, 3, coord_file);
    writeFile(h_particles_ptr->rs, h_particles_ptr->n_particles, 1, radius_file);
  } else if (test_type == 1) // load xyzs.txt and radius.txt
  {
    std::string coord_file = "infiles/xyzs.txt";
    std::string radius_file = "infiles/radius.txt";
    loadParticlesAfterDynamics(coord_file, radius_file, h_particles_ptr);
  } else
  {
    fprintf(stderr, "Unknown test type\n");
    exit(1);
  }

  cudaEventRecord(start_event, 0);
  // float_value_t z_zero = LAYER_HEIGHT;
  float_value_t z_zero;
  if (test_type == 2)
  {
    z_zero = 0.000030;
  } else
  {
    z_zero = 0.000120;
  }

  float_value_t x_zero = 0.0;
  float_value_t y_zero = 0.5*TABLE_WIDTH;

  laser_t* h_laser_ptr = initLaser(x_zero, y_zero, z_zero);
  printf("%s\n", "GPU simulation:");
  printf("number of particles = %d\n", h_particles_ptr->n_particles);
  printf("Laser intensity = %f\n", h_laser_ptr->intensity);
  printf("start simulation\n");

  laser_t* d_laser_ptr;
  cudaMalloc((laser_t **)& d_laser_ptr, sizeof(laser_t));
  cudaMemcpy(d_laser_ptr, h_laser_ptr, sizeof(laser_t), cudaMemcpyHostToDevice);

  float_value_t * d_ms;
  float_value_t * d_rs;
  float_value_t * d_xyzs;
  float_value_t * d_areas;
  float_value_t * d_Ts;

  const unsigned int bytes = N_PARTICLES*sizeof(float_value_t);
  const unsigned int bytes3 = 3*bytes;

  cudaMalloc((float_value_t **)& d_ms, bytes);
  cudaMalloc((float_value_t **)& d_rs, bytes);
  cudaMalloc((float_value_t **)& d_xyzs, bytes3);
  cudaMalloc((float_value_t **)& d_areas, bytes);

  cudaMalloc((float_value_t **)& d_Ts, bytes);

  cudaMemcpy(d_ms, h_particles_ptr->ms, bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_rs, h_particles_ptr->rs, bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_xyzs, h_particles_ptr->xyzs, bytes3, cudaMemcpyHostToDevice);
  cudaMemcpy(d_areas, h_particles_ptr->areas, bytes, cudaMemcpyHostToDevice);

  cudaMemcpy(d_Ts, h_particles_ptr->Ts, bytes, cudaMemcpyHostToDevice);

  particles_t * d_particles_ptr;
  cudaMalloc((particles_t**)&d_particles_ptr, sizeof(particles_t));
  cudaMemcpy(&(d_particles_ptr->n_particles), &(h_particles_ptr->n_particles), sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(&(d_particles_ptr->ms), &d_ms, sizeof(d_particles_ptr->ms), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_particles_ptr->rs), &d_rs, sizeof(d_particles_ptr->rs), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_particles_ptr->xyzs), &d_xyzs, sizeof(d_particles_ptr->xyzs), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_particles_ptr->areas), &d_areas, sizeof(d_particles_ptr->areas), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_particles_ptr->Ts), &d_Ts, sizeof(d_particles_ptr->Ts), cudaMemcpyHostToDevice);

  // thermal analysis
  float_value_t current_time = 0.0;
  float_value_t SCAN_LENGTH = TABLE_LENGTH * SCAN_LENGTH_RATIO;
  float_value_t end_time = SCAN_LENGTH/LASER_SCAN_SPEED;

  int count = 0;

  while(current_time < end_time)
  {
    current_time += THERMAL_TIME_STEP;
    float_value_t current_x = LASER_SCAN_SPEED*(current_time - 0.5*THERMAL_TIME_STEP);
    if (current_x > SCAN_LENGTH) break;
    int num_blocks = (N_PARTICLES + NUM_THREADS - 1)/NUM_THREADS;
    thermalkernel<<<num_blocks, NUM_THREADS>>>(current_x, y_zero, z_zero, d_particles_ptr, d_laser_ptr);

    cudaDeviceSynchronize();

    count += 1;
    // if (count == 2) break;
  }
  printf("number of steps = %d\n", count);

  cudaMemcpy(h_particles_ptr->Ts, d_Ts, bytes, cudaMemcpyDeviceToHost);

  cudaFree(d_ms);
  cudaFree(d_rs);
  cudaFree(d_xyzs);
  cudaFree(d_areas);
  cudaFree(d_Ts);
  cudaFree(d_particles_ptr);
  cudaFree(d_laser_ptr);

  printf("end simulation\n");

  std::string temp_file = "outfiles/temperatures_cuda.txt";
  writeFile(h_particles_ptr->Ts, h_particles_ptr->n_particles, 1, temp_file);

  // free memory
  freeParticles(h_particles_ptr);
  free(h_particles_ptr);
  h_particles_ptr = NULL;
  cudaEventRecord(end_event, 0);
  cudaEventSynchronize(end_event);
  float elapsedTime;
  cudaEventElapsedTime(&elapsedTime, start_event, end_event);
  printf("%0.6f\n", elapsedTime/1000.0);
  return 0;
}
