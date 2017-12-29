#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "particle.hpp"
#include "constants.hpp"
#include "dynamics.hpp"
#include "thermal.hpp"
#include "utility.hpp"

#define SCAN_LENGTH_RATIO (0.2)

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    fprintf(stderr, "Usage: >> ./cpu_thermal test_type\n");
    fprintf(stderr, "test type = 1: load xyzs.txt and radius.txt from infiles/ \n");
    fprintf(stderr, "test type = 2: initialize particle position and radius randomly \n");
    exit(1);
  }
  clock_t start_t, end_t;
  start_t = clock();

  float_value_t length = TABLE_LENGTH;
  float_value_t width = TABLE_WIDTH;
  particles_t* particles = (particles_t *) malloc(sizeof(particles_t));

  initParticlesAfterDynamics(length, width, particles); // initialize particles randomly

  int test_type = atoi(argv[1]);


  if (test_type == 2) // initialize particle randomly, write positions and radius to files
  {
    std::string coord_file = "outfiles/xyzs_after.txt";
    std::string radius_file = "outfiles/radius_after.txt";
    writeFile(particles->xyzs, particles->n_particles, 3, coord_file);
    writeFile(particles->rs, particles->n_particles, 1, radius_file);
  } else if (test_type == 1) // load xyzs.txt and radius.txt
  {
    std::string coord_file = "infiles/xyzs.txt";
    std::string radius_file = "infiles/radius.txt";
    loadParticlesAfterDynamics(coord_file, radius_file, particles);
  } else
  {
    fprintf(stderr, "Unknown test type\n");
    exit(1);
  }

  /*
  float_value_t max_z_temp = -100.0;
  for(int i = 0; i < particles->n_particles; i++)
  {
    if (particles->xyzs[i*3 + 2] > max_z_temp)
      max_z_temp = particles->xyzs[i*3 + 2];
  }
  */

  // float_value_t z_zero = LAYER_HEIGHT;
  // float_value_t z_zero = LAYER_HEIGHT;
  float_value_t z_zero;
  if (test_type == 2)
  {
    z_zero = 0.000030;
  } else
  {
    z_zero = 0.000120;
  }
  // printf("max_z = %f\n", max_z_temp);
  float_value_t x_zero = 0.0;
  float_value_t y_zero = 0.5*TABLE_WIDTH;

  laser_t* laser_ptr = initLaser(x_zero, y_zero, z_zero);
  printf("number of particles = %d\n", particles->n_particles);
  printf("intensity = %f\n", laser_ptr->intensity);
  printf("start simulation\n");

  // thermal analysis
  float_value_t current_time = 0.0;
  float_value_t SCAN_LENGTH = TABLE_LENGTH * SCAN_LENGTH_RATIO;
  float_value_t end_time = SCAN_LENGTH/LASER_SCAN_SPEED;
  int count = 0;
  float_value_t conducts[N_PARTICLES] = {0};
  float_value_t temp_Ts[N_PARTICLES] = {0};
  while (current_time < end_time)
  {
    current_time += THERMAL_TIME_STEP;
    float_value_t current_x = LASER_SCAN_SPEED*(current_time - 0.5*THERMAL_TIME_STEP);
    if (current_x > SCAN_LENGTH) break;
    updateLaser(laser_ptr, current_x, y_zero, z_zero);

    // loop over all pair of i and j
    for (int i = 0; i < particles->n_particles; i++)
    {
      for (int j = i + 1; j < particles->n_particles; j++)
      {
        float_value_t condij = computeConduction(particles, i, j);
        conducts[i] += condij;
        conducts[j] -= condij;
      }
    }
    // loop over all particles
    for (int i = 0; i < particles->n_particles; i++)
    {
      float_value_t Ti = particles->Ts[i];
      float_value_t mi = particles->ms[i];
      float_value_t ci = computeCapacity(Ti);
      float_value_t Hi = computeH(laser_ptr, particles, i);
      float_value_t Qi_conv = computeConvection(particles, i);
      float_value_t Qi_rad = computeRadiation(particles, i);
      float_value_t Qi_cond = conducts[i];
      float_value_t Qi_floor = computeConductionToFloor(particles, i); // add floor
      Qi_floor = 0.0; // not add floor;

      // float_value_t deltTi = THERMAL_TIME_STEP*(Hi - 4*Qi_conv - 4*Qi_rad - Qi_cond)/(mi*ci);
      float_value_t deltTi = THERMAL_TIME_STEP*(Hi - Qi_conv - Qi_rad - Qi_cond - Qi_floor)/(mi*ci);
      temp_Ts[i] = Ti + deltTi;
    }

    // zero conducts and update temperatures
    for(int i = 0; i < particles->n_particles; i++)
    {
      conducts[i] = 0;
      particles->Ts[i] = temp_Ts[i];
    }

    count += 1;
    // printf("count = %d\n", count);
    // if (count == 2) break;
  }
  printf("number of steps = %d\n", count);
  printf("end simulation\n");
  std::string temp_file = "outfiles/temperatures.txt";
  writeFile(particles->Ts, particles->n_particles, 1, temp_file);

  // free memory
  freeParticles(particles);
  free(particles);
  particles = NULL;
  end_t = clock();
  double total_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("%0.6f\n", total_time);
  return 0;
}
