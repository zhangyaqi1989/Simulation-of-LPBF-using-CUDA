#include<iostream>
#include<stdio.h>
#include<stdlib.h>

#include "particle.hpp"
#include "constants.hpp"
#include "dynamics.hpp"
#include "thermal.hpp"
#include "utility.hpp"

int main(void)
{
  float_value_t length = TABLE_LENGTH;
  float_value_t width = TABLE_WIDTH;
  particles_t* particles = (particles_t *) malloc(sizeof(particles_t));

  initParticles(length, width, particles);
  // write to file
  std::string coord_file = "xyzs.txt";
  std::string radius_file = "radius.txt";
  writeFile(particles->xyzs, particles->n_particles, 3, coord_file);
  writeFile(particles->rs, particles->n_particles, 1, radius_file);

  // simulation of the deposition of powder particles
  double current_time = 0.0f;
  while (current_time < 0.1)
  {
    updatePositions(particles);
    current_time += DYNAMIC_TIME_STEP;
  }
  // simulation of temperature evolution
  //
  // free memory
  freeParticles(particles);
  free(particles);
  particles = NULL;
  return 0;
}



