#include<stdlib.h>
#include<stdio.h>
#include<random>
#include<math.h>
#include "particle.hpp"
#include "constants.hpp"

/* memory operation */
// free memory
void freeParticles(particles_t * particles)
{
  delete[] particles->xyzs;
  particles->xyzs = NULL;
  delete[] particles->rs;
  particles->rs = NULL;
  delete[] particles->ms;
  particles->ms = NULL;
  delete[] particles->fs;
  particles->fs = NULL;
  delete[] particles->vs;
  particles->vs = NULL;
  delete[] particles->as;
  particles->as = NULL;
  delete[] particles->areas;
  particles->areas = NULL;
  delete[] particles->Ts;
  particles->Ts = NULL;
  /*
  free(particles->xyzs);
  particles->xyzs = NULL;
  free(particles->rs);
  particles->rs = NULL;
  free(particles->ms);
  particles->ms = NULL;
  free(particles->fs);
  particles->fs = NULL;
  free(particles->vs);
  particles->vs = NULL;
  free(particles->as);
  particles->as = NULL;
  free(particles->areas);
  particles->areas = NULL;
  free(particles->Ts);
  particles->Ts = NULL;
  */
}

void memgrowParticles(particles_t * particles)
{
  particles->n_particles = N_PARTICLES;
  particles->xyzs = allocFloatArray(N_PARTICLES*3);
  particles->rs = allocFloatArray(N_PARTICLES);
  particles->ms = allocFloatArray(N_PARTICLES);
  particles->fs = allocFloatArray(N_PARTICLES*3);
  particles->vs = allocFloatArray(N_PARTICLES*3);
  particles->as = allocFloatArray(N_PARTICLES*3);
  particles->areas = allocFloatArray(N_PARTICLES);
  particles->Ts = allocFloatArray(N_PARTICLES);
  /*
  particles->xyzs = initFloatArray(particles->n_particles*3, false);
  particles->rs = initFloatArray(particles->n_particles, false);
  particles->ms = initFloatArray(particles->n_particles, false);
  particles->fs = initFloatArray(particles->n_particles*3, false);
  particles->vs = initFloatArray(particles->n_particles*3, true); // zero init
  particles->as = initFloatArray(particles->n_particles*3, false);
  particles->areas = initFloatArray(particles->n_particles, true);
  particles->Ts = initFloatArray(particles->n_particles, true);
  */
}


/* allocate and init float array on the heap with new */
float_value_t * allocFloatArray(int num_elements)
{
  float_value_t * array = new float_value_t[num_elements]();
  return array;
}

/* allocate and init float array on the heap with malloc/calloc */
float_value_t* initFloatArray(int num_elements, bool zero)
{
  float_value_t * array = NULL;
  if (zero)
  {
    array = (float_value_t *) calloc(num_elements, sizeof(float_value_t));
  } else
  {
    array = (float_value_t *) malloc(sizeof(float_value_t) * num_elements);
  }
  if(array == NULL)
  {
    fprintf(stderr, "malloc fails\n");
    exit(1);
  } else
  {
    return array;
  }
}

/*
// random particle diameter
float_value_t randomParticleDiameter(std::default_random_engine generator, std::normal_distribution<float_value_t> diameter_distribution)
{
  float_value_t diameter;
  do
  {
    diameter = diameter_distribution(generator);
  } while (diameter < ELEMENT_DIAMETER_MIN || diameter > ELEMENT_DIAMETER_MAX)
  return diameter;
}
*/

/* initialize particles */
void initParticles(float_value_t length, float_value_t width, particles_t * particles)
{
  memgrowParticles(particles);
  float_value_t mu = ELEMENT_DIAMETER_MEAN;
  float_value_t sigma = ELEMENT_DIAMETER_STD;
  float_value_t min_height = MIN_HEIGHT;
  float_value_t max_height = MAX_HEIGHT;
  bool find_particle = false;
  float_value_t x, y, z, diameter, radius;
  srand(time(NULL));
  std::default_random_engine generator(rand());
  std::normal_distribution<float_value_t> diameter_distribution(mu, sigma);
  std::uniform_real_distribution<float_value_t> x_distribution(0.0, length);
  std::uniform_real_distribution<float_value_t> y_distribution(0.0, width);
  std::uniform_real_distribution<float_value_t> z_distribution(min_height, max_height);
  for (int i = 0; i < particles->n_particles; i++)
  {
    while(!find_particle)
    {
      while(true)
      {
        x = x_distribution(generator);
        y = y_distribution(generator);
        z = z_distribution(generator);
        do
        {
          diameter = diameter_distribution(generator);
        } while (diameter < ELEMENT_DIAMETER_MIN || diameter > ELEMENT_DIAMETER_MAX);
        radius = 0.5*diameter;
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
        if(dist < radius + particles->rs[j])
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

