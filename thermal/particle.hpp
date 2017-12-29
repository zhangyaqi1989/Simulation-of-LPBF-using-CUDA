#ifndef PARTICLE_HPP
#define PARTICLE_HPP
// #include<iostream>
typedef double float_value_t;

typedef struct
{
  int n_particles;
  float_value_t* xyzs; // position
  float_value_t* rs;   // radius
  float_value_t* ms;   // mass
  float_value_t* fs;   // force
  float_value_t* vs;   // velocity
  float_value_t* as;   // acceleration
  float_value_t* areas;// area A = pi*radius^2
  float_value_t* Ts; // temperature
} particles_t;

/* memory operations */
void freeParticles(particles_t * particles);
void memgrowParticles(particles_t * particles);

/* initialize particles */
float_value_t* initFloatArray(int num_elements, bool zero);
float_value_t* allocFloatArray(int num_elements);
void initParticles(float_value_t length, float_value_t width, particles_t * particles);


#endif

