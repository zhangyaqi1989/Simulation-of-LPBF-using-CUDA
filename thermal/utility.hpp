/* written by Yaqi Zhang
 * Dec 2017
 */
#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <iostream>
#include"particle.hpp"

/* utility functions */
void writeFile(float_value_t* array, int row, int col, std::string filename);


/* read file from an file */
void readFile(std::string filename, float_value_t * array, int nrows, int ncols);


/* load particles from the output coord_file and radius_file of dynamic analysis */
void loadParticlesAfterDynamics(std::string coord_file, std::string radius_file, particles_t* particles);

#endif

