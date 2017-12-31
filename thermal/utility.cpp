/* written by Yaqi Zhang
 * Dec 2017
 */
#include"utility.hpp"
#include "constants.hpp"
#include<fstream>
#include<math.h>

/* write an array to a file  according to the row and col */
void writeFile(float_value_t* array, int row, int col, std::string filename)
{
  std::ofstream myfile;
  myfile.open(filename);
  for(int i = 0; i < row; i++)
  {
    for(int j = 0; j < col; j++)
    {
      myfile << array[i * col + j] << " ";
    }
    myfile << "\n";
  }
  std::cout << "Writing data to " << filename << ".\n";
  myfile.close();
}


/* read file from an file */
void readFile(std::string filename, float_value_t * array, int nrows, int ncols)
{
  std::ifstream infile;
  infile.open(filename);
  int line = 0;
  while(!infile.eof())
  {
    for(int i = 0; i < ncols; i++)
    {
      infile >> array[line*ncols + i];
    }
    line++;
  }
  infile.close();
}


/* load particles from the output coord_file and radius_file of dynamic analysis */
void loadParticlesAfterDynamics(std::string coord_file, std::string radius_file, particles_t* particles)
{
  memgrowParticles(particles);
  int n_particles = particles->n_particles;
  readFile(coord_file, particles->xyzs, n_particles, 3);
  readFile(radius_file, particles->rs, n_particles, 1);
  float_value_t radius;
  for (int i = 0; i < n_particles; i++)
  {
    radius = particles->rs[i];
    particles->ms[i] = 4.0/3.0*M_PI*radius*radius*radius;
    particles->areas[i] = M_PI*radius*radius;
    particles->Ts[i] = PREHEAT_T;
  }
}
