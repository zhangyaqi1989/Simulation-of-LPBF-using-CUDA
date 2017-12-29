// deposition_sim.cpp : Defines the entry point for the console application.
// Xin Liu


#include <string.h> 
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <random>
#include <math.h>
#include <time.h>
#include "stopwatch.hpp"


#define C_bar (0.000027)
#define boundary_left (0)
#define boundary_right (0.001)
#define boundary_front (0)
#define boundary_rear (0.0005)
#define boundary_up (0.003)
#define boundary_bottom (0.000)

#define N_PARTICLES (512)
#define MAX_HEIGHT (0.00020)
#define MIN_HEIGHT (0.00014)
#define MU (0.000027)
#define SIGMA (0.00001/2.355)

#define zeta (0.1) 
#define mu_d (0.1)
#define mu_f (2.1e-3)
#define g_grav (-9.8)
#define rho (7952)

#define E_modulus (193e9)
#define nu (0.26)

#define length boundary_right
#define width boundary_rear

#define MAX_R 10

#define C_factor 0.0055  //sqrt(1/(4/3*rho*pi))
#define m_bar (4.0/3*M_PI*C_bar*C_bar*C_bar)
#define omega_bar (sqrt(E_modulus/(2-2*nu*nu)*C_bar/m_bar))

#define v0 (-00)
#define TIME_STEP (5e-8*omega_bar) //maximal limit 0.0000001
#define TIME 0.2
typedef double float_value_t;

typedef struct
{
	int n_particles;
	float_value_t* xyzs;
	float_value_t* rs;
	float_value_t* ms;
	float_value_t* fs;
	float_value_t* vs;
	float_value_t* as;
	float_value_t* E;
	float_value_t* mu;
} particles_t;

void writeFile(float_value_t* array, int row, int col, std::string filename);
void initParticles(particles_t * particles);
float_value_t* initFloatArray(int num_elements);
void freeParticles(particles_t * particles);
void updatePositions(particles_t *particles, float_value_t t, float_value_t delta_t);
// void updateAcceleration(particles_t * particles);
float_value_t norm(float_value_t* r)
{
	return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}
float_value_t dctn(float_value_t* r,int i)
{
	return r[i]/sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}
float_value_t dot(float_value_t* r, float_value_t* r1)
{
	return (r[0] * r1[0] + r[1] * r1[1] + r[2] * r1[2]);
}
void f(float_value_t* res, particles_t *part, float_value_t* x_input, float_value_t delta_t)
{
	for (int i = 0; i < N_PARTICLES; i++)
	{
		float_value_t R_i = part->rs[i];
		float_value_t E_i = part->E[i];
		float_value_t nu_i = part->mu[i];
		float_value_t m_i = part->ms[i];

		//float_value_t* x_i = &(part->xyzs[3 * i]);
		//float_value_t* x_i = &(x_input[3 * i]);
		float_value_t x_i[3];
		x_i[0] = x_input[3 * i]*C_bar;
		x_i[1] = x_input[3 * i+1] *C_bar;
		x_i[2] = x_input[3 * i+2] *C_bar;
		//float_value_t* dx_i = &(part->vs[3 * i]);
		//float_value_t* dx_i = &(x_input[3 * i + 3 * N_PARTICLES]);
		float_value_t dx_i[3];
		dx_i[0] = x_input[3 * i + 3 * N_PARTICLES]*omega_bar*C_bar;
		dx_i[1] = x_input[3 * i + 1 + 3 * N_PARTICLES] * omega_bar*C_bar;
		dx_i[2] = x_input[3 * i + 2 + 3 * N_PARTICLES] * omega_bar*C_bar;

		float_value_t F_env[3] = { 0 };
		float_value_t F_con[3] = { 0 };
		float_value_t F_fric[3] = { 0 };

		float_value_t F_con_iw[3] = { 0 };
		float_value_t R_aster, E_aster, delta_iw, v_delta_iw, d;
		float_value_t F_fric_iw[3] = { 0 };


		F_env[0] = -6*M_PI*mu_f*R_i*dx_i[0];
		F_env[1] = -6*M_PI*mu_f*R_i*dx_i[1];
		F_env[2] = -6*M_PI*mu_f*R_i*dx_i[2];
		
		if (x_i[0] - R_i <= boundary_left)
		{
			R_aster = R_i;
			E_aster = E_i / (2 * (1 - nu_i*nu_i));
			delta_iw = abs(x_i[0] - boundary_left - (R_i));
			v_delta_iw = abs(dx_i[0]);
			d = 2 * zeta*sqrt(2 * E_aster*m_i*sqrt(R_aster))*pow(abs(delta_iw), (0.25));
			F_con_iw[0] = 4.0 / 3 * sqrt(R_aster)*E_aster*pow(delta_iw, 1.5) + d*v_delta_iw;
			F_con_iw[1] = 0;
			F_con_iw[2] = 0;
		}
		if (x_i[0] + R_i >= boundary_right)
		{
			R_aster = R_i;
			E_aster = E_i / (2 * (1 - nu_i*nu_i));
			delta_iw = abs(x_i[0] - boundary_right + (R_i));
			v_delta_iw = abs(dx_i[0]);
			d = 2 * zeta*sqrt(2 * E_aster* m_i*sqrt(R_aster))*pow(abs(delta_iw), (0.25));
			F_con_iw[0] = -4.0 / 3 * sqrt(R_aster)*E_aster*pow(delta_iw, 1.5) - d*v_delta_iw;
			F_con_iw[1] = 0;
			F_con_iw[2] = 0;
		}
		if (x_i[1] - R_i <= boundary_front)
		{
			R_aster = R_i;
			E_aster = E_i / (2 * (1 - nu_i*nu_i));
			delta_iw = abs(x_i[1] - boundary_front - (R_i));
			v_delta_iw = abs(dx_i[1]);
			d = 2 * zeta*sqrt(2 * E_aster* m_i*sqrt(R_aster))*pow(abs(delta_iw), (0.25));
			F_con_iw[0] = 0;
			F_con_iw[1] = 4.0 / 3 * sqrt(R_aster)*E_aster*pow(delta_iw, 1.5) + d*v_delta_iw;
			F_con_iw[2] = 0;
		}
		if (x_i[1] + R_i >= boundary_rear)
		{
			R_aster = R_i;
			E_aster = E_i / (2 * (1 - nu_i*nu_i));
			delta_iw = abs(x_i[1] - boundary_rear + (R_i));
			v_delta_iw = abs(dx_i[1]);
			d = 2 * zeta*sqrt(2 * E_aster*  m_i*sqrt(R_aster))*pow(abs(delta_iw), (0.25));
			F_con_iw[0] = 0;
			F_con_iw[1] = -4.0 / 3 * sqrt(R_aster)*E_aster*pow(delta_iw, 1.5) - d*v_delta_iw;
			F_con_iw[2] = 0;
		}
		if (x_i[2] - R_i <= boundary_bottom)
		{
			R_aster = R_i;
			E_aster = E_i / (2.0 * (1 - nu_i*nu_i));
			delta_iw = abs(x_i[2] - boundary_bottom - (R_i));
			v_delta_iw = abs(dx_i[2]);
			d = 2 * zeta*sqrt(2 * E_aster* m_i*sqrt(R_aster))*pow(abs(delta_iw), (0.25));
			F_con_iw[0] = 0;
			F_con_iw[1] = 0;
			F_con_iw[2] = 4.0 / 3 * sqrt(R_aster)*E_aster*pow(delta_iw, 1.5) + d*v_delta_iw;

			//float_value_t normdxi = sqrt(dx_i[0] * dx_i[0] + dx_i[1] * dx_i[1] + dx_i[2] * dx_i[2])
		}
		// assume only collision with one wall
		if (norm(dx_i) <= 1e-30)
		{
			F_fric_iw[0] = 0;
			F_fric_iw[1] = 0;
			F_fric_iw[2] = 0;
		}
		else
		{
			F_fric_iw[0] = mu_d*norm(F_con_iw)*(-dx_i[0]) / norm(dx_i);
			F_fric_iw[1] = mu_d*norm(F_con_iw)*(-dx_i[1]) / norm(dx_i);
			F_fric_iw[2] = mu_d*norm(F_con_iw)*(-dx_i[2]) / norm(dx_i);
		}
		F_con[0] = F_con[0] + F_con_iw[0];
		F_con[1] = F_con[1] + F_con_iw[1];
		F_con[2] = F_con[2] + F_con_iw[2];
		F_fric[0] = F_fric[0] + F_fric_iw[0];
		F_fric[1] = F_fric[1] + F_fric_iw[1];
		F_fric[2] = F_fric[2] + F_fric_iw[2];
		int k = 0;
		for (int j = 0; j < N_PARTICLES; j++)
		{
			float_value_t R_j = part->rs[j];
			float_value_t E_j = part->E[j];
			float_value_t nu_j = part->mu[j];
			float_value_t m_j = part->ms[j];


			float_value_t R_aster = R_i*R_j / (R_i + R_j);
			float_value_t E_aster = E_i*E_j / (E_i*(1 - nu_i  *nu_i) + E_j*(1 - nu_j *nu_j));


			float_value_t x_j[3];
			x_j[0] = x_input[3 * j] *C_bar;
			x_j[1] = x_input[3 * j + 1] * C_bar;
			x_j[2] = x_input[3 * j + 2] * C_bar;
			//float_value_t* dx_i = &(part->vs[3 * i]);
			//float_value_t* dx_i = &(x_input[3 * i + 3 * N_PARTICLES]);
			float_value_t dx_j[3];
			dx_j[0] = x_input[3 * j + 3 * N_PARTICLES] * omega_bar*C_bar;
			dx_j[1] = x_input[3 * j + 1 + 3 * N_PARTICLES] * omega_bar*C_bar;
			dx_j[2] = x_input[3 * j + 2 + 3 * N_PARTICLES] * omega_bar*C_bar;


			float_value_t m_aster;

			float_value_t x_ij[3];
			x_ij[0] = x_j[0] - x_i[0];
			x_ij[1] = x_j[1] - x_i[1];
			x_ij[2] = x_j[2] - x_i[2];
			float_value_t dx_ij[3];
			dx_ij[0] = dx_j[0] - dx_i[0];
			dx_ij[1] = dx_j[1] - dx_i[1];
			dx_ij[2] = dx_j[2] - dx_i[2];


			float_value_t F_con_ij[3];
			float_value_t F_fric_ij[3];

			if (i != j)
			{
				float_value_t delta_ij = norm(x_ij) - (R_i + R_j);
				float_value_t v_delta_ij = dot(dx_ij, x_ij) / norm(x_ij);
				if (delta_ij < 0)
				{
					k++;
					//printf("%e\n", delta_ij);
					if (abs(delta_ij) > 1e-12)
					{
						//printf("%e\n", delta_ij);
						//delta_ij = -1e-14;
					}
					m_aster = m_i*m_j / (m_i + m_j);
					d = 2 * zeta*sqrt(2 * E_aster*m_aster*sqrt(R_aster))*pow(abs(delta_ij), 0.25);
					F_con_ij[0] =  -4.0 / 3 * sqrt(R_aster)*E_aster*pow(abs(delta_ij), 1.5)*dctn(x_ij, 0) + d*v_delta_ij*dctn(x_ij, 0);
					F_con_ij[1] = -4.0 / 3 * sqrt(R_aster)*E_aster*pow(abs(delta_ij), 1.5)*dctn(x_ij, 1) + d*v_delta_ij*dctn(x_ij, 1);
					F_con_ij[2] = -4.0 / 3 * sqrt(R_aster)*E_aster*pow(abs(delta_ij), 1.5)*dctn(x_ij, 2) + d*v_delta_ij*dctn(x_ij, 2);
					float_value_t vtij[3];
					vtij[0] = dx_ij[0] - (dot(dx_j, x_ij) - dot(dx_i, x_ij)) / ((R_i + R_j) + delta_ij) / ((R_i + R_j) + delta_ij)*x_ij[0];
					vtij[1] = dx_ij[1] - (dot(dx_j, x_ij) - dot(dx_i, x_ij)) / ((R_i + R_j) + delta_ij) / ((R_i + R_j) + delta_ij)*x_ij[1];
					vtij[2] = dx_ij[2] - (dot(dx_j, x_ij) - dot(dx_i, x_ij)) / ((R_i + R_j) + delta_ij) / ((R_i + R_j) + delta_ij)*x_ij[2];
					F_fric_ij[0] = mu_d*norm(F_con_ij)*dctn(vtij, 0);
					F_fric_ij[1] = mu_d*norm(F_con_ij)*dctn(vtij, 1);
					F_fric_ij[2] = mu_d*norm(F_con_ij)*dctn(vtij, 2);
					F_con[0] = F_con[0] + F_con_ij[0];
					F_con[1] = F_con[1] + F_con_ij[1];
					F_con[2] = F_con[2] + F_con_ij[2];

					F_fric[0] = F_fric[0] + F_fric_ij[0];
					F_fric[1] = F_fric[1] + F_fric_ij[1];
					F_fric[2] = F_fric[2] + F_fric_ij[2];


				}
			}
		}

		float_value_t r = norm(F_con) / (m_i * 9.8);
		if (r>MAX_R)
		{
			F_con[0] = F_con[0] / r * MAX_R;
			F_con[1] = F_con[1] / r * MAX_R;
			F_con[2] = F_con[2] / r * MAX_R;
			F_fric[0] = F_fric[0] / r * MAX_R;
			F_fric[1] = F_fric[1] / r * MAX_R;
			F_fric[2] = F_fric[2] / r * MAX_R;
		}
		part->as[3 * i] = (F_env[0] + F_con[0] + F_fric[0]) / (m_i*C_bar*omega_bar*omega_bar) ;
		part->as[3 * i + 1] = (F_env[1] + F_con[1] + F_fric[1]) / (m_i*C_bar*omega_bar*omega_bar);
		part->as[3 * i + 2] = (F_env[2] + F_con[2] + F_fric[2]+ m_i*g_grav) / (m_i*C_bar*omega_bar*omega_bar);
		part->vs[3 * i] = part->vs[3 * i] + delta_t*part->as[3 * i];
		part->vs[3 * i + 1] = part->vs[3 * i + 1] + delta_t*part->as[3 * i + 1];
		part->vs[3 * i + 2] = part->vs[3 * i + 2] + delta_t*part->as[3 * i + 2];
	}
	memcpy(res,  &x_input[3 * N_PARTICLES], 3 * sizeof(float_value_t)*N_PARTICLES);
	memcpy((&res[3 * N_PARTICLES]),  part->as, 3 * sizeof(float_value_t)*N_PARTICLES);
}
int main(void)
{
	stopwatch<std::milli, float> sw;
	particles_t* particles = (particles_t *)malloc(sizeof(particles_t));

	initParticles(particles);

	// simulation of the deposition of powder particles
	//printf("%d\n", 2<<int(log10f(600 / 2.0) / log10f(2.0)));

	float_value_t *x = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *y1 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *y2 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *y3 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *y4 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *fy1 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *fy2 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *fy3 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	float_value_t *fy4_1 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	//float_value_t *fy2_2 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	//float_value_t *fy3_2 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);




	std::string coord0_file = "xyzs0.txt";
	std::string vcord_file = "vxyzs0.txt";
	writeFile(particles->xyzs, particles->n_particles, 3, coord0_file);
	writeFile(&(x[3 * N_PARTICLES]), particles->n_particles, 3, vcord_file);
	for (int k = 0; k < N_PARTICLES; k++)
	{
		particles->xyzs[3 * k] = particles->xyzs[3 * k] / C_bar;
		particles->xyzs[3 * k + 1] = particles->xyzs[3 * k + 1] / C_bar;
		particles->xyzs[3 * k + 2] = particles->xyzs[3 * k + 2] / C_bar;
	}

	sw.start();
	memcpy(x, particles->xyzs, 3 * sizeof(float_value_t)*N_PARTICLES);
	memcpy(&(x[3 * N_PARTICLES]), particles->vs, 3 * sizeof(float_value_t)*N_PARTICLES);


	double current_time = 0.0f;
	int step = 0;
	double T_total = TIME*omega_bar;
	printf("%e, %e\n", T_total/ TIME_STEP, T_total / TIME_STEP/omega_bar);
//	getchar();
	int count = 0;
	while (current_time < T_total)
	{
		count++;
		for (int j = 0; j < 6 * N_PARTICLES; j++)
			y1[j] = x[j];
		f(fy1, particles, y1, 0);
		for (int j = 0; j < 6 * N_PARTICLES; j++)
		{
			y2[j] = x[j] + 0.5 * TIME_STEP*fy1[j];
		}
		f(fy2, particles, y2,  0.5 * TIME_STEP);
		for (int j = 0; j < 6 * N_PARTICLES; j++)
		{
			y3[j] = x[j] + 0.5 * TIME_STEP*fy2[j];
		}
		f(fy3, particles, y3,  0.5 * TIME_STEP);
		for (int j = 0; j < 6 * N_PARTICLES; j++)
		{
			y4[j] = x[j] +  TIME_STEP*fy3[j];
		}

		//f(fy2_2, y2, t + delta_t / 2);
		//f(fy3_2, y2, t + delta_t / 2);
		f(fy4_1, particles, y4, TIME_STEP);
		for (int j = 0; j < 6 * N_PARTICLES; j++)
		{
			x[j] = x[j] + 1.0/6 * TIME_STEP*(fy1[j] + 2 * fy2[j] + 2 * fy3[j] + fy4_1[j]);
		}

		//updatePositions(particles, current_time, TIME_STEP);
		current_time += TIME_STEP;
		if (current_time / T_total > step*0.01)
		{
			float_value_t min_z = boundary_up,e_t=0.0;
			for (int k = 0; k < N_PARTICLES; k++)
			{
				if (x[3 * k + 2] < min_z)
					min_z = x[3 * k + 2];
				e_t += 0.5*particles->ms[k] * (x[3 * k + 0 + 3 * N_PARTICLES] * x[3 * k + 0 + 3 * N_PARTICLES] + x[3 * k + 1 + 3 * N_PARTICLES] * x[3 * k + 1 + 3 * N_PARTICLES] + x[3 * k + 2 + 3 * N_PARTICLES] * x[3 * k + 2 + 3 * N_PARTICLES]) + particles->ms[k] * 9.8*x[3 * k + 2];
				/*if (abs(particles->xyzs[3 * k]) > 0.01 || abs(particles->xyzs[3 * k + 1]) > 0.01 || abs(particles->xyzs[3 * k + 2]) > 0.01)
				{
				printf("blow up!!! time: %f\n", current_time);
				b_blowup = true;
				break;
				}*/
				if (isnan(x[3 * k + 0]) || isnan(x[3 * k + 1]) || isnan(x[3 * k + 2]))
				{
					printf("blowup!\n");
					getchar();
					exit(1);
				}
				//printf("%d paritcle: z = %f, r = %f, v = [%e %e %e]\n", k, x[3 * k + 2] * C_bar, particles->rs[k] , x[3 * k + 0 + 3 * N_PARTICLES] * omega_bar*C_bar, x[3 * k + 1 + 3 * N_PARTICLES] * omega_bar*C_bar, x[3 * k + 2+3* N_PARTICLES] * omega_bar*C_bar);
			}
			printf("%d%% has been completed. total energy = %e\n", step, e_t);
			step++;

			memcpy(particles->xyzs, x, 3 * sizeof(float_value_t)*N_PARTICLES);
			memcpy((particles->vs), &(x[3 * N_PARTICLES]), 3 * sizeof(float_value_t)*N_PARTICLES);

			// write to file
			std::string coord_file = "xyzs.txt";
			std::string radius_file = "radius.txt";
			for (int k = 0; k < N_PARTICLES; k++)
			{
				particles->xyzs[3 * k] = particles->xyzs[3 * k] * C_bar;
				particles->xyzs[3 * k + 1] = particles->xyzs[3 * k + 1] * C_bar;
				particles->xyzs[3 * k + 2] = particles->xyzs[3 * k + 2] * C_bar;
				particles->vs[3 * k] = particles->vs[3 * k] * C_bar*omega_bar;
				particles->vs[3 * k + 1] = particles->vs[3 * k + 1] * C_bar*omega_bar;
				particles->vs[3 * k + 2] = particles->vs[3 * k + 2] * C_bar*omega_bar;
				particles->rs[k] = particles->rs[k];
			}
			std::string vs_file = "vxyz.txt";
			writeFile(particles->xyzs, particles->n_particles, 3, coord_file);
			writeFile(particles->rs, particles->n_particles, 1, radius_file);
			writeFile(particles->vs, particles->n_particles, 3, vs_file);
		}
	}
	sw.stop();
	printf("total time: %f\n", sw.count());
	// simulation of temperature evolution
	//
	// free memory

	freeParticles(particles);

	free(x);
	free(y1);
	free(y2);
	free(y3);
	free(y4);
	free(fy1);
	free(fy2);
	free(fy3);
	free(fy4_1);

	return 0;
}
void freeParticles(particles_t * particles)
{
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
	free(particles->E);
	particles->E = NULL;
	free(particles->mu);
	particles->mu = NULL;
	free(particles);
	particles = NULL;
}

void writeFile(float_value_t* array, int row, int col, std::string filename)
{
	std::ofstream myfile;
	myfile.open(filename);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			myfile << array[i * col + j] << " ";
		}
		myfile << "\n";
	}
	myfile.close();
}

void initParticles( particles_t * particles)
{
	particles->n_particles = N_PARTICLES;
	particles->xyzs = initFloatArray(particles->n_particles * 3);
	particles->rs = initFloatArray(particles->n_particles);
	particles->ms = initFloatArray(particles->n_particles);
	particles->fs = initFloatArray(particles->n_particles * 3);
	particles->vs = initFloatArray(particles->n_particles * 3); // zero init
	for (int i = 0; i < N_PARTICLES; i++)
	{
		(particles->vs)[i * 3] = 0;
		(particles->vs)[i * 3 + 1] = 0;
		(particles->vs)[i * 3 + 2] = v0;
	}
	particles->as = initFloatArray(particles->n_particles * 3);
	particles->E = initFloatArray(particles->n_particles * 3);
	particles->mu = initFloatArray(particles->n_particles * 3);
	float_value_t mu = MU;
	float_value_t sigma = SIGMA;
	float_value_t min_height = MIN_HEIGHT;
	float_value_t max_height = MAX_HEIGHT;
	bool find_particle = false;
	float_value_t x, y, z, radius;
	srand(time(NULL));
	std::default_random_engine generator(rand());
	std::normal_distribution<float_value_t> radius_distribution(mu, sigma);
	std::uniform_real_distribution<float_value_t> x_distribution(0.0, length);
	std::uniform_real_distribution<float_value_t> y_distribution(0.0, width);
	std::uniform_real_distribution<float_value_t> z_distribution(min_height, max_height);
	for (int i = 0; i < particles->n_particles; i++)
	{
		while (!find_particle)
		{
			while (true)
			{
				x = x_distribution(generator);
				y = y_distribution(generator);
				z = z_distribution(generator);
				radius = radius_distribution(generator);
				if (x + radius <= length && x - radius >= 0 && y + radius <= width && y - radius >= 0
            && radius >= 0.000010 && radius <= 0.000050) break;
			}

			find_particle = true;
			for (int j = 0; j < i; j++)
			{
				float_value_t dist = 0.0f;
				dist += (x - particles->xyzs[j * 3])*(x - particles->xyzs[j * 3]);
				dist += (y - particles->xyzs[j * 3 + 1])*(y - particles->xyzs[j * 3 + 1]);
				dist += (z - particles->xyzs[j * 3 + 2])*(z - particles->xyzs[j * 3 + 2]);
				dist = sqrt(dist);
				if (dist < radius + particles->rs[j])
				{
					find_particle = false;
					break;
				}
			}
		}

		find_particle = false;
		particles->xyzs[i * 3 + 0] = x;
		particles->xyzs[i * 3 + 1] = y;
		particles->xyzs[i * 3 + 2] = z;
		particles->rs[i] = radius;
		particles->ms[i] = 4.0/3.0*M_PI*radius*radius*radius*rho;
		particles->E[i] = E_modulus;
		particles->mu[i] = nu;
		// if(i == 0) printf("%0.15f\n", particles->ms[i]);
	}

}

float_value_t* initFloatArray(int num_elements)
{
	float_value_t * array = NULL;
	array = (float_value_t *)malloc(sizeof(float_value_t) * num_elements);
	if (array == NULL)
	{
		fprintf(stderr, "malloc fails\n");
		exit(1);
	}
	else
	{
		return array;
	}
}

