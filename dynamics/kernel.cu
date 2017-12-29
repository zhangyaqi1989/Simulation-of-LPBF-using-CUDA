
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
// deposition_sim.cpp : Defines the entry point for the console application.
//


#include<iostream>
#include<fstream>
#include<stdio.h>
#include<random>
#include<math.h>
#include <time.h>

#define C_bar (0.000027)
#define boundary_left (0)
//length
#define boundary_right (0.001)
//width
#define boundary_rear (0.0005)
#define boundary_front (0)
#define boundary_up (0.003)
#define boundary_bottom (0.000)

#define N_PARTICLES (512)

//initial height 
#define MAX_HEIGHT (0.00020)
#define MIN_HEIGHT (0.00004)

//mean radius
#define MU (0.000027)
#define BIND (0.000010)

#define zeta (0.1) 
#define mu_d (0.1)
#define mu_f (2.1e-3)
#define g_grav (-9.8)
#define rho (7800)

#define E_modulus (193e9)
#define nu (0.26)

#define length boundary_right
#define width boundary_rear

#define E_aster (E_modulus/(2-2*nu*nu))

#define norm_x_ij sqrt((x_j[0]-x_i[0])*(x_j[0]-x_i[0])+(x_j[1]-x_i[1])*(x_j[1]-x_i[1])+(x_j[2]-x_i[2])*(x_j[2]-x_i[2])) 
#define R_aster Rs[i]*Rs[j]/(Rs[i]+Rs[j])
//#define E_aster Es[i]*Es[j]/(Es[i]*(1-Mus[i]*Mus[i])+Es[j]*(1-Mus[j]*Mus[j]))
#define norm_dxij sqrt((dx_j[0]-dx_i[0])*(dx_j[0]-dx_i[0])+(dx_j[1]-dx_i[1])*(dx_j[1]-dx_i[1])+(dx_j[2]-dx_i[2])*(dx_j[2]-dx_i[2]))

#define C_factor 0.0055  
//sqrt(1/(4/3*rho*pi))
#define m_bar (4.0/3*M_PI*C_bar*C_bar*C_bar)
#define omega_bar (sqrt(E_modulus/(2-2*nu*nu)*C_bar/m_bar))
#define MAX_R 10

#define v0 (0)
#define TIME_STEP (5e-8*omega_bar) //maximal limit 0.0000001
#define TIME 0.1


#define BLOCKDIM 128
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
	// float_value_t* temperatures;
} particles_t;


void readParticles(particles_t * particles);
void writeFile(float_value_t* array, int row, int col, std::string filename);
void initParticles(particles_t * particles);
float_value_t* initFloatArray(int num_elements);
void freeParticles(particles_t * particles, particles_t * dp);
// void updateAcceleration(particles_t * particles);
void copyParticles(particles_t * particles, particles_t * dp)
{
}
__device__ float_value_t norm(float_value_t* r)
{
	return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}
__device__ float_value_t dot(float_value_t* r, float_value_t* r1)
{
	return (r[0] * r1[0] + r[1] * r1[1] + r[2] * r1[2]);
}
__device__ float_value_t dot3(float_value_t* r, float_value_t* r1, float_value_t* r2)
{
	return (r[0] * (r2[0]-r1[0]) + r[1] * (r2[1]-r1[1]) + r[2] * (r2[2]- r1[2]));
}
__global__ void vector_add(float_value_t* A, float_value_t *B, float_value_t *C, float_value_t factor, int N)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N)
	{
		C[i] = A[i] + factor*B[i];
	}
}
__global__ void vector_add3(float_value_t* x, float_value_t *fy2, float_value_t *fy3, float_value_t *fy4_1, int N)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N)
	{
		x[i] = x[i] + 2.0 / 6 * TIME_STEP*fy2[i] + 2.0 / 6 * TIME_STEP*fy3[i] + 1.0 / 6 * TIME_STEP*fy4_1[i];
	}
}
__global__ void vel_checker(float_value_t* x, float_value_t* y, int N)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N)
	{
		x[3 * i] = y[3 * i];
		x[3 * i + 1] = y[3 * i + 1];
		x[3 * i + 2] = y[3 * i + 2];

		float_value_t vx = x[3 * i + 3 * N_PARTICLES];
		float_value_t vy = x[3 * i + 1 + 3 * N_PARTICLES];
		float_value_t vz = x[3 * i + 2 + 3 * N_PARTICLES];
		float_value_t vn = 1 / 2 * (vx*vx + vy*vy + vz*vz)*omega_bar*C_bar*omega_bar*C_bar +(-g_grav)*x[3 * i + 2] * C_bar;

		float_value_t vx1 = y[3 * i + 3 * N_PARTICLES];
		float_value_t vy1 = y[3 * i + 1 + 3 * N_PARTICLES];
		float_value_t vz1 = y[3 * i + 2 + 3 * N_PARTICLES];
		float_value_t vn1 = 1 / 2 *(vx1*vx1 + vy1*vy1 + vz1*vz1)*omega_bar*C_bar *omega_bar*C_bar +(-g_grav)*y[3 * i + 2] * C_bar;
		if (abs(vn1) >abs(vn) && vn1>1e-6 )
		{
			x[3 * i + 3 * N_PARTICLES] = vx1*vn / vn1;
			x[3 * i + 1 + 3 * N_PARTICLES] = vy1*vn / vn1;
			x[3 * i + 2 + 3 * N_PARTICLES] = vz1*vn / vn1;
		}
		else
		{
			x[3 * i + 3 * N_PARTICLES] = vx1;
			x[3 * i + 1 + 3 * N_PARTICLES] = vy1;
			x[3 * i + 2 + 3 * N_PARTICLES] = vz1;
		}
	}
}
__global__ void f(float_value_t* res, float_value_t* res1, float_value_t coeff, float_value_t *Es, float_value_t *Rs, float_value_t *Ms, float_value_t *Mus, float_value_t *v3s, float_value_t* x_input, float_value_t delta_t)
{
	int j = threadIdx.x;//blockDim.x*blockIdx.x+
	int i = blockIdx.x;
	__shared__  float_value_t s_fcon[3 * N_PARTICLES];
	__shared__  float_value_t s_ffric[3 * N_PARTICLES];

	//float_value_t* x_i = &(part->xyzs[3 * i]);
	//float_value_t* x_i = &(x_input[3 * i]);
	float_value_t x_i[3];
	x_i[0] = x_input[3 * i]*C_bar;
	x_i[1] = x_input[3 * i + 1] * C_bar;
	x_i[2] = x_input[3 * i + 2] * C_bar;
	//float_value_t* dx_i = &(part->vs[3 * i]);
	//float_value_t* dx_i = &(x_input[3 * i + 3 * N_PARTICLES]);
	float_value_t dx_i[3];
	dx_i[0] = x_input[3 * i + 3 * N_PARTICLES] *omega_bar* C_bar;
	dx_i[1] = x_input[3 * i + 1 + 3 * N_PARTICLES] * omega_bar* C_bar;
	dx_i[2] = x_input[3 * i + 2 + 3 * N_PARTICLES] * omega_bar* C_bar;

	float_value_t F_env[3] = { 0 };
	float_value_t F_con[3] = { 0 };
	float_value_t F_fric[3] = { 0 };

	float_value_t F_con_iw[3] = { 0 };
	float_value_t F_fric_iw[3] = { 0 };


	F_env[0] = -6 * M_PI*mu_f*Rs[i] *dx_i[0];
	F_env[1] = -6 * M_PI*mu_f*Rs[i] *dx_i[1];
	F_env[2] = -6 * M_PI*mu_f*Rs[i] *dx_i[2];

	if (j == N_PARTICLES)
	{
		float_value_t delta_iw = 0.0;
		if (x_i[0] - Rs[i] <= boundary_left)
		{
			delta_iw = x_i[0] - boundary_left - (Rs[i]);
			F_con_iw[0] = 4.0 / 3 * sqrt(Rs[i]) * E_aster*pow(abs(delta_iw), 1.5) + 2 * zeta*sqrt(2 * E_aster* Ms[i] *sqrt(Rs[i]))*pow(abs(delta_iw), (0.25))* (-dx_i[0]);
			if (F_con_iw[0] < 0)
				F_con_iw[0] = 0;
			F_con_iw[1] = 0;
			F_con_iw[2] = 0;
		}
		if (x_i[0] + Rs[i] >= boundary_right)
		{
			delta_iw = x_i[0] - boundary_right + (Rs[i]);
			F_con_iw[0] = -4.0 / 3 * sqrt(Rs[i]) * E_aster*pow(abs(delta_iw), 1.5) - 2 * zeta*sqrt(2 * E_aster* Ms[i] *sqrt(Rs[i]))*pow(abs(delta_iw), (0.25))* dx_i[0] ;
			if (F_con_iw[0] > 0)
				F_con_iw[0] = 0;
			F_con_iw[1] = 0;
			F_con_iw[2] = 0;
		}
		if (x_i[1] - Rs[i] <= boundary_front)
		{
			delta_iw = x_i[1] - boundary_front - (Rs[i]);
			F_con_iw[0] = 0;
			F_con_iw[1] = 4.0 / 3 * sqrt(Rs[i]) * E_aster*pow(abs(delta_iw), 1.5) + 2 * zeta*sqrt(2 * E_aster* Ms[i] *sqrt(Rs[i]))*pow(abs(delta_iw), (0.25))* (-dx_i[1]);
			F_con_iw[2] = 0;
			if (F_con_iw[1] < 0)
				F_con_iw[1] = 0;
		}
		if (x_i[1] + Rs[i] >= boundary_rear)
		{
			delta_iw = x_i[1] - boundary_rear + (Rs[i]);
			F_con_iw[0] = 0;
			F_con_iw[1] = -4.0 / 3 * sqrt(Rs[i]) * E_aster*pow(abs(delta_iw), 1.5) - 2 * zeta*sqrt(2 * E_aster* Ms[i] *sqrt(Rs[i]))*pow(abs(delta_iw), (0.25))* dx_i[1] ;
			F_con_iw[2] = 0;
			if (F_con_iw[1] > 0)
				F_con_iw[1] = 0;
		}
		if (x_i[2] - Rs[i] <= boundary_bottom)
		{
			delta_iw = x_i[2] - boundary_bottom - (Rs[i]);
			F_con_iw[0] = 0;
			F_con_iw[1] = 0;
			F_con_iw[2] = 4.0 / 3 * sqrt(Rs[i]) * E_aster*pow(abs(delta_iw), 1.5) + 2 * zeta*sqrt(2 * E_aster* Ms[i] *sqrt(Rs[i]))*pow(abs(delta_iw), (0.25))* (-dx_i[2]);

			if (F_con_iw[2] < 0)
				F_con_iw[2] = 0;

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
	}
	if (j < N_PARTICLES)
	{
		if (i == j)
		{
			s_fcon[3 * j + 0] = 0;
			s_fcon[3 * j + 1] = 0;
			s_fcon[3 * j + 2] = 0;
			s_ffric[3 * j + 0] = 0;
			s_ffric[3 * j + 1] = 0;
			s_ffric[3 * j + 2] = 0;
		}
		else
		{
			float_value_t x_j[3];
			x_j[0] = x_input[3 * j] * C_bar;
			x_j[1] = x_input[3 * j + 1] * C_bar;
			x_j[2] = x_input[3 * j + 2] * C_bar;
			//float_value_t* x_j = &(x_input[3 * j]);
			//float_value_t* dx_j = &(x_input[3 * j + 3 * N_PARTICLES]);
			float_value_t dx_j[3];
			dx_j[0] = x_input[3 * j + 3 * N_PARTICLES] * omega_bar* C_bar;
			dx_j[1] = x_input[3 * j + 1 + 3 * N_PARTICLES] * omega_bar* C_bar;
			dx_j[2] = x_input[3 * j + 2 + 3 * N_PARTICLES] * omega_bar* C_bar;



			float_value_t delta_ij = norm_x_ij - (Rs[i] + Rs[j]);
			
			if (delta_ij < 0)
			{
				//if (abs(delta_ij) > 0.0001)
				//	delta_ij = 0;
				float_value_t v_delta_ij = ((x_j[0] - x_i[0])*(dx_j[0] - dx_i[0])+(x_j[1] - x_i[1])*(dx_j[1] - dx_i[1])+(x_j[2] - x_i[2])*(dx_j[2] - dx_i[2])) / (delta_ij + (Rs[i] + Rs[j]));

				s_fcon[3 * j] = -4.0 / 3 * sqrt(R_aster)*E_aster*pow(abs(delta_ij), 1.5)*(x_j[0] - x_i[0]) / (delta_ij+ (Rs[i] + Rs[j])) -2 * zeta*sqrt(2 * E_aster* Ms[i] * Ms[j] / (Ms[i] + Ms[j])*sqrt(R_aster))*pow(abs(delta_ij), 0.25)*v_delta_ij*(x_j[0] - x_i[0]) / (delta_ij + (Rs[i] + Rs[j]));
				if (s_fcon[3 * j] * (x_j[0] - x_i[0]) > 0)
					s_fcon[3 * j] = 0;
				s_fcon[3 * j + 1] = -4.0 / 3 * sqrt(R_aster)*E_aster*pow(abs(delta_ij), 1.5)*(x_j[1] - x_i[1]) / (delta_ij + (Rs[i] + Rs[j])) -2 * zeta*sqrt(2 * E_aster* Ms[i] * Ms[j] / (Ms[i] + Ms[j])*sqrt(R_aster))*pow(abs(delta_ij), 0.25)*v_delta_ij*(x_j[1] - x_i[1]) / (delta_ij + (Rs[i] + Rs[j]));
				if (s_fcon[3 * j + 1] * (x_j[1] - x_i[1]) > 0)
					s_fcon[3 * j + 1] = 0;
				s_fcon[3 * j + 2] = -4.0 / 3 * sqrt(R_aster)*E_aster*pow(abs(delta_ij), 1.5)*(x_j[2] - x_i[2]) / (delta_ij + (Rs[i] + Rs[j])) -2 * zeta*sqrt(2 * E_aster*Ms[i] * Ms[j] / (Ms[i] + Ms[j])*sqrt(R_aster))*pow(abs(delta_ij), 0.25)*v_delta_ij*(x_j[2] - x_i[2]) / (delta_ij + (Rs[i] + Rs[j]));
				if (s_fcon[3 * j + 2] * (x_j[2] - x_i[2]) > 0)
					s_fcon[3 * j + 2] = 0;
				float_value_t vtij[3];
				vtij[0] = x_j[0] - x_i[0] - ((dot3(dx_j, x_i, x_j) - dot3(dx_i, x_i, x_j)) / (delta_ij + (Rs[i] + Rs[j]))) *((x_j[0] - x_i[0]) / (delta_ij + (Rs[i] + Rs[j])));
				vtij[1] = x_j[1] - x_i[1] - ((dot3(dx_j, x_i, x_j) - dot3(dx_i, x_i, x_j)) / (delta_ij + (Rs[i] + Rs[j]))) *((x_j[1] - x_i[1]) / (delta_ij + (Rs[i] + Rs[j])));
				vtij[2] = x_j[2] - x_i[2] - ((dot3(dx_j, x_i, x_j) - dot3(dx_i, x_i, x_j)) / (delta_ij + (Rs[i] + Rs[j]))) *((x_j[2] - x_i[2]) / (delta_ij + (Rs[i] + Rs[j])));
				s_ffric[3 * j + 0] = mu_d*norm((float_value_t*)&(s_fcon[3 * j]))*(vtij[0] - vtij[0]) / sqrt(vtij[0] * vtij[0] + vtij[1] * vtij[1] + vtij[2] * vtij[2]);
				s_ffric[3 * j + 1] = mu_d*norm((float_value_t*)&(s_fcon[3 * j]))*(vtij[1] - vtij[1]) / sqrt(vtij[0] * vtij[0] + vtij[1] * vtij[1] + vtij[2] * vtij[2]);
				s_ffric[3 * j + 2] = mu_d*norm((float_value_t*)&(s_fcon[3 * j]))*(vtij[2] - vtij[2]) / sqrt(vtij[0] * vtij[0] + vtij[1] * vtij[1] + vtij[2] * vtij[2]);

			}
			else
			{
				s_fcon[3 * j + 0] = 0;
				s_fcon[3 * j + 1] = 0;
				s_fcon[3 * j + 2] = 0;
				s_ffric[3 * j + 0] = 0;
				s_ffric[3 * j + 1] = 0;
				s_ffric[3 * j + 2] = 0;
			}
		}
	}
	__syncthreads();
	//hard coding here. Using 300 because N_PARTICLES/2==300)
	if (j < N_PARTICLES / 2)
	{
		s_fcon[3 * j + 0] += s_fcon[3 * (j + N_PARTICLES / 2) + 0];
		s_fcon[3 * j + 1] += s_fcon[3 * (j + N_PARTICLES / 2) + 1];
		s_fcon[3 * j + 2] += s_fcon[3 * (j + N_PARTICLES / 2) + 2];
		s_ffric[3 * j + 0] += s_ffric[3 * (j + N_PARTICLES / 2) + 0];
		s_ffric[3 * j + 1] += s_ffric[3 * (j + N_PARTICLES / 2) + 1];
		s_ffric[3 * j + 2] += s_ffric[3 * (j + N_PARTICLES / 2) + 2];
	}
	else if (j<N_PARTICLES)
	{
		s_fcon[3 * j + 0] = 0;
		s_fcon[3 * j + 1] = 0;
		s_fcon[3 * j + 2] = 0;
		s_ffric[3 * j + 0] = 0;
		s_ffric[3 * j + 1] = 0;
		s_ffric[3 * j + 2] = 0;
	}

	__syncthreads(); //
	for (unsigned int s = (1 << int(log10f(N_PARTICLES / 2.0) / log10f(2.0))); s > 0; s >>= 1)
	{
		if (j < s)
		{
			s_fcon[3 * j + 0] += s_fcon[3 * (j + s) + 0];
			s_fcon[3 * j + 1] += s_fcon[3 * (j + s) + 1];
			s_fcon[3 * j + 2] += s_fcon[3 * (j + s) + 2];
			s_ffric[3 * j + 0] += s_ffric[3 * (j + s) + 0];
			s_ffric[3 * j + 1] += s_ffric[3 * (j + s) + 1];
			s_ffric[3 * j + 2] += s_ffric[3 * (j + s) + 2];
		}
		__syncthreads();
	}
	__syncthreads();
	if (j == N_PARTICLES)
	{
		F_con[0] = s_fcon[0] + F_con_iw[0];
		F_con[1] = s_fcon[1] + F_con_iw[1];
		F_con[2] = s_fcon[2] + F_con_iw[2];
		F_fric[0] = s_ffric[0] + F_fric_iw[0];
		F_fric[1] = s_ffric[1] + F_fric_iw[1];
		F_fric[2] = s_ffric[2] + F_fric_iw[2];

		float_value_t r = norm(F_con) /  (Ms[i] * 9.8);
		if (r>MAX_R)
		{
			F_con[0] = F_con[0] / r * MAX_R;
			F_con[1] = F_con[1] / r * MAX_R;
			F_con[2] = F_con[2] / r * MAX_R;
			F_fric[0] = F_fric[0] / r * MAX_R;
			F_fric[1] = F_fric[1] / r * MAX_R;
			F_fric[2] = F_fric[2] / r * MAX_R;
		}
		
		res[3 * i] = v3s[3 * i + 0];
		res[3 * i + 1] = v3s[3 * i + 1];
		res[3 * i + 2] = v3s[3 * i + 2];

		res[3 * i + 3 * N_PARTICLES] = (F_env[0] + F_con[0] + F_fric[0]) / (Ms[i] * C_bar* omega_bar *omega_bar);
		res[3 * i + 1 + 3 * N_PARTICLES] = (F_env[1] + F_con[1] + F_fric[1]) / (Ms[i] * C_bar* omega_bar *omega_bar);
		res[3 * i + 2 + 3 * N_PARTICLES] = (F_env[2] + F_con[2] + F_fric[2]+ g_grav*Ms[i]) / (Ms[i] * C_bar* omega_bar *omega_bar);

		//v3s[3 * i + 0] = v3s[3 * i] + delta_t*res[3 * i + 3 * N_PARTICLES];
		//v3s[3 * i + 1] = v3s[3 * i + 1] + delta_t*res[3 * i + 1 + 3 * N_PARTICLES];
		//v3s[3 * i + 2] = v3s[3 * i + 2] + delta_t*res[3 * i + 2 + 3 * N_PARTICLES];
	}
	__syncthreads();
	if (i < N_PARTICLES)
	{
		if (j < 3)
		{
			res1[3 * i + j] = coeff*res[3 * i + j] + x_input[3 * i + j];
		}
		else if(j>=3 && j<6)
			res1[3 * i + j-3 + 3 * N_PARTICLES] = coeff*res[3 * i + j - 3 + 3 * N_PARTICLES] + x_input[3 * i + j - 3 + 3 * N_PARTICLES];
	}
}
int main(void)
{
	particles_t* particles = (particles_t *)malloc(sizeof(particles_t));

	initParticles(particles);
	//printf("vel_max: %e\n",g_grav / (6 * M_PI*mu_f*MU/C_bar/C_bar / (4.0 / 3 * M_PI*MU*MU*MU)));
	// simulation of the deposition of powder particles

	float_value_t *x, *x1, *y2, *y3, *y4, *fy1, *fy2, *fy3, *fy4_1;
	particles_t dp[1];

	cudaMalloc((void**)&(dp->xyzs), 3 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->rs), 1 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->ms), 1 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->fs), 3 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->vs), 3 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->as), 3 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->E), 1 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->mu), 1 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&(dp->vs), 3 * sizeof(float_value_t)*N_PARTICLES);
	cudaMemcpy((dp->xyzs), (particles->xyzs), 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->vs), (particles->vs), 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->as), (particles->as), 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->fs), (particles->fs), 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->rs), (particles->rs), 1 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->E), (particles->E), 1 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->ms), (particles->ms), 1 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->mu), (particles->mu), 1 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy((dp->vs), (particles->vs), 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);


	cudaMalloc((void**)&x, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&x1, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&y2, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&y3, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&y3, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&y4, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&fy1, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&fy2, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&fy3, 6 * sizeof(float_value_t)*N_PARTICLES);
	cudaMalloc((void**)&fy4_1, 6 * sizeof(float_value_t)*N_PARTICLES);
	//float_value_t *fy2_2 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);
	//float_value_t *fy3_2 = (float_value_t*)malloc(6 * sizeof(float_value_t)*N_PARTICLES);

	std::string coord0_file = "xyzs0.txt";
	writeFile(particles->xyzs, particles->n_particles, 3, coord0_file);


	cudaEvent_t startEvent_inc, stopEvent_inc;
	cudaEventCreate(&startEvent_inc);
	cudaEventCreate(&stopEvent_inc);
	for (int k = 0; k < N_PARTICLES; k++)
	{
		particles->xyzs[3 * k] = particles->xyzs[3 * k] / C_bar;
		particles->xyzs[3 * k + 1] = particles->xyzs[3 * k+1] / C_bar;
		particles->xyzs[3 * k + 2] = particles->xyzs[3 * k + 2] / C_bar;
	}

	cudaEventRecord(startEvent_inc, 0); // starting timing for inclusive
	cudaMemcpy(x, particles->xyzs, 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy(&(x[3 * N_PARTICLES]), particles->vs, 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyHostToDevice);

	float elapsedTime_inc;

	printf("timestep: %e\nomega: %e\n", TIME_STEP, omega_bar);
	double current_time = 0.0f;
	int step = 0;
	double T_total = TIME*omega_bar;
	bool b_blowup = false;
	while (current_time < T_total &&!b_blowup)
	{

		f <<< N_PARTICLES, N_PARTICLES + 1 >>> (fy1, y2, 0.5*TIME_STEP, (dp->E), (dp->rs), (dp->ms), (dp->mu), &(x[3 * N_PARTICLES]), x, 0);
		//vector_add <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>> (x, fy1, y2, 0.5*TIME_STEP, 6 * N_PARTICLES);
		cudaError_t e = cudaGetLastError();
		if (e != cudaSuccess) {
			printf("Cuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e)); getchar();
			exit(0);
		}

		f <<< N_PARTICLES, N_PARTICLES+1 >>> (fy2, y3, 0.5*TIME_STEP, (dp->E), (dp->rs), (dp->ms), (dp->mu), &(x[3 * N_PARTICLES]), y2, 0.5 * TIME_STEP);
		//vector_add <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>> (x, fy2, y3, 0.5*TIME_STEP, 6 * N_PARTICLES);

		f <<< N_PARTICLES, N_PARTICLES + 1 >>> (fy3, y4, TIME_STEP, (dp->E), (dp->rs), (dp->ms), (dp->mu), &(x[3 * N_PARTICLES]), y3, 0.5 * TIME_STEP);
		//vector_add <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>> (x, fy3, y4, TIME_STEP, 6 * N_PARTICLES);


		//f(fy2_2, y2, t + delta_t / 2);
		//f(fy3_2, y2, t + delta_t / 2);
		f <<< N_PARTICLES, N_PARTICLES + 1 >>> (fy4_1, x1, 1.0 / 6 * TIME_STEP, (dp->E), (dp->rs), (dp->ms), (dp->mu), &(x[3 * N_PARTICLES]), y4, TIME_STEP);

		//vector_add <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>>(x, fy1, x, 1.0 / 6 * TIME_STEP, 6 * N_PARTICLES);
		//vector_add <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>>(x, fy2, x, 2.0 / 6 * TIME_STEP, 6 * N_PARTICLES);
		//vector_add <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>>(x, fy3, x, 2.0 / 6 * TIME_STEP, 6 * N_PARTICLES);
		//vector_add <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>>(x, fy4_1, x, 1.0 / 6 * TIME_STEP, 6 * N_PARTICLES);
		vector_add3 <<< (6 * N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>> (x, fy2, fy3, fy4_1, 6 * N_PARTICLES);
		//vel_checker <<< (N_PARTICLES + BLOCKDIM - 1) / BLOCKDIM, BLOCKDIM >>> (x, x1,  N_PARTICLES);
		//updatePositions(particles, current_time, TIME_STEP);
		current_time += TIME_STEP;
		
		if (current_time / T_total > step*0.01)
		{
			cudaMemcpy(particles->xyzs, x, 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyDeviceToHost);
			cudaMemcpy(particles->vs, &(x[3 * N_PARTICLES]), 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyDeviceToHost);

			float_value_t min_z = boundary_up, min_vz;
			for (int k = 0; k < N_PARTICLES; k++)
			{
				particles->xyzs[3 * k] = particles->xyzs[3 * k] * C_bar;
				particles->xyzs[3 * k + 1] = particles->xyzs[3 * k + 1] * C_bar;
				particles->xyzs[3 * k + 2] = particles->xyzs[3 * k + 2] * C_bar;
				particles->vs[3 * k] = particles->vs[3 * k] * C_bar*omega_bar;
				particles->vs[3 * k + 1] = particles->vs[3 * k + 1] * C_bar*omega_bar;
				particles->vs[3 * k + 2] = particles->vs[3 * k + 2] * C_bar*omega_bar;
				if (particles->xyzs[3 * k + 2] < min_z)
				{
					min_z = particles->xyzs[3 * k + 2];
					min_vz = particles->vs[3 * k + 2];
				}
				
				if (abs(particles->xyzs[3 * k]) > 0.01 || abs(particles->xyzs[3 * k + 1]) > 0.01 || abs(particles->xyzs[3 * k + 2]) > 0.01)
				{
					printf("blow up!!! time: %f\n", current_time);
					b_blowup = true;
					break;
				}
				
				if (isnan(particles->xyzs[3 * k]) || isnan(particles->xyzs[3 * k]) || isnan(particles->xyzs[3 * k]))
				{
					printf("blowup!\n");
					getchar();
					exit(1);
				}
				printf("%d paritcle: z = %e, r = %e, v = [%e %e %e]\n", k, particles->xyzs[3 * k + 2], particles->rs[k] , particles->vs[3 * k + 0], particles->vs[3 * k + 1], particles->vs[3 * k + 2]);
			}
			std::string coord_file = "xyzs.txt";
			std::string v_file = "vxyzs.txt";
			std::string radius_file = "radius.txt";
			writeFile(particles->xyzs, particles->n_particles, 3, coord_file);
			writeFile(particles->rs, particles->n_particles, 1, radius_file);
			writeFile(particles->vs, particles->n_particles, 3, v_file);
			printf("%d%% has been completed. min z = %e\t vz = %e\n", step,min_z, min_vz);
			step++;
		}
	}

	cudaMemcpy(particles->xyzs, x, 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyDeviceToHost);
	cudaMemcpy(particles->vs, &(x[3 * N_PARTICLES]), 3 * sizeof(float_value_t)*N_PARTICLES, cudaMemcpyDeviceToHost);

	cudaEventRecord(stopEvent_inc, 0);
	cudaEventSynchronize(stopEvent_inc);
	cudaEventElapsedTime(&elapsedTime_inc, startEvent_inc, stopEvent_inc);
	printf("total time: %f\n", elapsedTime_inc);

	// write to file

	// simulation of temperature evolution
	//
	// free memory

	freeParticles(particles, dp);

	cudaFree(x);
	cudaFree(x1);
	cudaFree(y2);
	cudaFree(y3);
	cudaFree(y4);
	cudaFree(fy1);
	cudaFree(fy2);
	cudaFree(fy3);
	cudaFree(fy4_1);

	return 0;
}
void freeParticles(particles_t * particles, particles_t * dp)
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

	cudaFree(dp->xyzs);
	cudaFree(dp->rs);
	cudaFree(dp->ms);
	cudaFree(dp->fs);
	cudaFree(dp->vs);
	cudaFree(dp->as);
	cudaFree(dp->E);
	cudaFree(dp->mu);
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

void initParticles(particles_t * particles)
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
		(particles->vs)[i * 3+1] = 0;
		(particles->vs)[i * 3 +2] = v0/omega_bar/C_bar;
	}
	particles->as = initFloatArray(particles->n_particles * 3);
	particles->E = initFloatArray(particles->n_particles * 3);
	particles->mu = initFloatArray(particles->n_particles * 3);
	bool find_particle = false;
	float_value_t x, y, z, radius;
	srand(time(NULL));
	for (int i = 0; i < particles->n_particles; i++)
	{
		while (!find_particle)
		{
			while (true)
			{
				x = (float_value_t)rand()/RAND_MAX*length;
				y = (float_value_t)rand() / RAND_MAX*width;
				z = (float_value_t)rand() / RAND_MAX*(MAX_HEIGHT - MIN_HEIGHT)+ MIN_HEIGHT;
				radius = (float_value_t)rand() / RAND_MAX*2* BIND+(MU- BIND);
				if (x + radius <= length && x - radius >= 0 && y + radius <= width && y - radius >= 0) break;
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
		particles->ms[i] = 4.0 / 3.0*M_PI*radius*radius*radius*rho;
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


void readParticles(particles_t * particles)
{
	FILE* fr = fopen("radius.txt", "r");
	FILE* fv = fopen("vxyzs.txt", "r");
	FILE* fx = fopen("xyzs.txt", "r");
	particles->n_particles = N_PARTICLES;
	particles->xyzs = initFloatArray(particles->n_particles * 3);
	particles->rs = initFloatArray(particles->n_particles);
	particles->ms = initFloatArray(particles->n_particles);
	particles->fs = initFloatArray(particles->n_particles * 3);
	particles->vs = initFloatArray(particles->n_particles * 3); // zero init
	particles->as = initFloatArray(particles->n_particles * 3);
	particles->E = initFloatArray(particles->n_particles * 3);
	particles->mu = initFloatArray(particles->n_particles * 3);
	for (int i = 0; i < N_PARTICLES; i++)
	{
		fscanf(fv, "%lf", &(particles->vs[i * 3]));
		fscanf(fv, "%lf", &(particles->vs[i * 3 + 1]));
		fscanf(fv, "%lf", &(particles->vs[i * 3 + 2]));
		(particles->vs)[i * 3] = (particles->vs)[i * 3]/omega_bar/C_bar;
		(particles->vs)[i * 3 + 1] = (particles->vs)[i * 3 + 1] / omega_bar/C_bar;
		(particles->vs)[i * 3 + 2] = (particles->vs)[i * 3 + 2] / omega_bar/C_bar;


		fscanf(fx, "%lf", &(particles->xyzs[i * 3]));
		fscanf(fx, "%lf", &(particles->xyzs[i * 3 + 1]));
		fscanf(fx, "%lf", &(particles->xyzs[i * 3 + 2]));
		(particles->xyzs)[i * 3] = (particles->xyzs)[i * 3] / C_bar;
		(particles->xyzs)[i * 3 + 1] = (particles->xyzs)[i * 3 + 1] / C_bar;
		(particles->xyzs)[i * 3 + 2] = (particles->xyzs)[i * 3 + 2] / C_bar;


		fscanf(fr, "%lf", &(particles->rs[i]));
		particles->ms[i] = 4.0 / 3.0*M_PI*particles->rs[i] * particles->rs[i] * particles->rs[i] *rho;
		particles->E[i] = E_modulus;
		particles->mu[i] = nu;
	}
	fclose(fr);
	fclose(fv);
	fclose(fx);
}
