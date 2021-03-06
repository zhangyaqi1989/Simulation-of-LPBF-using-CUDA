/* written by Yaqi Zhang
 * Dec 2017
 */
#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#define N_PARTICLES (512)
#define DYNAMIC_TIME_STEP (2E-8)
#define GRAVITY (9.81)
#define DENSITY (7952.0)
#define SCAN_LENGTH_RATIO (1.0)

// build platform dimension
#define TABLE_LENGTH (0.001)
#define TABLE_WIDTH (0.0005)

// change the table dimension if more particles is used
// #define TABLE_LENGTH (0.002)
// #define TABLE_WIDTH (0.001)

#define MAX_HEIGHT (0.0004)
#define MIN_HEIGHT (0.0002)
#define ELEMENT_DIAMETER_MEAN (0.000027)
#define ELEMENT_DIAMETER_STD (0.00001/2.355)
#define ELEMENT_DIAMETER_MAX (0.000050)
#define ELEMENT_DIAMETER_MIN (0.000010)
#define LAYER_HEIGHT (0.000120)

// thermal
#define CONDUCTIVITY (19.0)
#define HEAT_TRANSFER_COEFFICIENT (40.0)
#define EMMISIVITY (0.33)
#define STEFAN_BOLTZMANN (5.67E-8)
#define ABSORPTIVITY (0.33)
#define POROSITY (0.55)
#define LASER_POWER (100.0)
#define LASER_SCAN_SPEED (2.0)
#define LASER_SPOT_SIZE (0.000054)
#define PREHEAT_T (363.0)
#define TABLE_T (363.0)
#define ENVIORNMENT_T (363.0)
#define MELTING_T (1700.0)
#define BOILING_T (3130.0)
#define MELTING_LATENT_HEAT (2.99E5)
#define VAPORIZATION_LATENT_HEAT (6.09E6)
#define C_SOLID (452.0)
#define C_LIQUID (815.0)
#define C_GAS (815.0)
#define DELTA_T (100.0)
#define THERMAL_TIME_STEP (5E-9)
#define FLOOR_CONTACT_RATIO (0.1)

#endif
