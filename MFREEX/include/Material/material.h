#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "matrix.h"
#include "matrix2.h"
//#include "Material/Hyperelastic/hyperelastic_materials.h"


// stores the state of a material in terms of its deformation, 
// and material properties
typedef enum MATERIAL_TYPE{
	HYPERELASTIC,
	PLASTIC,
	BUCKLEY
}MATERIAL_TYPE;

typedef enum MATERIAL_FORMULATION
{
	MAT_ST_VENANT_KIRCHOFF=1,
	MAT_RIVLIN=2,
	MAT_YEOH=6,
	MAT_J2_PLASTICITY=3,
	MAT_BUCKLEY=4

}MATERIAL_FORMULATION;

typedef enum HYPERELASTIC_LAW
{
	MOONEY_RIVLIN,
	NEO_HOOKEAN,
	CUBIC
}HYPERELASTIC_LAW;

typedef enum PLASTIC_LAW
{
	J2,
	DRUCKER_PRAGER

}PLASTIC_LAW;



typedef struct MATERIAL
{
	MATERIAL_TYPE material_type;
	void * MATERIAL_LAW;
	MATERIAL_FORMULATION material;
	
	VEC * params;




}MATERIAL;


typedef struct MOONEY_RIVLIN_MATERIAL
{
	// FUNCTION POINTER
	//int (*GET_STRESS)(VEC*,MAT*,VEC*) = &mooneyRivlin;
	// CONSTANTS
	double c1,c2,kappa;


}MOONEY_RIVLIN_MATERIAL;


typedef struct state_variables
{


	/* ------------------------------------------*/
	/* -----------------Deformation--------------*/
	/* ------------------------------------------*/
	
	// Deformation gradient
	MAT * F;
	MAT * invF;
	MAT * Fdot;

	// Determinant of F
	double Jacobian;

	// Deformation tensors
	MAT * B;
	MAT * C;

	// Polar Decomposition
	MAT * U;
	MAT * V;
	MAT * R;
	MAT * delta_R;
	MAT * Vdot;



	// Rate measures
	MAT * L;
	MAT * D;
	MAT * W;
	MAT * Omega;
	double div_v;


	// Stress
	MAT * sigma;
	MAT * sigma_R;


	// rotated tensors
	MAT * d;
	MAT * T;
	MAT * omega; 



	// temperature
	double temperature;

	/* ------------------------------------------*/
	/* -----------------Plasticity--------------*/
	/* ------------------------------------------*/
	MAT * alpha;
	MAT * back_stress;
	MAT * eta;
	MAT * d_el;
	MAT * d_pl;
	MAT * S_trial;

	double Deps;


	/* ------------------------------------------*/
	/* -----------------Buckley------------------*/
	/* ------------------------------------------*/

	// Eigen values
	VEC * eigValDBar;
	VEC * eigValVBar;
	MAT * eigVecDBar;
	MAT * eigVecVBar;
	VEC * lambdaDot;
	MAT * dbar;
	MAT * EP_bar;
	// Isochoric components
	MAT * Dbar;
	MAT * Lbar;
	MAT * Fbar;
	MAT * Vbar;
	MAT * Ubar;
	MAT * delta_Ubar;
	MAT * delta_ep_bar;

	// Conformational branch network left cauchy green tensor
	MAT * Bbar;
	MAT * Dn;
	MAT * ep_n;
	// Network stretch
	VEC * lambdaNBar;


	// mean stress
	double mSigma;
	// max stretch
	double lambdaNMax;

	// currently conformational gamma
	double gamma; 

	// critical network stretch
	double critLambdaBar;
	// relaxation time
	double tau;

	// Buckley output stresses
	MAT * Sb;
	MAT * Sc;

	// rotated stresses
	MAT * Sb_R;
	MAT * Sc_R;

	/* ------------------------------------------*/
	/* -----------------Workspace----------------*/
	/* ------------------------------------------*/

	// Workspace matricies
	MAT * m_temp1;
	MAT * m_temp2;
	MAT * m_temp3;
	MAT * m_temp4;

	// Workspace vectors
	VEC * v_temp1;
	VEC * v_temp2;
	VEC * v_temp3;
	VEC * v_temp4;



} state_variables;

state_variables ** new_material_states(double * temperatures, int num_nodes, int is_buckley,
 int is_plastic, int dim, int is_AXI);


MATERIAL * create_new_material(MATERIAL_TYPE material_type, void * MATERIAL_LAW, VEC * params);



state_variables * new_material_state( double temperature, MATERIAL_TYPE mat_type, int dim, int is_AXI);



#endif