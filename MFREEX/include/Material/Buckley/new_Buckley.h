#ifndef NEW_BUCKLEY_MATERIAL_H_
#define NEW_BUCKLEY_MATERIAL_H_

typedef struct Buckley_Material{


	// bond
	double Vs;
	double Vp;
	double mu_star;
	double T_inf;
	double T_star;
	double Tf_star;
	double Cv;
	double H0;
	double R;
	double Kb;
	double Gb;
	double Tf;

	// conformational
	double alpha_c;
	double eta_c;
	double Ns_c;
	double boltzmann;

	// slippage
	double lambdaCrit;
	double Ts;
	double gamma0_ref;
	double Cs;
	double T_inf_s;
	double C1_s;
	double C2_s;
	double beta_;




}Buckley_Material;

int new_Buckley_material();


#endif
