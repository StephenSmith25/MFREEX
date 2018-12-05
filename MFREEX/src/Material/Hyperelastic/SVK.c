/*
 * =====================================================================================
 *
 *       Filename:  SVK.c
 *
 *    Description:  Saint Venant-Kirchhoff Material model 
 *
 *        Version:  1.0
 *        Created:  19/02/18 11:45:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  STEPHEN SMITH
 *   Organization:  QUEEN'S UNIVERSITY BELFAST
 *
 * =====================================================================================
 */
#include "Material/Hyperelastic/SVK.h"

int SVK(VEC * stressVoigt, MAT * defGrad, VEC * params){





	// MAT * strainGL = m_get(defGrad->m,defGrad->m);
	// MAT * ident = m_get(defGrad->m,defGrad->m);
	// MAT * term1 = m_get(defGrad->m,defGrad->m);
	// MAT * term2 = m_get(defGrad->m,defGrad->m);
	// MAT * stress1PKF = m_get(defGrad->m,defGrad->m);
	// MAT * stress2PKF = m_get(defGrad->m,defGrad->m);


	// m_ident(ident);

	int dim = defGrad->m;


	if (dim == 2)
	{
		double F[4];
		double E[3];
	}else{
		double F[9];
		double E[6];
	}


	// Green lagrange strain tensor = 1/2 ( F'F - I);
	double E11 = 0, E12 = 0, E13 = 0, E21 = 0, E22 = 0, E23 = 0, E31 = 0, E32 = 0, E33 = 0;


	// deformation gradient components
	double F11 = 0, F12 = 0, F13 = 0, F21 = 0, F22 = 0, F23 = 0, F31 = 0, F32 = 0, F33 = 0;

	// Components of green lagrange tensor
	if ( dim == 2)
	{
		// row 1
		F11 = defGrad->me[0][0];F12 = defGrad->me[0][1];
		// row 2
		F21 = defGrad->me[1][0]; F22 = defGrad->me[1][1];


		/*Green lagrange strain tensor */
		// row 1
		E11 = 0.5*(F11*F11 + F21*F21 -1 );  E12 = 0.5*(F11*F12 + F21*F22);
		// row 2
		E21 = E12 ; E22 = 0.5*(F12*F12 + F22*F22 - 1);

	}else{

		// row 1
		F11 = defGrad->me[0][0];
		F12 = defGrad->me[0][1];
		F13 = defGrad->me[0][2];
		// row2
		F21 = defGrad->me[1][0];
		F22 = defGrad->me[1][1];
		F23 = defGrad->me[1][2];
		// row 3
		F31 = defGrad->me[2][0];
		F32 = defGrad->me[2][1];
		F33 = defGrad->me[2][2];


		/*Green lagrange strain tensor */
		// row 1
		E11 = 0.5*(F11*F11 + F12*F12 + F13*F13 -1 ); E12 = 0.5*(F11*F21 + F12*F22 + F13*F23 ); E13 = 0.5*(F11*F31 + F12*F32 + F13*F33);
		// row 2
		E21 = E12 ; E22 = 0.5*(F21*F21 + F22*F22 + F23*F23 -1 ); E23 = 0.5*(F21*F31+F22*F32+F23*F33);
		// row 3
		E31 = E13 ; E32 = E23 ; E33 = 0.5*(F31*F31 + F32*F32 + F33*F33-1);
	}




	// mtrm_mlt(defGrad,defGrad,strainGL);
	
	// m_sub(strainGL,ident,strainGL);
	// sm_mlt(0.5000,strainGL,strainGL);

	double traceE = E11 + E22 + E33;



	// second piola kirchoff stress
	double S11 = 0, S12 = 0, S13 = 0, S21 = 0, S22 = 0, S23 = 0, S31 = 0, S32 = 0, S33 = 0;


	// First piola kirchoff stress
	double P11 = 0, P12 = 0, P13 = 0, P21 = 0, P22 = 0, P23 = 0, P31 = 0, P32 = 0, P33 = 0;



	// Find second piola kirchoff stress based on st venant kirchoff material model

	if ( dim == 2)
	{
		// row 1
		S11 = traceE*params->ve[0] + 2*params->ve[1]*E11;
		S12 = 2*params->ve[1]*E12;
		// row 2
		S21 = S12;
		S22 = traceE*params->ve[0] + 2*params->ve[1]*E22;

	}else{
		// row 1
		S11 = traceE*params->ve[0] + 2*params->ve[1]*E11;
		S12 = 2*params->ve[1]*E12;
		S13 = 2*params->ve[1]*E13;
		// row 2
		S21 = 2*params->ve[1]*E21;
		S22 = traceE*params->ve[0] + 2*params->ve[1]*E22;
		S23 = 2*params->ve[1]*E23;
		// row 3
		S31 = 2*params->ve[1]*E31;
		S32 = 2*params->ve[1]*E32;
		S33 = traceE*params->ve[0] + 2*params->ve[1]*E33;

	}



	int dim_s = stressVoigt->max_dim;


	char * returnType = "1PKF";

	if ( strcmp(returnType, "1PKF") == 0 )
	{
		// first first piola stress

		if ( dim == 2)
		{
			// row 1
			P11 = F11*S11 + F12*S21;
			P12 = F11*S12 + F12*S22;
			// row 2
			P21 = F21*S11 + F22*S21;
			P22 = F21*S12 + F22*S22;

		}else{
			// row 1
			P11 = F11*S11 + F12*S21 + F13*S31;
			P12 = F11*S12 + F12*S22 + F13*S32;
			P13 = F11*S13 + F12*S23 + F13*S33;
			// row 2
			P21 = F21*S11 + F22*S21 + F23*S31;
			P22 = F21*S12 + F22*S22 + F23*S32;
			P23 = F21*S13 + F22*S23 + F23*S33;
			// row 3
			P31 = F31*S11 + F32*S21 + F33*S31;
			P32 = F31*S12 + F32*S22 + F33*S32;
			P33 = F31*S13 + F32*S23 + F33*S33;

		}

		if ( dim_s == 4)
		{
			stressVoigt->ve[0] = P11;
			stressVoigt->ve[1] = P22;
			stressVoigt->ve[2] = P12;
			stressVoigt->ve[3] = P21;

		}else if ( dim_s == 5)
		{
			stressVoigt->ve[0] = P11;
			stressVoigt->ve[1] = P22;
			stressVoigt->ve[2] = P12;
			stressVoigt->ve[3] = P21;
			stressVoigt->ve[4] = P33;

		}else if ( dim_s == 9) {


			stressVoigt->ve[0] = P11;
			stressVoigt->ve[1] = P22;
			stressVoigt->ve[2] = P33;
			stressVoigt->ve[3] = P23;
			stressVoigt->ve[4] = P13;
			stressVoigt->ve[5] = P12;
			stressVoigt->ve[6] = P21;
			stressVoigt->ve[7] = P31;
			stressVoigt->ve[8] = P32;
		}else{
			fprintf(stderr,"stress error, dimension of stress voigt not set \n");
		}

	}



	// sm_mlt(params->ve[0]*traceE,ident,term1);
	// sm_mlt(2*params->ve[1],strainGL,term2);

	// m_add(term1,term2,stress2PKF);

	// m_mlt(defGrad,stress2PKF,stress1PKF);




	// if ( AXI == 1){
	// stressVoigt->ve[4] = stress1PKF->me[2][2];
	// }
	






	// m_free(strainGL);
	// m_free(ident);
	// m_free(term1);
	// m_free(term2);
	// m_free(stress1PKF);
	// m_free(stress2PKF);



	return 0;







}

