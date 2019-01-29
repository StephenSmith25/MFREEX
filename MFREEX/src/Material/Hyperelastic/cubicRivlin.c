
#include "Material/Hyperelastic/cubicRivlin.h"
#include "m_inverse_small.h"
int cubicRivlin(VEC * stressVoigt, MAT * defGrad, VEC * matParams){


	/*  Initialsie matricies */
	MAT * C = m_get(defGrad->m,defGrad->n);
	MAT * invC = m_get(defGrad->m,defGrad->n);
	MAT * ident = m_get(defGrad->m,defGrad->n);
	m_ident(ident);
	MAT * stress1PKF = m_get(defGrad->m,defGrad->n);
	MAT * stress2PKF = m_get(defGrad->m,defGrad->n);


	MAT * term1 = m_get(defGrad->m,defGrad->n);
	MAT * term2 = m_get(defGrad->m,defGrad->n);
	int dim_stress = stressVoigt->max_dim;


	mtrm_mlt(defGrad,defGrad,C);
	/*  Find inv C */
	m_inverse_small(C,invC);

	/*  Find I1 */
	double I1 = trace(C);
	/*  Find I3 */
	double I3 = determinant(C);
	
	/*  Find I1bar and J */
	double J = sqrt(I3);
	double I1bar = I1*pow(I3,-1.00/3.00);

	/*  Find K1 = dWd/dI1 */
	double K1 = matParams->ve[0] + 2*matParams->ve[1]*(I1bar - 3) + 3*matParams->ve[2]*pow(I1bar-3,2);
	/*  Find P = dWv/dJ */
	double P = matParams->ve[3]*(J-1);


	/*  Distortional part */
	double c1 = 2*K1*pow(I3,-1.00/3.00);
	sm_mlt(I1/3.00,invC,term1);
	m_sub(ident,term1,term2);
	sm_mlt(c1,term2,stress2PKF);

	/*  Volumetric part */
	sm_mlt(P*J,invC,term1);
	m_add(stress2PKF,term1,stress2PKF);



	m_mlt(defGrad,stress2PKF,stress1PKF);

	if ( dim_stress == 5){

		stressVoigt->ve[0] = stress1PKF->me[0][0];
		stressVoigt->ve[1] = stress1PKF->me[1][1];
		stressVoigt->ve[2] = stress1PKF->me[0][1];
		stressVoigt->ve[3] = stress1PKF->me[1][0];
		stressVoigt->ve[4] = stress1PKF->me[2][2];


	}else{

		stressVoigt->ve[0] = stress1PKF->me[0][0];
		stressVoigt->ve[1] = stress1PKF->me[1][1];
		stressVoigt->ve[2] = stress1PKF->me[0][1];
		stressVoigt->ve[3] = stress1PKF->me[1][0];




	}



	/*  free memory */
	M_FREE(C);
	M_FREE(invC);
	M_FREE(ident);
	M_FREE(stress1PKF);
	M_FREE(stress2PKF);
	M_FREE(term1);
	M_FREE(term2);



	return 0;
}
