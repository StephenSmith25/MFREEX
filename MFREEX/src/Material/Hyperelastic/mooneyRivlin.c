
#include "Material/Hyperelastic/mooneyRivlin.h"
#include <math.h>
#include "m_inverse_small.h"

int mooneyRivlin(VEC * stressVoigt, MAT * defGrad, VEC * params){


	/*  Deformation tensors */
	MAT * C = m_get(defGrad->m, defGrad->n);
	MAT * invC = m_get(defGrad->m, defGrad->n);
	/*  Stress tensors */
	MAT * stress1PKF = m_get(defGrad->m, defGrad->n);
	MAT * stress2PKF = m_get(defGrad->m, defGrad->n);

	MAT * ident = m_get(defGrad->m,defGrad->n);
	m_ident(ident);
	/*  Terms */
	MAT * term1 = m_get(defGrad->m,defGrad->n);
	MAT * term2 = m_get(defGrad->m,defGrad->n);
	MAT * term3 = m_get(defGrad->m,defGrad->n);

	int dim = C->m;


	/*  Material constants */
	double c1 = params->ve[0];
	double c2 = params->ve[1];
	double lambda = params->ve[2];
	double p0 = c1 + 2*c2;

	/*  Find C */
	mtrm_mlt(defGrad,defGrad,C);
	/*  Find C^{-1} */
	m_inverse_small(C,invC);

	/*  I1 and I3 */
	double I1 = trace(C);

	if ( dim == 2)
	{
		I1 = I1 + 1;
	}
	double I3 = determinant(C);



	/*  Stress computation */
	/*  2*c1*d_ij */
	sm_mlt(2*c1,ident,term1);
	/*  I1*ident */
	sm_mlt(I1,ident,term2);
	m_sub(term2,C,term2);
	/*  2*c2*(I1 - C) */
	sm_mlt(2*c2,term2,term2);
	/*  Pressure term */
	sm_mlt(-2*(p0-lambda*log(I3)),invC,term3);

	/*  Find stress */
	m_add(term1,term2,stress2PKF);
	m_add(stress2PKF,term3,stress2PKF);

	m_mlt(defGrad,stress2PKF,stress1PKF);


	// if ( AXI == 1){
	// stressVoigt->ve[4] = stress1PKF->me[2][2];
	// }
	stressVoigt->ve[0] = stress1PKF->me[0][0];
	stressVoigt->ve[1] = stress1PKF->me[1][1];
	stressVoigt->ve[2] = stress1PKF->me[0][1];
	stressVoigt->ve[3] = stress1PKF->me[1][0];


	M_FREE(C);
	M_FREE(invC);
	M_FREE(term1);
	M_FREE(term2);
	M_FREE(term3);
	M_FREE(ident);
	M_FREE(stress2PKF);
	M_FREE(stress1PKF);





	return 0;
}
