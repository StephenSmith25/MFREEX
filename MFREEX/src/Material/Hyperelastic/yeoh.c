
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  yeoh
 *  Description: Yeoh Material model ( Hyperelastic Law )  
 * =====================================================================================
 */
#include "Material/Hyperelastic/yeoh.h"

#include "determinant.h"
#include "trace.h"
#include <math.h>



int yeoh (VEC * stressVoigt, MAT * defGrad, VEC * params  )
{




	// Initialise matricies

	// Find B

	MAT * B = mmtr_mlt(defGrad,defGrad,MNULL);



	// Find Jacobian ( Determinant of the deformation gradient )


	double jacobian = determinant(defGrad);



	// Find B*, B* = 	
	MAT * Bstar = sm_mlt(1.00/pow(jacobian,(2.00/3.00)),B,MNULL);


	// Find trace of B*
	MAT * IDENT = m_get(3,3);
	IDENT = m_ident(IDENT);


	// Find Volumetric and deviatoric parts of Bstar
	double traceBstar = Bstar->me[0][0] + Bstar->me[1][1] + Bstar->me[2][2] ; 
	MAT * volBstar = sm_mlt((double) (1.0/3.0)*traceBstar,IDENT,MNULL);	
	MAT * devBstar = m_sub(Bstar,volBstar,MNULL ) ; 

	double I1star = Bstar->me[0][0] + Bstar->me[1][1] + Bstar->me[2][2];


	m_ident(IDENT);
	// Find Cauchy Stress
	MAT * stressCauchyVol = sm_mlt(params->ve[3]*(jacobian-1),IDENT,MNULL);
	MAT * stressCauchyDev = sm_mlt((2.00/jacobian)*(params->ve[0] + 2*params->ve[1]*(I1star - 3) + 3*params->ve[2]*pow((I1star-3),2)),devBstar,MNULL);
	MAT * stressCauchy = m_add(stressCauchyVol,stressCauchyDev,MNULL);


	// Stress measure requires pullback, find inv(defGrad)

	MAT * invDefGrad = m_inverse(defGrad,MNULL);
	// find Kirchkoff stress ( volumetrically scaled )
	MAT * stressKirchkoff = sm_mlt(jacobian,stressCauchy,MNULL);
	MAT * stress1PKF = mmtr_mlt(stressKirchkoff,invDefGrad,MNULL);
	MAT * stress2PKF = m_mlt(invDefGrad,stress1PKF,MNULL);


	stressVoigt->ve[0] = stress2PKF->me[0][0];
	stressVoigt->ve[1] = stress2PKF->me[1][1];
	stressVoigt->ve[2] = stress2PKF->me[0][1];
	stressVoigt->ve[3] = stress2PKF->me[2][2];



	M_FREE(B);
	M_FREE(Bstar);
	M_FREE(IDENT);
	M_FREE(volBstar);
	M_FREE(devBstar);
	M_FREE(stressCauchyVol);
	M_FREE(stressCauchyDev);
	M_FREE(stressCauchy);
	M_FREE(invDefGrad);
	M_FREE(stressKirchkoff);
	M_FREE(stress1PKF);
	M_FREE(stress2PKF);


	return 0;


}		/* -----  end of function yeoh  ----- */
