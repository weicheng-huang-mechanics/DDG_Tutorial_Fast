#include "timeStepper.h"

timeStepper::timeStepper(elasticPlate &m_plate)
{
	plate = &m_plate;

	total_dof = plate->ndof;
	uncon_dof = plate->uncons;
	con_dof = plate->ncons;

	GlobalForceVec.setZero(uncon_dof, 1);
    GlobalMotionVec.setZero(uncon_dof, 1);
    GlobalJacobianMax.setZero(uncon_dof, uncon_dof);

    n = uncon_dof;
}

timeStepper::~timeStepper()
{
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
}

void timeStepper::addForce(int ind, double p)
{
	if (plate->getIfConstrained(ind) == 0) // free dof
	{
		mappedInd = plate->fullToUnconsMap[ind];
		GlobalForceVec[mappedInd] = GlobalForceVec[mappedInd] + p; // subtracting elastic force
	}
}

void timeStepper::addJacobian(int ind1, int ind2, double p)
{
	mappedInd1 = plate->fullToUnconsMap[ind1];
	mappedInd2 = plate->fullToUnconsMap[ind2];

	if (plate->getIfConstrained(ind1) == 0 && plate->getIfConstrained(ind2) == 0) // both are free
	{
		GlobalJacobianMax(mappedInd1, mappedInd2) = GlobalJacobianMax(mappedInd1, mappedInd2) + p;
	}
}

void timeStepper::setZero()
{
	GlobalForceVec.setZero(uncon_dof, 1);
    GlobalMotionVec.setZero(uncon_dof, 1);

    for (MKL_INT i = 0; i < total_num; i++)
	{
		GlobalJacobianMax(num1(i), num2(i)) = 0.0;
	}
}

void timeStepper::integrator()
{
	pardisoSolver();
}

void timeStepper::first_time_PARDISO_setup()
{
	ia = new int[n+1];
	ia[0] = 1;

	int temp = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			if (GlobalJacobianMax(i,j) != 0 || i == j)
			{
				temp = temp + 1;
			}
		}
		ia[i+1] = temp+1;
	}

	// prepare for ja and a
	total_num = ia[n];
	ja = new int[total_num];
	a = new double[total_num];

	num1 = VectorXi::Zero(total_num);
	num2 = VectorXi::Zero(total_num);

	int temp2 = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			if( GlobalJacobianMax(i,j) != 0 || i == j)
			{
				ja[temp2] = j + 1;
				a[temp2] = GlobalJacobianMax(i,j);

				num1(temp2) = i;
				num2(temp2) = j;

				temp2 = temp2 + 1;
			}
		}
	}

    mtype = -2;       /* Real unsymmetric matrix */

    /* RHS and solution vectors. */
    b = new double[n];
    x = new double[n];

    nrhs = 1;     /* Number of right hand sides. */
   
	/* Setup Pardiso control parameters. */
    for (MKL_INT i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }

    iparm[0] = 0;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */


    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information  */
    error = 0;            /* Initialize error flag */

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
    for (MKL_INT i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: " IFORMAT, error);
        exit (1);
    }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = " IFORMAT, iparm[17]);
    printf ("\nNumber of factorization MFLOPS = " IFORMAT, iparm[18]);

    cout << endl;

    cout << "total num: " << total_num << "; n: " << n << endl;
}

void timeStepper::pardisoSolver()
{
    for (MKL_INT i = 0; i < total_num; i++)
	{
		a[i] = GlobalJacobianMax(num1(i), num2(i));
	}

	for (MKL_INT i = 0; i < n; i++)
	{
		b[i] = GlobalForceVec(i);
	}

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: " IFORMAT, error);
        exit (2);
    }
    //printf ("\nFactorization completed ... ");

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
    phase = 33;

  	//  Loop over 3 solving steps: Ax=b, AHx=b and ATx=b
  	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
  			&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);

    for (int i = 0; i < n; i++)
    {
    	GlobalMotionVec(i) = x[i];
    }
}