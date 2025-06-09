#ifndef GLOBALSTEPPER_H
#define GLOBALSTEPPER_H

#include "elasticPlate.h"
#include "elasticRod.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <GL/glut.h>

// Define the format to printf MKL_INT values
#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

class GlobalStepper
{
public:

	GlobalStepper(int m_uncons_dof_mapped, VectorXi m_isConstrained_global, VectorXi m_fullToUnconsMap_global);
	~GlobalStepper();

	VectorXd GlobalForceVec;
	MatrixXd GlobalJacobianMax;
	VectorXd GlobalMotionVec;

	void setZero();
	void addForce(int ind, double p);
	void addJacobian(int ind1, int ind2, double p);
	void integrator();

	void first_time_PARDISO_setup();

private:

	int total_dof;
    int uncon_dof;
    int con_dof;

    VectorXi isConstrained_global;
    VectorXi fullToUnconsMap_global;

    int mappedInd, mappedInd1, mappedInd2;

    int n;
    int *ia;
    int *ja;
    double *a;

    MKL_INT mtype;

    double *b;
    double *x;

    MKL_INT nrhs;

    void *pt[64]; 

    MKL_INT iparm[64];
    double dparm[64];

    MKL_INT maxfct, mnum, phase, error, msglvl;

    double ddum;

    MKL_INT idum;

    int total_num;

    void pardisoSolver();

    VectorXi num1;
    VectorXi num2;
};

#endif
