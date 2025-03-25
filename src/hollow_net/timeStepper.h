#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "elasticPlate.h"

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

class timeStepper
{
public:

	timeStepper(elasticPlate &m_plate);
	~timeStepper();

	VectorXd GlobalForceVec;
	MatrixXd GlobalJacobianMax;
	VectorXd GlobalMotionVec;

    SparseMatrix<double> sm1;

	void setZero();
	void addForce(int ind, double p);
	void addJacobian(int ind1, int ind2, double p);
	void integrator();

	void first_time_PARDISO_setup();

private:
	elasticPlate *plate;

	int total_dof;
	int uncon_dof;
	int con_dof;

	int mappedInd, mappedInd1, mappedInd2;

    int n;
    int *ia;
    int *ja;
    double *a;

    int *ia1;
    int *ja1;

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

    VectorXi num11;
    VectorXi num22;
};

#endif
