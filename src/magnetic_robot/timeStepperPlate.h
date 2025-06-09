#ifndef TIMESTEPPERPLATE_H
#define TIMESTEPPERPLATE_H

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

class timeStepperPlate
{
public:

	timeStepperPlate(elasticPlate &m_plate);
	~timeStepperPlate();

	VectorXd ForceVector;
    MatrixXd JacobianMatrix;

	void setZero();
	void addForce(int ind, double p);
	void addJacobian(int ind1, int ind2, double p);

private:
	elasticPlate *plate;
};

#endif
