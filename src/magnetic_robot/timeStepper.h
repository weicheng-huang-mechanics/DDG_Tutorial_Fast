#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "elasticRod.h"

class timeStepper
{
public:
	timeStepper(elasticRod &m_rod);
	~timeStepper();

	void setZero();
	void addForce(int ind, double p);
	void addJacobian(int ind1, int ind2, double p);
	
	VectorXd ForceVector;
	MatrixXd JacobianMatrix;

private:
	elasticRod *rod;
};

#endif
