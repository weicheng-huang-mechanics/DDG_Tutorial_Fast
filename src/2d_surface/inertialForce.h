#ifndef INERTIALFORCE_H
#define INERTIALFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class inertialForce
{
public:
	inertialForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~inertialForce();
	void computeFi();
	void computeJi();

	void setFirstJacobian();

	VectorXd TotalForceVec;

private:
	elasticPlate *plate;
	timeStepper *stepper;
    			
    int ind1, ind2, mappedInd1, mappedInd2;	
    double f, jac;
};

#endif
