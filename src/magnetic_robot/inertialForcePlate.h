#ifndef INERTIALFORCEPLATE_H
#define INERTIALFORCEPLATE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepperPlate.h"

class inertialForcePlate
{
public:
	inertialForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper);
	~inertialForcePlate();

	void computeFi();
	void computeJi();

	void setFirstJacobian();

	VectorXd TotalForceVec;

private:
	elasticPlate *plate;
	timeStepperPlate *stepper;
    			
    int ind1, ind2, mappedInd1, mappedInd2;	
    double f, jac;
};

#endif
