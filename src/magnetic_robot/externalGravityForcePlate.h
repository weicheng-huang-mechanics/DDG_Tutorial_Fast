#ifndef EXTERNALGRAVITYFORCEPLATE_H
#define EXTERNALGRAVITYFORCEPLATE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepperPlate.h"

class externalGravityForcePlate
{
public:
	externalGravityForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper, Vector3d m_gVector);
	~externalGravityForcePlate();
	
	void computeFg();
	void computeJg();

	void setGravity();

	Vector3d gVector;
	
private:
	elasticPlate *plate;
	timeStepperPlate *stepper;
	
    VectorXd massGravity;
};

#endif
