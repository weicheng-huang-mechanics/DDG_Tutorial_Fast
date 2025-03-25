#ifndef EXTERNALGRAVITYFORCE_H
#define EXTERNALGRAVITYFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalGravityForce
{
public:
	externalGravityForce(elasticPlate &m_plate, timeStepper &m_stepper, Vector3d m_gVector);
	~externalGravityForce();
	void computeFg();
	void computeJg();

	void setGravity();

	Vector3d gVector;
	
private:
	elasticPlate *plate;
	timeStepper *stepper;
	
    VectorXd massGravity;
};

#endif
