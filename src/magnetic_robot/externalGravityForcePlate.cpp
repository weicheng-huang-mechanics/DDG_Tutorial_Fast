#include "externalGravityForcePlate.h"

externalGravityForcePlate::externalGravityForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper, Vector3d m_gVector)
{
	plate = &m_plate;
	stepper = &m_stepper;
	gVector = m_gVector;
	
	setGravity();
}

externalGravityForcePlate::~externalGravityForcePlate()
{
	;
}

void externalGravityForcePlate::computeFg()
{
	for (int i=0; i < plate->ndof; i++)
	{
		stepper->addForce(i, -massGravity[i]); // subtracting gravity force
	}	
}

void externalGravityForcePlate::computeJg()
{
	;
}

void externalGravityForcePlate::setGravity()
{
	massGravity = VectorXd::Zero(plate->ndof);

	for (int i = 0; i < plate->nv; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			int ind = 3 * i + k;
			massGravity[ind] = gVector[k] * plate->massArray[ind];
		}
	}
}
