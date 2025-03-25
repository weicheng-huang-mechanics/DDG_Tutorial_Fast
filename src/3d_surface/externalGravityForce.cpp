#include "externalGravityForce.h"

externalGravityForce::externalGravityForce(elasticPlate &m_plate, timeStepper &m_stepper, Vector3d m_gVector)
{
	plate = &m_plate;
	stepper = &m_stepper;
	gVector = m_gVector;
	setGravity();
}

externalGravityForce::~externalGravityForce()
{
	;
}

void externalGravityForce::computeFg()
{

	/*
	massGravity = VectorXd::Zero(plate->ndof);
	gVector[2] = 100 * sin(15 * m_time) - 100.0;

	for (int i = 0; i < plate->nv; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			int ind = 3 * i + k;
			massGravity[ind] = gVector[k] * plate->massArray[ind];
		}
	}
	*/

	for (int i=0; i < plate->ndof; i++)
	{
		stepper->addForce(i, -massGravity[i]); // subtracting gravity force
	}	
}

void externalGravityForce::computeJg()
{
	;
}

void externalGravityForce::setGravity()
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
