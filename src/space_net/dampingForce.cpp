#include "dampingForce.h"

dampingForce::dampingForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_viscosity)
{
	plate = &m_plate;
    stepper = &m_stepper;

	viscosity = m_viscosity;
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	for (int i = 0; i < plate->nv; i++)
	{
		u1 = plate->getVelocity(i);

		f1 = - u1 * viscosity * plate->massArray[i];

		for (int k = 0; k < 3; k++)
		{
			ind = 3 * i + k;
			
			stepper->addForce(ind, - f1[k]);
		}
	}
}

void dampingForce::computeJd()
{
	for (int i = 0; i < plate->nv; i++)
	{
		jac(0) = - viscosity * plate->massArray[i] / plate->dt;
		jac(1) = - viscosity * plate->massArray[i] / plate->dt;
		jac(2) = - viscosity * plate->massArray[i] / plate->dt;
		
		for (int j = 0; j < 3; j++)
		{
			ind = 3 * i + j;
			stepper->addJacobian(ind, ind, - jac(j));
		}
	}
}

void dampingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->nv; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ind = 3 * i + j;
			stepper->addJacobian(ind, ind, 1);
		}
	}
}
