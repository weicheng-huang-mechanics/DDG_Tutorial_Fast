#include "dampingForcePlate.h"

dampingForcePlate::dampingForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper, double m_viscosity)
{
	plate = &m_plate;
    stepper = &m_stepper;

	viscosity = m_viscosity;
}

dampingForcePlate::~dampingForcePlate()
{
	;
}

void dampingForcePlate::computeFd()
{
	for (int i = 0; i < plate->nv; i++)
	{
		u = plate->getVelocity(i);

		f = - 2 * u * viscosity * plate->massArray(3*i);

		for (int k = 0; k < 3; k++)
		{
			ind = 3 * i + k;
			
			stepper->addForce(ind, - f[k]);
		}
	}
}

void dampingForcePlate::computeJd()
{
	for (int i = 0; i < plate->nv; i++)
	{
		jac(0) = - 2 * viscosity * plate->massArray(3*i) / plate->dt;
		jac(1) = - 2 * viscosity * plate->massArray(3*i) / plate->dt;
		jac(2) = - 2 * viscosity * plate->massArray(3*i) / plate->dt;
		
		for (int j = 0; j < 3; j++)
		{
			ind = 3 * i + j;
			stepper->addJacobian(ind, ind, - jac(j));
		}
	}
}

void dampingForcePlate::setFirstJacobian()
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
