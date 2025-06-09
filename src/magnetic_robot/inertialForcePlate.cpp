#include "inertialForcePlate.h"

inertialForcePlate::inertialForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;
}

inertialForcePlate::~inertialForcePlate()
{
	;
}

void inertialForcePlate::computeFi()
{
	for (int i = 0; i < plate->ndof; i++)
	{
		f = plate->massArray[i] * (plate->x[i] - plate->x0[i]) / ((plate->dt) *(plate->dt))
				- (plate->massArray[i] * plate->u[i])/(plate->dt);
		
		stepper->addForce(i, f);
	}
}

void inertialForcePlate::computeJi()
{
	for (int i=0; i<plate->ndof; i++)
    {
		jac = plate->massArray(i)/ ((plate->dt) *(plate->dt));
		stepper->addJacobian(i, i, jac);
	}
}

void inertialForcePlate::setFirstJacobian()
{
	for (int i=0; i<plate->ndof; i++)
	{
		stepper->addJacobian(i, i, 1);
	}
}
