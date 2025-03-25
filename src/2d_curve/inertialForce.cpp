#include "inertialForce.h"

inertialForce::inertialForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;

	TotalForceVec = VectorXd::Zero(plate->ndof);
}

inertialForce::~inertialForce()
{
	;
}

void inertialForce::computeFi()
{
	TotalForceVec = VectorXd::Zero(plate->ndof);

	for (int i = 0; i < plate->ndof; i++)
	{
		f = plate->massArray[i] * (plate->x[i] - plate->x0[i]) / ((plate->dt) *(plate->dt))
				- (plate->massArray[i] * plate->u[i])/(plate->dt);

		TotalForceVec(i) = TotalForceVec(i) + f;
		
		stepper->addForce(i, f);
	}
}

void inertialForce::computeJi()
{
	for (int i=0; i<plate->ndof; i++)
    {
		jac = plate->massArray(i)/ ((plate->dt) *(plate->dt));
		stepper->addJacobian(i, i, jac);
	}
}

void inertialForce::setFirstJacobian()
{
	for (int i=0; i<plate->ndof; i++)
	{
		stepper->addJacobian(i, i, 1);
	}
}
