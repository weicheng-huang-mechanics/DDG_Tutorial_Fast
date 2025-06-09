#include "inertialForce.h"

inertialForce::inertialForce(elasticRod &m_rod, timeStepper &m_stepper)
{
	rod = &m_rod;
	stepper = &m_stepper;
}

inertialForce::~inertialForce()
{
	;
}

void inertialForce::computeFi()
{
	// TODO: we should not need to compute this at every iteration.
	// We should compute and store it in iteration 1 and then reuse it.
	for (int i=0; i<rod->ndof; i++)
	{
		f = rod->massArray[i] * (rod->x[i] - rod->x0[i]) / ((rod->dt) *(rod->dt))
				- (rod->massArray[i] * rod->u[i])/(rod->dt);
		stepper->addForce(i, f);
	}
}

void inertialForce::computeJi()
{
	for (int i=0; i<rod->ndof; i++)
    {
		jac = rod->massArray(i)/ ((rod->dt) *(rod->dt));
		stepper->addJacobian(i, i, jac);
	}
}

void inertialForce::setFirstJacobian()
{
	for (int i=0; i<rod->ndof; i++)
    {
		stepper->addJacobian(i, i, 1);
	}
}