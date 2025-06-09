#include "dampingForce.h"
#include <iostream>

dampingForce::dampingForce(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	for (int i=0;i<rod->nv;i++)
	{
		u = rod->getVelocity(i);

		f = - 2 * u * viscosity * rod->massArray(4*i);
		
		for (int k = 0; k < 3; k++)
		{
			ind = 4*i + k;
			stepper->addForce(ind, - f[k]); // subtracting external force
		}
	}
}

void dampingForce::computeJd()
{
	// Remember that dF/dx = 1/dt * dF/dv 
	for (int i=0;i<rod->nv;i++)
	{
		u = rod->getVelocity(i);

		jac = - 2 * viscosity * rod->massArray(4*i) / rod->dt;
		
		for (int kx = 0; kx < 3; kx++)
		{
			ind = 4 * i + kx;

			stepper->addJacobian(ind, ind, - jac); // subtracting external force
		}
	}
}
