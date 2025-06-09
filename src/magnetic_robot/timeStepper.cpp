#include "timeStepper.h"

timeStepper::timeStepper(elasticRod &m_rod)
{
	rod = &m_rod;

	ForceVector.setZero(rod->ndof, 1);
    JacobianMatrix.setZero(rod->ndof, rod->ndof);

    //cout << "rod " << rod->ndof << endl;
}

timeStepper::~timeStepper()
{
	;
}

void timeStepper::addForce(int ind, double p)
{
	ForceVector(ind) = ForceVector(ind) + p;
}

void timeStepper::addJacobian(int ind1, int ind2, double p)
{
	JacobianMatrix(ind1, ind2) = JacobianMatrix(ind1, ind2) + p;
}

void timeStepper::setZero()
{
	ForceVector.setZero(rod->ndof, 1);
    JacobianMatrix.setZero(rod->ndof, rod->ndof);
}
