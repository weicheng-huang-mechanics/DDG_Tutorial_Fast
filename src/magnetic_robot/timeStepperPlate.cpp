#include "timeStepperPlate.h"

timeStepperPlate::timeStepperPlate(elasticPlate &m_plate)
{
	plate = &m_plate;

	ForceVector.setZero(plate->ndof, 1);
    JacobianMatrix.setZero(plate->ndof, plate->ndof);

    //cout << " plate " << plate->ndof << endl;
}

timeStepperPlate::~timeStepperPlate()
{
	;
}

void timeStepperPlate::addForce(int ind, double p)
{
	ForceVector[ind] = ForceVector[ind] + p; // subtracting elastic force
}

void timeStepperPlate::addJacobian(int ind1, int ind2, double p)
{
	JacobianMatrix(ind1, ind2) = JacobianMatrix(ind1, ind2) + p;
}

void timeStepperPlate::setZero()
{
    ForceVector.setZero(plate->ndof, 1);
    JacobianMatrix.setZero(plate->ndof, plate->ndof);
}