#ifndef DAMPINGFORCEPLATE_H
#define DAMPINGFORCEPLATE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepperPlate.h"

class dampingForcePlate
{
public:
	dampingForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper, double m_viscosity);
	~dampingForcePlate();

	void computeFd();
	void computeJd();

	void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepperPlate *stepper;

    double viscosity;

    Vector3d u, f;

    Vector3d jac;

    double ind;

};

#endif
