#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class dampingForce
{
public:
	dampingForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_viscosity);
	~dampingForce();
	void computeFd();
	void computeJd();

	void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double viscosity;

    Vector2d u1;
    Vector2d u2;
    Vector2d f1;
    Vector2d f2;

    Vector2d jac;

    double ind;

    int index1;
    int index2;

    double edgeLength;

    double sectionArea;

};

#endif
