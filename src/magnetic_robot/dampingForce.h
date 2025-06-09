#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class dampingForce
{
public:
	dampingForce(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity);
	~dampingForce();
	void computeFd();
	void computeJd();

private:
	elasticRod *rod;
	timeStepper *stepper;
	double viscosity;
	
    Vector3d u, f;
    int ind;
    double jac;
};

#endif
