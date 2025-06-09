#ifndef EXTERNALCONATCTFORCE_H
#define EXTERNALCONATCTFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticRod &m_rod, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar, double m_mu, double m_epsilonV);
	~externalContactForce();

	void computeFc();
	void computeJc();

private:
	elasticRod *rod;
    timeStepper *stepper;

    int ind;

    double dEnergydD;
    double d2EnergydD2;

    Vector3d f;

    double stiffness;
    double dBar;

    double mu;
    double epsilonV;

    double dt;

    double fVelocity;
    Vector3d tK;

    Vector3d friction;
    Matrix3d frictionJacobian;
    Vector3d dfVelocity;
    Matrix3d dtK;

    Matrix3d Id3;
    Matrix3d IdG;
};

#endif
