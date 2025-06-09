#ifndef EXTERNALCONATCTFORCEPLATE_H
#define EXTERNALCONATCTFORCEPLATE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepperPlate.h"

class externalContactForcePlate
{
public:
	externalContactForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper, 
        double m_stiffness, double m_dBar, double m_mu, double m_epsilonV);
	~externalContactForcePlate();

	void computeFc();
	void computeJc();

private:
	elasticPlate *plate;
    timeStepperPlate *stepper;

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
