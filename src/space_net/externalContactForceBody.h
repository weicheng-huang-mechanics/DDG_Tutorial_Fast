#ifndef EXTERNALCONATCTFORCEBODY_H
#define EXTERNALCONATCTFORCEBODY_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalContactForceBody
{
public:
	externalContactForceBody(elasticPlate &m_plate, timeStepper &m_stepper, double m_boxSize,
		double m_stiffness, double m_dBar);
	~externalContactForceBody();

	void computeFc();
	void computeJc();

    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double boxSize;

    int ind;

    Vector3d dDdEdge;
    Matrix3d Id3;

    double dEnergydD;
    double d2EnergydD2;

    Matrix3d d2DdEdge2;
    Matrix3d d2EdEdge2;

    Vector3d f;

    double stiffness;
    double dBar;

    VectorXi indexArray;
    MatrixXd jacobian;
};

#endif
