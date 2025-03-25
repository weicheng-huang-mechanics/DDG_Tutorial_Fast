#ifndef ELASTICSTRETCHINGFORCE_H
#define ELASTICSTRETCHINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticStretchingForce
{
public:
	elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticStretchingForce();
	void computeFs();
	void computeJs();
    
    void setFirstJacobian();
    
private:
	elasticPlate *plate;
    timeStepper *stepper;

    VectorXi localDOF;

    VectorXd flocal;
    MatrixXd Jss;

    double EA;

    Matrix3d Id3;

    void localForce();
    void localJacobian();
};

#endif
