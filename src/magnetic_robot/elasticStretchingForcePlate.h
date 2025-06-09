#ifndef ELASTICSTRETCHINGFORCEPLATE_H
#define ELASTICSTRETCHINGFORCEPLATE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepperPlate.h"

class elasticStretchingForcePlate
{
public:
	elasticStretchingForcePlate(elasticPlate &m_plate, timeStepperPlate &m_stepper);
	~elasticStretchingForcePlate();
    
	void computeFs();
	void computeJs();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepperPlate *stepper;
	
	double len, refLength;
    double epsX;
    Vector3d u;
    Vector3d dxx;
    Vector3d f;
    Matrix3d Id3;
    Matrix3d M0;
    Matrix<double,1,3> v;
    Matrix<double,6,6> Jss;

    VectorXi arrayNum;
    Vector3d tangentVec;

    Vector3d x_1;
    Vector3d x_2;
    
    double EA;
    int ind, ind1, ind2;	
};

#endif
