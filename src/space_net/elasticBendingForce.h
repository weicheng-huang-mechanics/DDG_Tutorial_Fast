#ifndef ELASTICBENDINGFORCE_H
#define ELASTICBENDINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticBendingForce
{
public:
	elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticBendingForce();

	void computeFb();
	void computeJb();
    void setFirstJacobian();
    
    VectorXd TotalForceVec;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    Matrix3d Id3;
    double EI;

    Vector3d x_1;
    Vector3d x_2;
    Vector3d x_3;

    Vector3d e_1;
    Vector3d e_2;

    Vector3d t_1;
    Vector3d t_2;

    double norm_1;
    double norm_2;

    double voroni;

    Vector3d nBar;

    Matrix3d gradN1;
    Matrix3d gradN2;

    Vector3d dEde1;
    Vector3d dEde2;
  
    Matrix3d d2Ede12;
    Matrix3d d2Ede22;

    Matrix3d d2Ede1de2;
    Matrix3d d2Ede2de1;

    Matrix3d hessionMatrix1;
    Matrix3d hessionMatrix2;
    Matrix3d hessionMatrix3;

    VectorXi arrayNum;

    VectorXd forceVec;
    MatrixXd Jbb;

    int ind, ind1, ind2;

    double phi;
    double dKdPhi;
};

#endif
