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
    Vector3d x_4;

    Vector3d e_1;
    Vector3d e_2;
    Vector3d e_3;
    Vector3d e_4;

    Vector3d n_1;
    Vector3d n_2;

    Vector3d norm_1;
    Vector3d norm_2;

    Vector3d nBar;

    Matrix3d gradN1;
    Matrix3d gradN2;

    Vector3d dEde1;
    Vector3d dEde2;
    Vector3d dEde3;
    Vector3d dEde4;

    Matrix3d d2Ede12;
    Matrix3d d2Ede22;
    Matrix3d d2Ede32;
    Matrix3d d2Ede42;

    Matrix3d d2Ede1de2;
    Matrix3d d2Ede2de1;
    Matrix3d d2Ede1de3;
    Matrix3d d2Ede3de1;
    Matrix3d d2Ede1de4;
    Matrix3d d2Ede4de1;
    Matrix3d d2Ede2de3;
    Matrix3d d2Ede3de2;
    Matrix3d d2Ede2de4;
    Matrix3d d2Ede4de2;
    Matrix3d d2Ede3de4;
    Matrix3d d2Ede4de3;

    VectorXi arrayNum;

    VectorXd forceVec;
    MatrixXd Jbb;

    int ind, ind1, ind2;

    Matrix3d skemMatrix(Vector3d a);
    Matrix3d hessionMatrix_1(Vector3d a);
    Matrix3d hessionMatrix_2(Vector3d a);
    Matrix3d hessionMatrix_3(Vector3d a);
};

#endif
