#ifndef EXTERNALCONATCTFORCE_H
#define EXTERNALCONATCTFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper,
        double m_stiffness, double m_dBar);
	~externalContactForce();

	void computeFc();
	void computeJc();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    VectorXi arrayNum;

    int ind, ind1, ind2;

    double d;
    double dBar;
    double stiffness;

    Matrix3d Id3;

    Vector3d dDdEdge;
    Matrix3d d2DdEdge2;

    double dEnergydD;
    double d2EnergydD2;

    Matrix3d M0;

    Vector3d f;
    VectorXd fVec;

    int nv_0;
    int nv_1;
    int nv_2;

    int edgeLocal;

    Vector3d x_0;
    Vector3d x_1;
    Vector3d x_2;

    Vector3d v_0;
    Vector3d v_1;
    Vector3d v_2;
    Vector3d v_r;
    Vector3d v_rt;

    double mu;
    double epsilonV;

    Vector3d friction;

    Vector3d tK;

    double fVelocity;

    double minDistance;
    int contactIndex;

    Vector3d e_1, e_2, e_3, e_4;

    double normE1;
    double normE2;
    double normE3;
    double normE4;

    Vector3d t_1;
    Vector3d t_2;
    Vector3d t_3;
    Vector3d t_4;

    double d1, d2, d3;

    Vector3d de1;
    Vector3d de2;
    Vector3d de3;

    Matrix3d crossMat(Vector3d a);

    Matrix3d gradT3;
    Matrix3d gradT4;

    Matrix3d de1de1;
    Matrix3d de2de2;
    Matrix3d de3de3;
    Matrix3d de1de2;
    Matrix3d de1de3;
    Matrix3d de2de3;
    Matrix3d de2de1;
    Matrix3d de3de1;
    Matrix3d de3de2;

    VectorXd dx;
    MatrixXd dxdx;

    MatrixXd jacobian;

    Vector3d normDirection;

    VectorXi contactRingIndex;
    VectorXi boundaryRodIndex;

    int startPoint;
    int endPoint;
};

#endif
