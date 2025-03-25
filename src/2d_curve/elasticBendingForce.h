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

private:

	elasticPlate *plate;
    timeStepper *stepper;

    double xkm1, ykm1, xk, yk, xkp1, ykp1, l_k, e0e1, curvature0;
    VectorXd flocal;
    VectorXd f;
    MatrixXd Jbb;
    int ind, ind0, ind1, ind2;
    Vector2d p0, p, p1;
    double EI_local;

    Matrix2d Id2;

    VectorXi localDOF;

    Vector2d t_1, t_2;

    double norm_1, norm_2;

    Vector2d nBar;

    Matrix2d gradN1;
    Matrix2d gradN2;

    Vector2d dEde1;
    Vector2d dEde2;
  
    Matrix2d d2Ede12;
    Matrix2d d2Ede22;

    Matrix2d d2Ede1de2;
    Matrix2d d2Ede2de1;

    Matrix2d hessionMatrix1;
    Matrix2d hessionMatrix2;
    
    VectorXd computeBendingForce(double xa, double ya, double xi, double yi, double xb, double yb, 
        double lBar, double kappaBar);
    MatrixXd computeBendingJacobian(double xa, double ya, double xi, double yi, double xb, double yb, 
        double lBar, double kappaBar);
    VectorXd ListVec(double a1, double a2, double a3, 
        double a4, double a5, double a6);
    MatrixXd ListMat(VectorXd a1, VectorXd a2, VectorXd a3, 
        VectorXd a4, VectorXd a5, VectorXd a6);

};

#endif
