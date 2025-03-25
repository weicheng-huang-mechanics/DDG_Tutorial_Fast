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

    VectorXd totalForce;

private:

	elasticPlate *plate;
    timeStepper *stepper;

    double xkm1, ykm1, xk, yk, xkp1, ykp1;
    double l_k, curvature1, curvature2;
    double localRadius;

    VectorXd flocal;
    MatrixXd Jbb;

    int ind0, ind1, ind2;
    Vector2d p0, p, p1;
    double EI, mu;

    VectorXi localDOF;

    VectorXd computeBendingForce(double xa, double ya, double xi, double yi, double xb, double yb, 
        double mu, double lBar, double kappaBar1, double kappaBar2);
    MatrixXd computeBendingJacobian(double xa, double ya, double xi, double yi, double xb, double yb, 
        double mu, double lBar, double kappaBar1, double kappaBar2);
    VectorXd ListVec(double a1, double a2, double a3, 
        double a4, double a5, double a6);
    MatrixXd ListMat(VectorXd a1, VectorXd a2, VectorXd a3, 
        VectorXd a4, VectorXd a5, VectorXd a6);
};

#endif
