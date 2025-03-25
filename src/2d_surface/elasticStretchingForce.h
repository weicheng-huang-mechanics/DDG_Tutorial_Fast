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

    VectorXd totalForce;
    
private:
	elasticPlate *plate;
    timeStepper *stepper;
	
	double EA;
    double mu;
    int ndof;
    
    double xk, yk, xkp1, ykp1, l_k;
    double rBar1, rBar2;

    VectorXd flocal;
    MatrixXd Jss;

    int ind1, ind2;

    Vector2d p, p1;

    VectorXi localDOF;

    VectorXd computeStretchingForce(double xa, double ya, double xb, double yb, 
        double mu, double lBar, double rBar1, double rBar2);
    MatrixXd computeStretchingJacobian(double xa, double ya, double xb, double yb, 
        double mu, double lBar, double rBar1, double rBar2);
    VectorXd ListVec(double a1, double a2, double a3, double a4);
    MatrixXd ListMat(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4);
};

#endif
