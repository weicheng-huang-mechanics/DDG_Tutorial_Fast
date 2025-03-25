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
	
	double EA;
    int ndof;
    
    double xk, yk, xkp1, ykp1, l_k;
    VectorXd flocal;
    VectorXd f;
    MatrixXd Jss;
    int ind1, ind2;
    Vector2d p, p1;

    VectorXi localDOF;

    void localForce();
    void localJacobian();
};

#endif
