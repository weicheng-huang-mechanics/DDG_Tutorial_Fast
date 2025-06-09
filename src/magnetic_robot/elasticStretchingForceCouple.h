#ifndef ELASTICSTRETCHINGFORCECOUPLE_H
#define ELASTICSTRETCHINGFORCECOUPLE_H

#include "eigenIncludes.h"
#include "GlobalStepper.h"

class elasticStretchingForceCouple
{
public:
	elasticStretchingForceCouple(GlobalStepper &m_stepper, double m_EA, 
    int m_rodDOF, int m_plateDOF);
	~elasticStretchingForceCouple();

	void computeFs(int index1, int index2, Vector3d x1, Vector3d x2, double refLen);
	void computeJs(int index1, int index2, Vector3d x1, Vector3d x2, double refLen);
    void setFirstJacobian(int index1, int index2);

private:

    GlobalStepper *Globalstepper;

	double len, refLength;
    double epsX;
    Vector3d tangent;
    Vector3d u;
    Vector3d dxx;
    Vector3d f;
    Matrix3d Id3;
    Matrix3d M0;
    Matrix<double,1,3> v;
    Matrix<double,6,6> Jss;
    
    double EA;
    int ind, ind1, ind2;

    int rodDOF;
    int plateDOF;

    VectorXi indexVec;
};

#endif
