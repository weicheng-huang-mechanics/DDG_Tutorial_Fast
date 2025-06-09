#ifndef EXTERNALMAGNETICFORCE_H
#define EXTERNALMAGNETICFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class externalMagneticForce
{
public:
	externalMagneticForce(elasticRod &m_rod, timeStepper &m_stepper, 
    Vector3d m_bAVector, Vector3d m_bRVector, double m_muZero, double m_omega);
	~externalMagneticForce();

	void computeFm(double m_currentTime, double m_rol, double m_row, int m_cut);
	void computeJm();
	
private:
	elasticRod *rod;
	timeStepper *stepper;

    Vector3d baVector_ref;
    Vector3d brVector_ref;

    Vector3d baVector;
    Vector3d brVector;
    double muZero;
    double omega;

    Vector3d m1_current;
    Vector3d m2_current;
    Vector3d m3_current;

    Vector3d m1_start;
    Vector3d m2_start;
    Vector3d m3_start;

    Matrix3d Id3;

    Matrix3d dm1de;
    Matrix3d dm2de;
    Matrix3d dm3de;

    Vector3d dm1dtheta;
    Vector3d dm2dtheta;
    Vector3d dm3dtheta;

    Vector3d x1;
    Vector3d x2;

    double edge;

    Matrix3d dMde;
    Vector3d dEde;

    Vector3d dMdtheta;
    double dEdtheta;

    Vector3d tempM3;
    Vector3d tempM2;
    Vector3d tempM1;

    double tempMTheta;

    Vector3d d2m1dtheta2;
    Vector3d d2m2dtheta2;
    Vector3d d2m3dtheta2;

    Matrix3d d2Temp3de2;
    Matrix3d d2Temp1de2;
    Matrix3d d2Temp2de2;

    Vector3d d2Mdtheta2;

    Matrix3d d2Ede2;
    Vector3d d2Ededtheta;
    double d2Edtheta2;

    VectorXd force;
    MatrixXd jacob;
    
    Matrix3d gradientBa;

    Vector3d Mag;
};

#endif
