#ifndef ELASTICBENDINGFORCECOUPLE_H
#define ELASTICBENDINGFORCECOUPLE_H

#include "eigenIncludes.h"
#include "GlobalStepper.h"

class elasticBendingForceCouple
{
public:
	elasticBendingForceCouple(GlobalStepper &m_stepper, double m_EI, 
    int m_rodDOF, int m_plateDOF);
	~elasticBendingForceCouple();

	void computeFb(int m_plate_nv_1, int m_plate_nv_2, int m_plate_nv_3, 
        Vector3d m_plate_v1, Vector3d m_plate_v2, Vector3d m_plate_v3,
        int m_rodIndex, Vector3d m_rod_v1, Vector3d m_rod_v2, Vector3d m_nBar);

	void computeJb(int m_plate_nv_1, int m_plate_nv_2, int m_plate_nv_3, 
        Vector3d m_plate_v1, Vector3d m_plate_v2, Vector3d m_plate_v3,
        int m_rodIndex, Vector3d m_rod_v1, Vector3d m_rod_v2, Vector3d m_nBar);

    void setFirstJacobian(int m_plate_nv_1, int m_plate_nv_2, int m_plate_nv_3, int m_rodIndex);

private:

    GlobalStepper *Globalstepper;

    Matrix3d Id3;
    double EI;

    Vector3d x_1;
    Vector3d x_2;
    Vector3d x_3;
    Vector3d x_4;
    Vector3d x_5;

    Vector3d e_1;
    Vector3d e_2;
    Vector3d e_3;

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

    Matrix3d d2Ede12;
    Matrix3d d2Ede22;
    Matrix3d d2Ede32;

    Matrix3d d2Ede1de2;
    Matrix3d d2Ede2de1;
    Matrix3d d2Ede1de3;
    Matrix3d d2Ede3de1;

    Matrix3d d2Ede2de3;
    Matrix3d d2Ede3de2;

    VectorXi indexVec;

    VectorXd forceVec;
    MatrixXd Jbb;

    int ind, ind1, ind2;

    int rodDOF;
    int plateDOF;

    Matrix3d skemMatrix(Vector3d a);

    Matrix3d hession1;
    Matrix3d hession2;
    Matrix3d hession3;
    Matrix3d hession;
};

#endif
