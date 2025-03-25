#include "elasticStretchingForce.h"
#include <iostream>

elasticStretchingForce::elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;
	
	Jss.setZero(6, 6);
	flocal = VectorXd::Zero(6);
	EA = plate->EA;

	localDOF = VectorXi::Zero(6);

	Id3<<1,0,0,
	     0,1,0,
	     0,0,1;
}

elasticStretchingForce::~elasticStretchingForce()
{
	;
}

void elasticStretchingForce::computeFs()
{
	for (int k = 0; k < plate->edgeNum; k++)
	{
		flocal = VectorXd::Zero(6);

		int nv_1 = plate->v_edgeElement[k].nv_1;
		int nv_2 = plate->v_edgeElement[k].nv_2;

		Vector3d p = plate->getVertex(nv_1);
		Vector3d p1 = plate->getVertex(nv_2);

		Vector3d tangent = (p1 - p) / (p1 - p).norm();

		double edgeLen = (p1 - p).norm();

		double epsX = edgeLen / plate->v_edgeElement[k].refLength - 1.0;

		Vector3d f = EA * tangent * epsX;

		localDOF = plate->v_edgeElement[k].arrayNum;

		for(int i = 0; i < 3; i++)
		{
			int ind1 = localDOF(i);
			stepper->addForce(ind1, - f[i]);
		}

		for(int i = 0; i < 3; i++)
		{
			int ind2 = localDOF(i+3);
			stepper->addForce(ind2, f[i]); // subtracting elastic force
		}

	}
	
}

void elasticStretchingForce::computeJs()
{
	
	for (int k = 0; k < plate->edgeNum; k++)
	{
		Jss.setZero(6,6);

		int nv_1 = plate->v_edgeElement[k].nv_1;
		int nv_2 = plate->v_edgeElement[k].nv_2;

		Vector3d p = plate->getVertex(nv_1);
		Vector3d p1 = plate->getVertex(nv_2);

		Vector3d dxx = p1 - p;

		Vector3d u;
		u = dxx;

		Matrix<double,1,3> v;
		v = u.transpose();

		double refLength = plate->v_edgeElement[k].refLength;
		double len = (p1 - p).norm();

		Matrix3d M0 = EA * ((1/refLength - 1/len) * Id3 + (1/len) * (u*v) / (u.norm() * u.norm()));
		
		Jss.block(0,0,3,3) =  - M0;
		Jss.block(3,3,3,3) =  - M0;
		Jss.block(3,0,3,3) =    M0;
		Jss.block(0,3,3,3) =    M0;

		localDOF = plate->v_edgeElement[k].arrayNum;

		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				stepper->addJacobian(localDOF(i), localDOF(j), - Jss(i,j));
			}
		}
	}
	
}

void elasticStretchingForce::setFirstJacobian()
{
	for (int k = 0; k < plate->edgeNum; k++)
	{
		localDOF = plate->v_edgeElement[k].arrayNum;

		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				stepper->addJacobian(localDOF(i), localDOF(j), 1);
			}
		}
	}
}