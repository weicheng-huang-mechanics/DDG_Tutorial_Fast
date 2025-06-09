#include "elasticStretchingForceCouple.h"

elasticStretchingForceCouple::elasticStretchingForceCouple(GlobalStepper &m_stepper, double m_EA, 
	int m_rodDOF, int m_plateDOF)
{
	Globalstepper = &m_stepper;

	f.setZero(3);
	Jss.setZero(6, 6);

	indexVec = VectorXi::Zero(6);

	Id3<<1,0,0,
	     0,1,0,
	     0,0,1;

	EA = m_EA;

	rodDOF = m_rodDOF;
	plateDOF = m_plateDOF;
}

elasticStretchingForceCouple::~elasticStretchingForceCouple()
{
	;
}

void elasticStretchingForceCouple::computeFs(int index1, int index2, Vector3d x1, Vector3d x2, double refLen)
{
	indexVec(0) = plateDOF + index1 * rodDOF + 0;
	indexVec(1) = plateDOF + index1 * rodDOF + 1;
	indexVec(2) = plateDOF + index1 * rodDOF + 2;

	indexVec(3) = 3 * index2 + 0;
	indexVec(4) = 3 * index2 + 1;
	indexVec(5) = 3 * index2 + 2;

	epsX = (x2 - x1).norm() / refLen - 1.0;
	tangent = (x2 - x1) / (x2 - x1).norm();
	f = EA * tangent * epsX;

	for (int k = 0; k < 3; k++)
	{
		ind = indexVec(k);
		Globalstepper->addForce(ind, - f[k]); // subtracting elastic force

		ind = indexVec(k+3);
		Globalstepper->addForce(ind, f[k]); // adding elastic force
	}
}

void elasticStretchingForceCouple::computeJs(int index1, int index2, Vector3d x1, Vector3d x2, double refLen)
{
	indexVec(0) = plateDOF + index1 * rodDOF + 0;
	indexVec(1) = plateDOF + index1 * rodDOF + 1;
	indexVec(2) = plateDOF + index1 * rodDOF + 2;

	indexVec(3) = 3 * index2 + 0;
	indexVec(4) = 3 * index2 + 1;
	indexVec(5) = 3 * index2 + 2;

	len = (x2 - x1).norm();
	refLength = refLen;

	dxx(0) = x2(0) - x1(0);
	dxx(1) = x2(1) - x1(1);
	dxx(2) = x2(2) - x1(2);

	u = dxx;
	v = u.transpose();
	M0= EA * ((1/refLength - 1/len) * Id3 + (1/len) * (u*v) / (u.norm() * u.norm()));

	Jss.block(0,0,3,3) =  - M0;
	Jss.block(3,3,3,3) =  - M0;
	Jss.block(3,0,3,3) =    M0;
	Jss.block(0,3,3,3) =    M0;

	for (int j=0;j<6;j++)
	{
		for (int k=0;k<6;k++)
		{
			ind1 = indexVec(j);
			ind2 = indexVec(k);
			Globalstepper->addJacobian(ind1, ind2, - Jss(j, k));
		}
	}
}

void elasticStretchingForceCouple::setFirstJacobian(int index1, int index2)
{
	indexVec(0) = plateDOF + index1 * rodDOF + 0;
	indexVec(1) = plateDOF + index1 * rodDOF + 1;
	indexVec(2) = plateDOF + index1 * rodDOF + 2;

	indexVec(3) = 3 * index2 + 0;
	indexVec(4) = 3 * index2 + 1;
	indexVec(5) = 3 * index2 + 2;

	for (int j=0;j<6;j++)
	{
		for (int k=0;k<6;k++)
		{
			ind1 = indexVec(j);
			ind2 = indexVec(k);
			Globalstepper->addJacobian(ind1, ind2, 1);
		}
	}
}