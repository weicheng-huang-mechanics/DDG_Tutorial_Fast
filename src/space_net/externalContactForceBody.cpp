#include "externalContactForceBody.h"

externalContactForceBody::externalContactForceBody(elasticPlate &m_plate, timeStepper &m_stepper, double m_boxSize,
	double m_stiffness, double m_dBar)
{
	plate = &m_plate;
    stepper = &m_stepper;

    boxSize = m_boxSize;
    stiffness = m_stiffness;
    dBar = m_dBar;

    indexArray = VectorXi::Zero(6);
    jacobian = MatrixXd::Zero(6, 6);

    Id3<<1,0,0,
         0,1,0,
         0,0,1;
}

externalContactForceBody::~externalContactForceBody()
{
	;
}

void externalContactForceBody::computeFc()
{
	Vector3d xEnd = plate->getVertex(plate->nv - 1);

	for(int i = 0; i < plate->nv - 1; i++)
	{
		Vector3d xCurrent = plate->getVertex(i);

		double d = (xEnd - xCurrent).norm();

		if (d <= dBar)
		{
			dDdEdge = (xEnd - xCurrent) / (xEnd - xCurrent).norm();

			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			f = stiffness * dEnergydD * dDdEdge;

			for (int k = 0; k < 3; k++)
			{
				ind = 3 * i + k;
				stepper->addForce(ind, - f[k]);
			}

			for (int k = 0; k < 3; k++)
			{
				ind = 3 * (plate->nv-1) + k;
				stepper->addForce(ind, f[k]);
			}
		}
	}
}

void externalContactForceBody::computeJc()
{
	Vector3d xEnd = plate->getVertex(plate->nv - 1);

	for(int i = 0; i < plate->nv - 1; i++)
	{
		indexArray(0) = 3 * (plate->nv - 1) + 0;
		indexArray(1) = 3 * (plate->nv - 1) + 1;
		indexArray(2) = 3 * (plate->nv - 1) + 2;

		indexArray(3) = 3 * i + 0;
		indexArray(4) = 3 * i + 1;
		indexArray(5) = 3 * i + 2;

		Vector3d xCurrent = plate->getVertex(i);

		double d = (xEnd - xCurrent).norm();

		if (d <= dBar)
		{
			dDdEdge = (xEnd - xCurrent) / (xEnd - xCurrent).norm();

			d2DdEdge2 = (Id3 - dDdEdge * dDdEdge.transpose()) / (xEnd - xCurrent).norm();

			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			d2EnergydD2 = - 2 * log(d / dBar) - 2 * (d - dBar) / d - 2 * (d - dBar) / d + (d - dBar) * (d - dBar) / (d * d);

			d2EdEdge2 = stiffness * ( d2EnergydD2 * dDdEdge * dDdEdge.transpose() + dEnergydD * d2DdEdge2 );

			jacobian.block(0,0,3,3) =   d2EdEdge2;
			jacobian.block(3,3,3,3) =   d2EdEdge2;
			jacobian.block(3,0,3,3) = - d2EdEdge2;
			jacobian.block(0,3,3,3) = - d2EdEdge2;

			for (int j = 0; j < 6; j++)
			{
				for (int k = 0; k < 6; k++)
				{
					stepper->addJacobian(indexArray(j), indexArray(k), jacobian(k,j));
				}
			}

		}
	}
}

void externalContactForceBody::setFirstJacobian()
{
	for(int i = 0; i < plate->nv - 1; i++)
	{
		indexArray(0) = 3 * (plate->nv - 1) + 0;
		indexArray(1) = 3 * (plate->nv - 1) + 1;
		indexArray(2) = 3 * (plate->nv - 1) + 2;

		indexArray(3) = 3 * i + 0;
		indexArray(4) = 3 * i + 1;
		indexArray(5) = 3 * i + 2;

		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				stepper->addJacobian(indexArray(j), indexArray(k), 1);
			}
		}
	}
}