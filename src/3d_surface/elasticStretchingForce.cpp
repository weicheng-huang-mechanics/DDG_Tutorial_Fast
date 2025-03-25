#include "elasticStretchingForce.h"
#include <iostream>

elasticStretchingForce::elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;
	
	f.setZero(3);
	Jss.setZero(6,6);

	arrayNum.setZero(6);

	Id3<<1,0,0,
	     0,1,0,
	     0,0,1;

	TotalForceVec = VectorXd::Zero(plate->ndof);
}

elasticStretchingForce::~elasticStretchingForce()
{
	;
}

void elasticStretchingForce::computeFs()
{
	TotalForceVec = VectorXd::Zero(plate->ndof);

	sumStretch = 0.0;

	for (int i = 0; i < plate->edgeNum; i++)
	{
		arrayNum = plate->v_edgeElement[i].arrayNum;

		x_1 = plate->v_edgeElement[i].x_1;
		x_2 = plate->v_edgeElement[i].x_2;

		len = ( x_2 - x_1 ).norm();

		refLength = plate->v_edgeElement[i].refLength;

		epsX = len / refLength - 1.0;
		
		tangentVec = ( x_2 - x_1 ) / len;

		f = plate->EA * refLength * tangentVec * epsX;

		sumStretch = sumStretch + epsX;

		for (int k = 0; k < 3; k++)
		{
			ind = arrayNum(k);

			TotalForceVec(ind) = TotalForceVec(ind) + f[k];
			stepper->addForce(ind, -f[k]);

			ind = arrayNum(k+3);

			TotalForceVec(ind) = TotalForceVec(ind) - f[k];
			stepper->addForce(ind, f[k]);
		}
	}

	sumStretch = sumStretch / plate->edgeNum;
}

void elasticStretchingForce::computeJs()
{
	for (int i = 0; i < plate->edgeNum; i++)
	{
		arrayNum = plate->v_edgeElement[i].arrayNum;

		x_1 = plate->v_edgeElement[i].x_1;
		x_2 = plate->v_edgeElement[i].x_2;

		len = ( x_2 - x_1 ).norm();

		refLength = plate->v_edgeElement[i].refLength;	

		dxx(0) = x_2(0) - x_1(0);
		dxx(1) = x_2(1) - x_1(1);
		dxx(2) = x_2(2) - x_1(2);	

		u = dxx;
		v = u.transpose();
		M0= plate->EA * refLength * ((1/refLength - 1/len) * Id3 + (1/len) * (u*v) / (u.norm() * u.norm()));

		Jss.block(0,0,3,3) =  - M0;
		Jss.block(3,3,3,3) =  - M0;
		Jss.block(3,0,3,3) =    M0;
		Jss.block(0,3,3,3) =    M0;

		for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
				ind1 = arrayNum(j);
				ind2 = arrayNum(k);

				stepper->addJacobian(ind1, ind2, - Jss(k,j));
            }
        }
	}
}

void elasticStretchingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->edgeNum; i++)
	{
		arrayNum = plate->v_edgeElement[i].arrayNum;

		for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
				ind1 = arrayNum(j);
				ind2 = arrayNum(k);

				stepper->addJacobian(ind1, ind2, 1);
            }
        }
	}
}
