#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper,
	double m_stiffness, double m_dBar)
{
	plate = &m_plate;
    stepper = &m_stepper;

    stiffness = m_stiffness;
    dBar = m_dBar;

    arrayNum.setZero(9);

    fVec.setZero(9);

    boundaryRodIndex.setZero(plate->boundaryRod.size());
    contactRingIndex.setZero(plate->ringIndex.size());

    Id3<<1,0,0,
         0,1,0,
         0,0,1;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{

	for (int ijk = 0; ijk < 6; ijk++)
	{
		contactRingIndex = plate->ringIndexMatrix.col(ijk);

		boundaryRodIndex = plate->boundaryRodMatrix.col(ijk);

		startPoint = boundaryRodIndex(0);
		endPoint = boundaryRodIndex(plate->boundaryRod.size() - 1);


		for (int i = 0; i < plate->ringIndex.size(); i++)
		{
			nv_0 = contactRingIndex(i);
			x_0 = plate->getVertexOld(nv_0);

			minDistance = 100000.00;
			contactIndex = 0;

			for (int j = startPoint; j < endPoint - 1; j++)
			{
				nv_1 = j;
				nv_2 = j + 1;

				x_1 = plate->getVertexOld(nv_1);
				x_2 = plate->getVertexOld(nv_2);

				e_1 = x_0 - x_1;
				e_2 = x_0 - x_2;
				e_3 = x_2 - x_1;

				e_4 = e_1.cross(e_2);

				d1 = ( x_0 - x_1 ).norm();
				d2 = ( x_0 - x_2 ).norm();
				d3 = e_4.norm() / e_3.norm();

				d = d1 + d2;

				if (d < minDistance)
				{
					minDistance = d;
					contactIndex = j;
				}
			}

			//cout << minDistance << " " << nv_0 << " " << contactIndex << endl;

			nv_0 = contactRingIndex(i);
			x_0 = plate->getVertexOld(nv_0);

			nv_1 = contactIndex;
			nv_2 = contactIndex + 1;

			x_1 = plate->getVertexOld(nv_1);
			x_2 = plate->getVertexOld(nv_2);

			e_1 = x_0 - x_1;
			e_2 = x_0 - x_2;
			e_3 = x_2 - x_1;

			e_4 = e_1.cross(e_2);

			normE1 = e_1.norm();
			normE2 = e_2.norm();
			normE3 = e_3.norm();
			normE4 = e_4.norm();

			t_1 = e_1 / normE1;
			t_2 = e_2 / normE2;
			t_3 = e_3 / normE3;
			t_4 = e_4 / normE4;

			d1 = ( x_0 - x_1 ).norm();
			d2 = ( x_0 - x_2 ).norm();
			d3 = e_4.norm() / e_3.norm();

			//	if (disType == 3)
			{
				arrayNum(0) = 3 * nv_0 + 0;
				arrayNum(1) = 3 * nv_0 + 1;
				arrayNum(2) = 3 * nv_0 + 2;
				arrayNum(3) = 3 * nv_1 + 0;
				arrayNum(4) = 3 * nv_1 + 1;
				arrayNum(5) = 3 * nv_1 + 2;
				arrayNum(6) = 3 * nv_2 + 0;
				arrayNum(7) = 3 * nv_2 + 1;
				arrayNum(8) = 3 * nv_2 + 2;

				d = e_4.norm() / e_3.norm();

				dEnergydD = - 2 * d * log( (-d+dBar) / dBar ) - d * d  / (d - dBar);

				de1 =   crossMat(e_2) * t_4 / normE3;
				de2 = - crossMat(e_1) * t_4 / normE3;
				de3 = - t_3 * normE4 / (normE3 * normE3);

				fVec.segment(0,3) =   de1 + de2;
				fVec.segment(3,3) = - de1 - de3;
				fVec.segment(6,3) = - de2 + de3;

				fVec = - fVec * stiffness * dEnergydD;

				for (int kk = 0; kk < 9; kk++)
				{
					ind = arrayNum(kk);
					stepper->addForce(ind, -fVec[kk]);
				}
			}
		}

	}

	/*
	for (int i = 0; i < plate->ringIndex.size(); i++)
	{
		nv_0 = plate->ringIndex[i];
		x_0 = plate->getVertexOld(nv_0);

		minDistance = 100000.00;
		contactIndex = 0;

		for (int j = plate->edgeStart; j < plate->edgeEnd; j++)
		{
			nv_1 = plate->v_edgeElement[j].nv_1;
			nv_2 = plate->v_edgeElement[j].nv_2;

			x_1 = plate->getVertexOld(nv_1);
			x_2 = plate->getVertexOld(nv_2);

			e_1 = x_0 - x_1;
			e_2 = x_0 - x_2;
			e_3 = x_2 - x_1;

			e_4 = e_1.cross(e_2);

			d1 = ( x_0 - x_1 ).norm();
			d2 = ( x_0 - x_2 ).norm();
			d3 = e_4.norm() / e_3.norm();

			d = d1 + d2;

			if (d < minDistance)
			{
				minDistance = d;
				contactIndex = j;
			}
		}

		//cout << contactIndex << endl;

		nv_0 = plate->ringIndex[i];
		x_0 = plate->getVertexOld(nv_0);

		nv_1 = plate->v_edgeElement[contactIndex].nv_1;
		nv_2 = plate->v_edgeElement[contactIndex].nv_2;

		x_1 = plate->getVertexOld(nv_1);
		x_2 = plate->getVertexOld(nv_2);

		e_1 = x_0 - x_1;
		e_2 = x_0 - x_2;
		e_3 = x_2 - x_1;

		e_4 = e_1.cross(e_2);

		normE1 = e_1.norm();
		normE2 = e_2.norm();
		normE3 = e_3.norm();
		normE4 = e_4.norm();

		t_1 = e_1 / normE1;
		t_2 = e_2 / normE2;
		t_3 = e_3 / normE3;
		t_4 = e_4 / normE4;

		d1 = ( x_0 - x_1 ).norm();
		d2 = ( x_0 - x_2 ).norm();
		d3 = e_4.norm() / e_3.norm();

		//	if (disType == 3)
		{
			arrayNum(0) = 3 * nv_0 + 0;
			arrayNum(1) = 3 * nv_0 + 1;
			arrayNum(2) = 3 * nv_0 + 2;
			arrayNum(3) = 3 * nv_1 + 0;
			arrayNum(4) = 3 * nv_1 + 1;
			arrayNum(5) = 3 * nv_1 + 2;
			arrayNum(6) = 3 * nv_2 + 0;
			arrayNum(7) = 3 * nv_2 + 1;
			arrayNum(8) = 3 * nv_2 + 2;

			d = e_4.norm() / e_3.norm();

			dEnergydD = - 2 * d * log( (-d+dBar) / dBar ) - d * d  / (d - dBar);

			de1 =   crossMat(e_2) * t_4 / normE3;
			de2 = - crossMat(e_1) * t_4 / normE3;
			de3 = - t_3 * normE4 / (normE3 * normE3);

			fVec.segment(0,3) =   de1 + de2;
			fVec.segment(3,3) = - de1 - de3;
			fVec.segment(6,3) = - de2 + de3;

			fVec = - fVec * stiffness * dEnergydD;

			for (int kk = 0; kk < 9; kk++)
			{
				ind = arrayNum(kk);
				stepper->addForce(ind, -fVec[kk]);
			}
		}
	}
	*/
}

void externalContactForce::computeJc()
{
	/*
	for (int i = 0; i < plate->ringIndex.size(); i++)
	{
		nv_0 = plate->ringIndex[i];
		x_0 = plate->getVertex(nv_0);

		minDistance = 100000.00;
		contactIndex = 0;

		for (int j = plate->edgeStart; j < plate->edgeEnd; j++)
		{
			nv_1 = plate->v_edgeElement[j].nv_1;
			nv_2 = plate->v_edgeElement[j].nv_2;

			x_1 = plate->getVertex(nv_1);
			x_2 = plate->getVertex(nv_2);

			e_1 = x_0 - x_1;
			e_2 = x_0 - x_2;
			e_3 = x_2 - x_1;

			e_4 = e_1.cross(e_2);

			d1 = ( x_0 - x_1 ).norm();
			d2 = ( x_0 - x_2 ).norm();
			d3 = e_4.norm() / e_3.norm();

			d = d1 + d2;

			if (d < minDistance)
			{
				minDistance = d;
				contactIndex = j;
			}
		}

		//cout << contactIndex << endl;

		nv_0 = plate->ringIndex[i];
		x_0 = plate->getVertex(nv_0);

		nv_1 = plate->v_edgeElement[contactIndex].nv_1;
		nv_2 = plate->v_edgeElement[contactIndex].nv_2;

		x_1 = plate->getVertex(nv_1);
		x_2 = plate->getVertex(nv_2);

		e_1 = x_0 - x_1;
		e_2 = x_0 - x_2;
		e_3 = x_2 - x_1;

		e_4 = e_1.cross(e_2);

		normE1 = e_1.norm();
		normE2 = e_2.norm();
		normE3 = e_3.norm();
		normE4 = e_4.norm();

		t_1 = e_1 / normE1;
		t_2 = e_2 / normE2;
		t_3 = e_3 / normE3;
		t_4 = e_4 / normE4;

		d1 = ( x_0 - x_1 ).norm();
		d2 = ( x_0 - x_2 ).norm();
		d3 = e_4.norm() / e_3.norm();

		//	if (disType == 3)
		{
			arrayNum(0) = 3 * nv_0 + 0;
			arrayNum(1) = 3 * nv_0 + 1;
			arrayNum(2) = 3 * nv_0 + 2;
			arrayNum(3) = 3 * nv_1 + 0;
			arrayNum(4) = 3 * nv_1 + 1;
			arrayNum(5) = 3 * nv_1 + 2;
			arrayNum(6) = 3 * nv_2 + 0;
			arrayNum(7) = 3 * nv_2 + 1;
			arrayNum(8) = 3 * nv_2 + 2;

			d = e_4.norm() / e_3.norm();

			gradT3 = (Id3 - t_3 * t_3.transpose()) / normE3;
			gradT4 = (Id3 - t_4 * t_4.transpose()) / normE4;

			de1de1 = - crossMat(e_2) * gradT4 * crossMat(e_2) / normE3;
			de2de2 = - crossMat(e_1) * gradT4 * crossMat(e_1) / normE3;
			de3de3 = - (gradT3 * normE3 - 2 * t_3 * t_3.transpose()) * normE4 / (normE3 * normE3 * normE3);

			de1de2 = (   crossMat(t_4) + crossMat(e_2) * gradT4 * crossMat(e_1) ) / normE3;
			de1de3 = ( - crossMat(e_2) * t_4 * t_3.transpose() / ( normE3 * normE3) ).transpose();
			de2de3 = (   crossMat(e_1) * t_4 * t_3.transpose() / ( normE3 * normE3) ).transpose();

			de2de1 = de1de2.transpose();
			de3de1 = de1de3.transpose();
			de3de2 = de2de3.transpose();

			de1 =   crossMat(e_2) * t_4 / normE3;
			de2 = - crossMat(e_1) * t_4 / normE3;
			de3 = - t_3 * normE4 / (normE3 * normE3);

			dx.segment(0,3) =   de1 + de2;
			dx.segment(3,3) = - de1 - de3;
			dx.segment(6,3) = - de2 + de3;

			dxdx.block(0,0,3,3) = de1de1 + de2de2 + de1de2 + de2de1;
			dxdx.block(3,3,3,3) = de1de1 + de3de3 + de1de3 + de3de1;
			dxdx.block(6,6,3,3) = de2de2 + de3de3 - de2de3 - de3de2;

			dxdx.block(0,3,3,3) = - de1de1 - de1de3 - de2de1 - de2de3;
        	dxdx.block(0,6,3,3) = - de2de2 + de1de3 + de2de3 - de1de2;
        	dxdx.block(3,6,3,3) = - de3de3 - de1de3 + de1de2 + de3de2;

        	dxdx.block(3,0,3,3) = dxdx.block(0,3,3,3).transpose();
        	dxdx.block(6,0,3,3) = dxdx.block(0,6,3,3).transpose();
        	dxdx.block(6,3,3,3) = dxdx.block(3,6,3,3).transpose();

        	dEnergydD = - 2 * d * log( (-d+dBar) / dBar ) - d * d  / (d - dBar);
        	d2EnergydD2 = d * d / ((-d+dBar) * (-d+dBar)) + 4 * d / (-d+dBar) - 2 * log( (-d+dBar) / dBar );

        	jacobian = stiffness * (d2EnergydD2 * dx * dx.transpose() + dEnergydD * dxdx);

			for (int jj = 0; jj < 9; jj++)
        	{
            	for (int kk = 0; kk < 9; kk++)
            	{
					ind1 = arrayNum(jj);
					ind2 = arrayNum(kk);

					stepper->addJacobian(ind1, ind2, jacobian(jj,kk));
            	}
        	}
		}
	}
	*/
}

void externalContactForce::setFirstJacobian()
{
	/*
	for (int i = 0; i < plate->ringIndex.size(); i++)
	{
		for (int j = plate->edgeStart; j < plate->edgeEnd; j++)
		{
			nv_0 = plate->ringIndex[i];

			nv_1 = plate->v_edgeElement[j].nv_1;
			nv_2 = plate->v_edgeElement[j].nv_2;

			arrayNum(0) = 3 * nv_0 + 0;
			arrayNum(1) = 3 * nv_0 + 1;
			arrayNum(2) = 3 * nv_0 + 2;
			arrayNum(3) = 3 * nv_1 + 0;
			arrayNum(4) = 3 * nv_1 + 1;
			arrayNum(5) = 3 * nv_1 + 2;
			arrayNum(6) = 3 * nv_2 + 0;
			arrayNum(7) = 3 * nv_2 + 1;
			arrayNum(8) = 3 * nv_2 + 2;

			for (int j = 0; j < 9; j++)
        	{
            	for (int k = 0; k < 9; k++)
            	{
					ind1 = arrayNum(j);
					ind2 = arrayNum(k);

					stepper->addJacobian(ind1, ind2, 1);
            	}
        	}
		}
	}
	*/
}

Matrix3d externalContactForce::crossMat(Vector3d a)
{
	Matrix3d b;

	b<<0,-a(2),a(1),
	a(2),0,-a(0),
	-a(1),a(0),0;

	return b;
}
