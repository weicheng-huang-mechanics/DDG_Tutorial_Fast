#include "elasticStretchingForce.h"
#include <iostream>

elasticStretchingForce::elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;
	
	Jss.setZero(4, 4);
	flocal = VectorXd::Zero(4);
	EA = plate->EA;

	localDOF = VectorXi::Zero(4);
}

elasticStretchingForce::~elasticStretchingForce()
{
	;
}

void elasticStretchingForce::computeFs()
{
	for (int k = 0; k < plate->edgeNum; k++)
	{
		flocal = VectorXd::Zero(4);

		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		p = plate->v_edgeElement[k].x_1;
		p1 = plate->v_edgeElement[k].x_2;

		xk = p[0];
		yk = p[1];
		xkp1 = p1[0];
		ykp1 = p1[1];

		l_k = plate->v_edgeElement[k].refLength;

		localForce(); // populate flocal: first derivative of (axial stretch)^2
		flocal = (- 0.5 * EA * l_k) * flocal; // scale with stiffness

		localDOF = VectorXi::Zero(4);

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			stepper->addForce(localDOF(i), - flocal(i));
		}
	}
}

void elasticStretchingForce::computeJs()
{
	for (int k = 0; k < plate->edgeNum; k++)
	{
		Jss.setZero(4,4);

		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		p = plate->v_edgeElement[k].x_1;
		p1 = plate->v_edgeElement[k].x_2;

		xk = p[0];
		yk = p[1];
		xkp1 = p1[0];
		ykp1 = p1[1];

		l_k = plate->v_edgeElement[k].refLength;

		localJacobian();
		Jss = (0.5 * EA * l_k) * Jss; // scale with stiffness

		localDOF = VectorXi::Zero(4);

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				stepper->addJacobian(localDOF(i), localDOF(j), Jss(i,j));
			}
		}
	}
}

void elasticStretchingForce::setFirstJacobian()
{
	for (int k = 0; k < plate->edgeNum; k++)
	{
		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		localDOF = VectorXi::Zero(4);

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				stepper->addJacobian(localDOF(i), localDOF(j), 1);
			}
		}
	}
}

void elasticStretchingForce::localForce()
{
	flocal(0) = -(0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk);
	flocal(1) = -(0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k * (-0.2e1 * ykp1 + 0.2e1 * yk);
	flocal(2) = -(0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k * (0.2e1 * xkp1 - 0.2e1 * xk);
	flocal(3) = -(0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k * (0.2e1 * ykp1 - 0.2e1 * yk);
}

void elasticStretchingForce::localJacobian()
{
	Jss(0,0) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * pow(-0.2e1 * xkp1 + 0.2e1 * xk, 0.2e1) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * pow(-0.2e1 * xkp1 + 0.2e1 * xk, 0.2e1) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k;
	Jss(0,1) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * (-0.2e1 * ykp1 + 0.2e1 * yk) * (-0.2e1 * xkp1 + 0.2e1 * xk) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk) * (-0.2e1 * ykp1 + 0.2e1 * yk) / 0.2e1;
	Jss(0,2) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * (0.2e1 * xkp1 - 0.2e1 * xk) * (-0.2e1 * xkp1 + 0.2e1 * xk) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk) * (0.2e1 * xkp1 - 0.2e1 * xk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k;
	Jss(0,3) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * (0.2e1 * ykp1 - 0.2e1 * yk) * (-0.2e1 * xkp1 + 0.2e1 * xk) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk) * (0.2e1 * ykp1 - 0.2e1 * yk) / 0.2e1;
	Jss(1,0) = Jss(0,1);
	Jss(1,1) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * pow(-0.2e1 * ykp1 + 0.2e1 * yk, 0.2e1) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * pow(-0.2e1 * ykp1 + 0.2e1 * yk, 0.2e1) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k;
	Jss(1,2) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * (0.2e1 * xkp1 - 0.2e1 * xk) * (-0.2e1 * ykp1 + 0.2e1 * yk) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * (-0.2e1 * ykp1 + 0.2e1 * yk) * (0.2e1 * xkp1 - 0.2e1 * xk) / 0.2e1;
	Jss(1,3) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * (0.2e1 * ykp1 - 0.2e1 * yk) * (-0.2e1 * ykp1 + 0.2e1 * yk) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * (-0.2e1 * ykp1 + 0.2e1 * yk) * (0.2e1 * ykp1 - 0.2e1 * yk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k;
	Jss(2,0) = Jss(0,2);
	Jss(2,1) = Jss(1,2);
	Jss(2,2) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * pow(0.2e1 * xkp1 - 0.2e1 * xk, 0.2e1) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * pow(0.2e1 * xkp1 - 0.2e1 * xk, 0.2e1) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k;
	Jss(2,3) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * (0.2e1 * ykp1 - 0.2e1 * yk) * (0.2e1 * xkp1 - 0.2e1 * xk) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * (0.2e1 * xkp1 - 0.2e1 * xk) * (0.2e1 * ykp1 - 0.2e1 * yk) / 0.2e1;
	Jss(3,0) = Jss(0,3);
	Jss(3,1) = Jss(1,3);
	Jss(3,2) = Jss(2,3);
	Jss(3,3) = 0.1e1 / (pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) * pow(l_k, -0.2e1) * pow(0.2e1 * ykp1 - 0.2e1 * yk, 0.2e1) / 0.2e1 + (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.3e1 / 0.2e1) / l_k * pow(0.2e1 * ykp1 - 0.2e1 * yk, 0.2e1) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1)) / l_k) * pow(pow(xkp1 - xk, 0.2e1) + pow(ykp1 - yk, 0.2e1), -0.1e1 / 0.2e1) / l_k;
}