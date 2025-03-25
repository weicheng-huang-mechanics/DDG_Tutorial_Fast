#include "elasticStretchingForce.h"
#include <iostream>

elasticStretchingForce::elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;
	
	Jss.setZero(4, 4);
	flocal = VectorXd::Zero(4);
	EA = plate->EA;
	mu = plate->Possion;

	localDOF = VectorXi::Zero(4);
}

elasticStretchingForce::~elasticStretchingForce()
{
	;
}

void elasticStretchingForce::computeFs()
{
	totalForce = VectorXd::Zero(plate->ndof);

	for (int k = 0; k < plate->edgeNum; k++)
	{
		flocal = VectorXd::Zero(4);

		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		p = plate->getVertex(ind1);
		p1 = plate->getVertex(ind2);

		xk = p[0];
		yk = p[1];
		xkp1 = p1[0];
		ykp1 = p1[1];

		l_k = plate->v_edgeElement[k].refLength;

		rBar1 = plate->v_edgeElement[k].x_1_start(0);
		rBar2 = plate->v_edgeElement[k].x_2_start(0);

		flocal = computeStretchingForce(xk, yk, xkp1, ykp1, mu, l_k, rBar1, rBar2);

		flocal = (- 0.5 * EA * l_k * M_PI * (rBar1 + rBar2) ) * flocal;

		localDOF = VectorXi::Zero(4);

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			stepper->addForce(localDOF(i), - flocal(i));

			totalForce(localDOF(i)) = totalForce(localDOF(i)) + flocal(i);
		}
	}
}

void elasticStretchingForce::computeJs()
{
	for (int k = 0; k < plate->edgeNum; k++)
	{
		Jss.setZero(4, 4);

		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		p = plate->getVertex(ind1);
		p1 = plate->getVertex(ind2);

		xk = p[0];
		yk = p[1];
		xkp1 = p1[0];
		ykp1 = p1[1];

		l_k = plate->v_edgeElement[k].refLength;

		rBar1 = plate->v_edgeElement[k].x_1_start(0);
		rBar2 = plate->v_edgeElement[k].x_2_start(0);

		Jss = computeStretchingJacobian(xk, yk, xkp1, ykp1, mu, l_k, rBar1, rBar2);
		Jss = ( 0.5 * EA * l_k * M_PI * (rBar1 + rBar2) ) * Jss; // scale with stiffness

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

VectorXd elasticStretchingForce::computeStretchingForce(double xa, double ya, double xb, double yb, 
  double mu, double lBar, double rBar1, double rBar2)
{
  VectorXd vecResult;

  vecResult = ListVec(((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)/(2.*rBar1) - 
    (mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
     (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
    (mu*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/(lBar*rBar1) - 
    (2*(-xa + xb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
     (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
   -((mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
       (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2)))) - 
    (2*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
     (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
   ((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)/(2.*rBar2) + 
    (mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
     (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
    (mu*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/(lBar*rBar2) + 
    (2*(-xa + xb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
     (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
   (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
     (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
    (2*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
     (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))));

  return vecResult;
}


MatrixXd elasticStretchingForce::computeStretchingJacobian(double xa, double ya, double xb, double yb, 
  double mu, double lBar, double rBar1, double rBar2)
{
  MatrixXd matResult;

  matResult = ListMat(ListVec(1/(2.*pow(rBar1,2)) - (mu*pow(-xa + xb,2)*
        ((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) + 
     (2*pow(-xa + xb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*mu*(-xa + xb))/(lBar*rBar1*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*pow(-xa + xb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) + 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
    -((mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
        (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))) + 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*(-ya + yb))/(lBar*rBar1*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)),
    1/(2.*rBar1*rBar2) + (mu*pow(-xa + xb,2)*
        ((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*pow(-xa + xb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*(-xa + xb))/(lBar*rBar1*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*(-xa + xb))/(lBar*rBar2*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*pow(-xa + xb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
    (mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*(-ya + yb))/(lBar*rBar1*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))),
   ListVec(-((mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
        (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))) + 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*(-ya + yb))/(lBar*rBar1*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)),
    -((mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*pow(-ya + yb,2))/
        (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))) + 
     (2*pow(-ya + yb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*pow(-ya + yb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) + 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
    (mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*(-ya + yb))/(lBar*rBar2*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)),
    (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*pow(-ya + yb,2))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*pow(-ya + yb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*pow(-ya + yb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2)))),
   ListVec(1/(2.*rBar1*rBar2) + (mu*pow(-xa + xb,2)*
        ((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*pow(-xa + xb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*(-xa + xb))/(lBar*rBar1*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*(-xa + xb))/(lBar*rBar2*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*pow(-xa + xb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
    (mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*(-ya + yb))/(lBar*rBar2*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)),
    1/(2.*pow(rBar2,2)) - (mu*pow(-xa + xb,2)*
        ((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) + 
     (2*pow(-xa + xb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*mu*(-xa + xb))/(lBar*rBar2*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*pow(-xa + xb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) + 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
    -((mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
        (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))) + 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*(-ya + yb))/(lBar*rBar2*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))),
   ListVec((mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*(-ya + yb))/(lBar*rBar1*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)),
    (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*pow(-ya + yb,2))/
      (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*pow(-ya + yb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (2*pow(-ya + yb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) - 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))),
    -((mu*(-xa + xb)*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*(-ya + yb))/
        (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))) + 
     (2*(-xa + xb)*(-ya + yb))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*(-ya + yb))/(lBar*rBar2*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*(-xa + xb)*(-ya + yb)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)),
    -((mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2)*pow(-ya + yb,2))/
        (lBar*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5))) + 
     (2*pow(-ya + yb,2))/(pow(lBar,2)*(pow(-xa + xb,2) + pow(-ya + yb,2))) + 
     (mu*((-rBar1 + xa)/rBar1 + (-rBar2 + xb)/rBar2))/
      (lBar*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))) - 
     (2*pow(-ya + yb,2)*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*pow(pow(-xa + xb,2) + pow(-ya + yb,2),1.5)) + 
     (2*(-lBar + sqrt(pow(-xa + xb,2) + pow(-ya + yb,2))))/
      (pow(lBar,2)*sqrt(pow(-xa + xb,2) + pow(-ya + yb,2)))));

  return matResult;
}

VectorXd elasticStretchingForce::ListVec(double a1, double a2, double a3, double a4)
{
  VectorXd vecResult;

  vecResult.setZero(4, 1);

  vecResult(0) = a1;
  vecResult(1) = a2;
  vecResult(2) = a3;
  vecResult(3) = a4;

  return vecResult;
}

MatrixXd elasticStretchingForce::ListMat(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4)
{
  MatrixXd matResult;

  matResult.setZero(4, 4);

  matResult.col(0) = a1;
  matResult.col(1) = a2;
  matResult.col(2) = a3;
  matResult.col(3) = a4;

  return matResult;
}