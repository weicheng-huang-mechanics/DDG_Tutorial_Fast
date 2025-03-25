#include "elasticBendingForce.h"

elasticBendingForce::elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
    plate = &m_plate;
    stepper = &m_stepper;

    Jbb = MatrixXd::Zero(6,6);
    flocal = VectorXd::Zero(6);

    localDOF = VectorXi::Zero(6);
}

elasticBendingForce::~elasticBendingForce()
{
	;
}

void elasticBendingForce::computeFb()
{
    for(int k = 0; k < plate->bendingNum; k++)
    {
        flocal = VectorXd::Zero(6);

        ind0 = plate->v_bendingElement[k].nv_1;
        ind1 = plate->v_bendingElement[k].nv_2;
        ind2 = plate->v_bendingElement[k].nv_3;

        p0 = plate->v_bendingElement[k].x_1;
        p  = plate->v_bendingElement[k].x_2;
        p1 = plate->v_bendingElement[k].x_3;

        xkm1 = p0[0];
        ykm1 = p0[1];
        xk = p[0];
        yk = p[1];
        xkp1 = p1[0];
        ykp1 = p1[1];

        l_k = plate->v_bendingElement[k].voroniLength;
        curvature0 = plate->v_bendingElement[k].nBar;
 
        EI_local = plate->v_bendingElement[k].EI_local;

        flocal = - EI_local * computeBendingForce(xkm1, ykm1, xk, yk, xkp1, ykp1, l_k, curvature0);

        localDOF = VectorXi::Zero(6);

        localDOF(0) = 2 * ind0 + 0;
        localDOF(1) = 2 * ind0 + 1;
        localDOF(2) = 2 * ind1 + 0;
        localDOF(3) = 2 * ind1 + 1;
        localDOF(4) = 2 * ind2 + 0;
        localDOF(5) = 2 * ind2 + 1;

        for (int i = 0; i < 6; i++)
        {
            stepper->addForce(localDOF(i), - flocal(i));
        }
    }
}

void elasticBendingForce::computeJb()
{
    for(int k=0; k < plate->bendingNum; k++)
    {
        Jbb = MatrixXd::Zero(6,6);

        ind0 = plate->v_bendingElement[k].nv_1;
        ind1 = plate->v_bendingElement[k].nv_2;
        ind2 = plate->v_bendingElement[k].nv_3;

        p0 = plate->v_bendingElement[k].x_1;
        p  = plate->v_bendingElement[k].x_2;
        p1 = plate->v_bendingElement[k].x_3;

        xkm1 = p0[0];
        ykm1 = p0[1];
        xk = p[0];
        yk = p[1];
        xkp1 = p1[0];
        ykp1 = p1[1];

        l_k = plate->v_bendingElement[k].voroniLength;
        curvature0 = plate->v_bendingElement[k].nBar;

        EI_local = plate->v_bendingElement[k].EI_local;

        Jbb = EI_local * computeBendingJacobian(xkm1, ykm1, xk, yk, xkp1, ykp1, l_k, curvature0);

        localDOF = VectorXi::Zero(6);

        localDOF(0) = 2 * ind0 + 0;
        localDOF(1) = 2 * ind0 + 1;
        localDOF(2) = 2 * ind1 + 0;
        localDOF(3) = 2 * ind1 + 1;
        localDOF(4) = 2 * ind2 + 0;
        localDOF(5) = 2 * ind2 + 1;

        for (int i = 0; i < 6; i++)
        {
            for(int j = 0; j < 6; j++)
            {
                stepper->addJacobian(localDOF(i), localDOF(j), Jbb(i,j));
            }
        }
    }
}

void elasticBendingForce::setFirstJacobian()
{
    for(int k=0; k < plate->bendingNum; k++)
    {
        ind0 = plate->v_bendingElement[k].nv_1;
        ind1 = plate->v_bendingElement[k].nv_2;
        ind2 = plate->v_bendingElement[k].nv_3;

        localDOF = VectorXi::Zero(6);

        localDOF(0) = 2 * ind0 + 0;
        localDOF(1) = 2 * ind0 + 1;
        localDOF(2) = 2 * ind1 + 0;
        localDOF(3) = 2 * ind1 + 1;
        localDOF(4) = 2 * ind2 + 0;
        localDOF(5) = 2 * ind2 + 1;

        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                stepper->addJacobian(localDOF(i), localDOF(j), 1);
            }
        }
    }
}

VectorXd elasticBendingForce::computeBendingForce(double xa, double ya, double xi, double yi, double xb, double yb, 
    double lBar, double kappaBar)
{
    VectorXd vecResult;

    vecResult = ListVec(1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
           ((-xa + xi)*(yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
           (xb - xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
         (((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
      (2*((pow(-xa + xi,2)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
           ((xb - xi)*(-xa + xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
           (yb - yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
    (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
   1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
           ((yb - yi)*pow(-ya + yi,2))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
           (yb - yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
         (((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
      (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
           ((xb - xi)*pow(-ya + yi,2))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
           (xb - xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
    (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
   1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
         (-(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
           ((-xa + xi)*(yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
           (pow(xb - xi,2)*(-xa + xi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           (xb - xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           (-xa + xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((xb - xi)*(yb - yi)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
      (2*(-((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
           ((xb - xi)*(-xa + xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
           ((xb - xi)*(-xa + xi)*(yb - yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           (yb - yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           (pow(xb - xi,2)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           (-ya + yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
    (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
   1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
    ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
           ((xb - xi)*pow(-ya + yi,2))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
           (xb - xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           (-xa + xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((-xa + xi)*pow(yb - yi,2))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(yb - yi)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
      (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
         (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
           ((yb - yi)*pow(-ya + yi,2))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
           ((xb - xi)*(-xa + xi)*(yb - yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           (yb - yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           (-ya + yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           (pow(yb - yi,2)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
   1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
         (-((pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
           (-xa + xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(yb - yi)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
      (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
           (pow(xb - xi,2)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           (-ya + yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
    (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
   1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
    ((2*((-xa + xi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((-xa + xi)*pow(yb - yi,2))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((xb - xi)*(yb - yi)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
      (2*(((-xa + xi)*(yb - yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           ((xb - xi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
         (-(((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
           (-ya + yi)/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
           (pow(yb - yi,2)*(-ya + yi))/
            (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
       (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
           ((yb - yi)*(-ya + yi))/
            (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
              sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))));

    return vecResult;
}

MatrixXd elasticBendingForce::computeBendingJacobian(double xa, double ya, double xi, double yi, double xb, double yb, 
    double lBar, double kappaBar)
{
    MatrixXd matResult;

    matResult = ListMat(ListVec(1.*lBar*pow((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))),2) + 
     1.*lBar*((4*pow(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((3*(xb - xi)*pow(-xa + xi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*pow(-xa + xi,3)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((-3*(xb - xi)*pow(-xa + xi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (pow(xb - xi,2)*pow(-xa + xi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*pow(-xa + xi,3)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*((-3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-((pow(xb - xi,2)*pow(-xa + xi,2))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))),
   ListVec(1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*pow((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))),2) + 
     1.*lBar*((4*pow(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (4*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(yb - yi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((-3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*((-3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(yb - yi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(yb - yi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*(xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-((pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(yb - yi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))),
   ListVec(1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((-3*(xb - xi)*pow(-xa + xi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (pow(xb - xi,2)*pow(-xa + xi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*pow(-xa + xi,3)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((-3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*pow((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))),2) + 
     1.*lBar*((4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           pow(-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((3*(xb - xi)*pow(-xa + xi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (2*pow(xb - xi,2)*pow(-xa + xi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*(xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (2*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (2*(xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*pow(xb - xi,3)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (2*pow(xb - xi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             2/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (4*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*pow(-xa + xi,3)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (2*(xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*(-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (2*pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (2*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (2*(xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (2*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (2*(xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((pow(xb - xi,2)*pow(-xa + xi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,3)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(pow(-xa + xi,2)/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))),
   ListVec(1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*((-3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*((-3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(yb - yi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(yb - yi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((-3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(xb - xi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*(xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*pow(-xa + xi,2)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((3*(xb - xi)*pow(-xa + xi,2)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (2*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (2*(xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*pow((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)),2) + 
     1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((3*(-xa + xi)*(yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             (3*(xb - xi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (3*(xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (2*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (2*(-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (2*(xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (2*(xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(-xa + xi)*pow(yb - yi,3))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (4*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           pow(-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((3*(xb - xi)*(-xa + xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) + 
             (3*(yb - yi)*pow(-ya + yi,3))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),2.5)) - 
             ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (2*(xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (2*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (2*pow(yb - yi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             2/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (2*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(yb - yi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(-xa + xi)*pow(yb - yi,3))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(yb - yi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(yb - yi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))),
   ListVec(1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-((pow(xb - xi,2)*pow(-xa + xi,2))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-((pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))) + 
     1.*lBar*((4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((pow(xb - xi,2)*pow(-xa + xi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,3)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(xb - xi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             pow(xb - xi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((pow(xb - xi,2)*(-xa + xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*pow((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))),2) + 
     1.*lBar*((4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           pow(-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((3*pow(xb - xi,3)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (4*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      (-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))),
   ListVec(1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(pow(-xa + xi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*pow(-xa + xi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*((pow(-xa + xi,2)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (pow(yb - yi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) + 
        (4*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(pow(-xa + xi,2)/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             (pow(-xa + xi,2)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*pow(-xa + xi,2)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*pow(-xa + xi,2))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((-xa + xi)*(yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(xb - xi,2)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-((pow(-xa + xi,2)*(yb - yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*(-xa + xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))*
      ((2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*(-(((-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((-xa + xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             ((xb - xi)*(yb - yi)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(-xa + xi)*pow(yb - yi,3))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(-(((-xa + xi)*(yb - yi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) + 
             ((xb - xi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             (xb - xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) - 
        (2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(-ya + yi))/
                (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                  pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5))) - 
             ((yb - yi)*pow(-ya + yi,2))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (yb - yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (((xb - xi)*(-xa + xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) - 
             pow(-ya + yi,2)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             (pow(yb - yi,2)*pow(-ya + yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                pow(pow(-xa + xi,2) + pow(-ya + yi,2),1.5)) + 
             ((xb - xi)*(-xa + xi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             1/(sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             pow(yb - yi,2)/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(yb - yi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((-2*(-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*pow(xb - xi,2)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           ((3*pow(xb - xi,2)*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))) + 
     1.*lBar*((-2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-((pow(xb - xi,2)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (2*(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (pow(xb - xi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))),
    1.*lBar*pow((2*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)),2) + 
     1.*lBar*(-kappaBar + (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))))*
      ((2*((-3*(-xa + xi)*(yb - yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*(-xa + xi)*pow(yb - yi,3))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(xb - xi)*pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))) - 
        (4*((-xa + xi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((xb - xi)*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)) + 
        (4*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           pow(-(((xb - xi)*(-xa + xi)*(yb - yi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (-ya + yi)/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (pow(yb - yi,2)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),3)) - 
        (2*(((-xa + xi)*(yb - yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             ((xb - xi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))))*
           (-(((xb - xi)*(-xa + xi))/
                (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                  sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))) + 
             (3*(xb - xi)*(-xa + xi)*pow(yb - yi,2))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) - 
             (3*(yb - yi)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),1.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             (3*pow(yb - yi,3)*(-ya + yi))/
              (pow(pow(xb - xi,2) + pow(yb - yi,2),2.5)*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2)))))/
         (lBar*pow(1 + ((xb - xi)*(-xa + xi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))) + 
             ((yb - yi)*(-ya + yi))/
              (sqrt(pow(xb - xi,2) + pow(yb - yi,2))*
                sqrt(pow(-xa + xi,2) + pow(-ya + yi,2))),2)))));

    return matResult;
}

VectorXd elasticBendingForce::ListVec(double a1, double a2, double a3, 
    double a4, double a5, double a6)
{
    VectorXd vecResult;

    vecResult.setZero(6, 1);

    vecResult(0) = a1;
    vecResult(1) = a2;
    vecResult(2) = a3;
    vecResult(3) = a4;
    vecResult(4) = a5;
    vecResult(5) = a6;

    return vecResult;
}

MatrixXd elasticBendingForce::ListMat(VectorXd a1, VectorXd a2, VectorXd a3, 
    VectorXd a4, VectorXd a5, VectorXd a6)
{
    MatrixXd matResult;

    matResult.setZero(6, 6);

    matResult.col(0) = a1;
    matResult.col(1) = a2;
    matResult.col(2) = a3;
    matResult.col(3) = a4;
    matResult.col(4) = a5;
    matResult.col(5) = a6;

    return matResult;
}