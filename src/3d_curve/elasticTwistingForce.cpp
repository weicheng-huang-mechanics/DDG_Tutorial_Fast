#include "elasticTwistingForce.h"

elasticTwistingForce::elasticTwistingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;
	
    gradTwist = MatrixXd::Zero(plate->bendingNum, 11);

    DDtwist.setZero(11,11);
    Jtt.setZero(11,11);
    gradTwistLocal.setZero(11);
    f.setZero(11);

    GJ = plate->GJ;    

    localDOF = VectorXi::Zero(11);
}

elasticTwistingForce::~elasticTwistingForce()
{
	;
}

void elasticTwistingForce::computeFt()
{
    for (int i = 0; i < plate->bendingNum; i++)
    {
        deltam = plate->v_bendingElement[i].theta_2 - plate->v_bendingElement[i].theta_1;

        norm_e = plate->v_bendingElement[i].norm_1;
        norm_f = plate->v_bendingElement[i].norm_2;

        gradTwist.row(i).segment(0,3) = -0.5 / norm_e * plate->v_bendingElement[i].kb;
        gradTwist.row(i).segment(8,3) = 0.5 / norm_f * plate->v_bendingElement[i].kb;
        gradTwist.row(i).segment(4,3) = -(gradTwist.row(i).segment(0,3)+gradTwist.row(i).segment(8,3));
        gradTwist(i, 3) = -1;
        gradTwist(i, 7) =  1;

        if ( plate->v_bendingElement[i].sign_1 < 0 )
        {
            gradTwist(i, 3) = - gradTwist(i, 3);
            
        }

        if ( plate->v_bendingElement[i].sign_2 < 0 )
        {
            gradTwist(i, 7) = - gradTwist(i, 7);
        }
        
    }

    for (int i = 0; i < plate->bendingNum; i++)
    {
        deltam = plate->v_bendingElement[i].theta_2 - plate->v_bendingElement[i].theta_1;

        value = GJ / plate->v_bendingElement[i].voroniLength * (deltam + plate->v_bendingElement[i].refTwist);
        
        f = - value * gradTwist.row(i);
        localDOF = plate->v_bendingElement[i].arrayNum;

        for (int k = 0; k < 11; k++)
		{
			ind = localDOF(k);
			stepper->addForce(ind, -f[k]); // subtracting elastic force
		}
    }
}

void elasticTwistingForce::computeJt()
{
    for (int i = 0; i < plate->bendingNum; i++)
    {
        norm_e = plate->v_bendingElement[i].norm_1;
        norm_f = plate->v_bendingElement[i].norm_2;
        te = plate->v_bendingElement[i].t_1;
        tf = plate->v_bendingElement[i].t_2;

        norm2_e=norm_e*norm_e;
        norm2_f=norm_f*norm_f;

        kbLocal = plate->v_bendingElement[i].kb;

        chi=1.0+te.dot(tf);
        tilde_t=(te+tf)/chi;

        crossMat(te,teMatrix);

        D2mDe2 = -0.25 / norm2_e * (kbLocal * (te+tilde_t).transpose()
            + (te+tilde_t) * kbLocal.transpose());
        D2mDf2 = -0.25 / norm2_f * (kbLocal * (tf+tilde_t).transpose()
            + (tf+tilde_t) * kbLocal.transpose());
        D2mDeDf = 0.5  / (norm_e*norm_f) * (2.0 / chi * teMatrix
            - kbLocal*tilde_t.transpose());
        D2mDfDe = D2mDeDf.transpose();

        DDtwist.block(0,0,3,3) = D2mDe2;
        DDtwist.block(0,4,3,3) =-D2mDe2 + D2mDeDf;
        DDtwist.block(4,0,3,3) =-D2mDe2 + D2mDfDe;
        DDtwist.block(4,4,3,3) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
        DDtwist.block(0,8,3,3) =-D2mDeDf;
        DDtwist.block(8,0,3,3) =-D2mDfDe;
        DDtwist.block(8,4,3,3) = D2mDfDe - D2mDf2;
        DDtwist.block(4,8,3,3) = D2mDeDf - D2mDf2;
        DDtwist.block(8,8,3,3) = D2mDf2;

        gradTwistLocal = gradTwist.row(i);
        
        milen = -1 / plate->v_bendingElement[i].voroniLength;

        deltam = plate->v_bendingElement[i].theta_2 - plate->v_bendingElement[i].theta_1;

        Jtt = GJ * milen * ((deltam + plate->v_bendingElement[i].refTwist) 
            * DDtwist + gradTwistLocal * gradTwistLocal.transpose());

        localDOF = plate->v_bendingElement[i].arrayNum;

        for (int j=0;j<11;j++)
        {
            for (int k=0;k<11;k++)
            {
				ind1 = localDOF(j);
				ind2 = localDOF(k);
				stepper->addJacobian(ind1, ind2, - Jtt(k,j));
            }
        }
    }
}

// Utility
void elasticTwistingForce::crossMat(const Vector3d &a,Matrix3d &b)
{
	b<<0,-a(2),a(1),
	a(2),0,-a(0),
	-a(1),a(0),0;
}

void elasticTwistingForce::setFirstJacobian()
{
    for (int k = 0; k < plate->bendingNum; k++)
    {
        localDOF = plate->v_bendingElement[k].arrayNum;

        for (int i = 0; i < 11; i++)
        {
            for (int j = 0; j < 11; j++)
            {
                stepper->addJacobian(localDOF(i), localDOF(j), 1);
            }
        }
    }
}