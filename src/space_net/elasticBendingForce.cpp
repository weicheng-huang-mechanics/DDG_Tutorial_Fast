#include "elasticBendingForce.h"

elasticBendingForce::elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;

	Id3<<1,0,0,
         0,1,0,
         0,0,1;

	double EI = plate->EI;

    arrayNum = VectorXi::Zero(9);

    forceVec = VectorXd::Zero(9);
    Jbb = MatrixXd::Zero(9, 9);

    TotalForceVec = VectorXd::Zero(plate->ndof);
}

elasticBendingForce::~elasticBendingForce()
{
	;
}

void elasticBendingForce::computeFb()
{

    TotalForceVec = VectorXd::Zero(plate->ndof);

    for (int i = 0; i < plate->bendingNum; i++)
    {
        arrayNum = plate->v_bendingElement[i].arrayNum;

        x_1 = plate->v_bendingElement[i].x_1;
        x_2 = plate->v_bendingElement[i].x_2;
        x_3 = plate->v_bendingElement[i].x_3;

        nBar = plate->v_bendingElement[i].nBar;

        e_1 = plate->v_bendingElement[i].e_1;
        e_2 = plate->v_bendingElement[i].e_2;

        t_1 = plate->v_bendingElement[i].t_1;
        t_2 = plate->v_bendingElement[i].t_2;

        norm_1 = plate->v_bendingElement[i].norm_1;
        norm_2 = plate->v_bendingElement[i].norm_2;

        voroni = plate->v_bendingElement[i].voroniLength;
        
        gradN1 = ( Id3 - t_1 * t_1.transpose() ) / norm_1;
        gradN2 = ( Id3 - t_2 * t_2.transpose() ) / norm_2;

        // force
        dEde1 = gradN1 * (t_1 - t_2 + nBar);
        dEde2 = gradN2 * (t_2 - t_1 - nBar);

        forceVec = VectorXd::Zero(9);

        forceVec.segment(0,3) = - dEde1;
        forceVec.segment(3,3) =   dEde1 - dEde2;
        forceVec.segment(6,3) =   dEde2;

        // force scale by EI and delta L

        forceVec = - forceVec * plate->EI / voroni;

        for (int j = 0; j < 9; j++)
        {
            ind = arrayNum(j);

            TotalForceVec(ind) = TotalForceVec(ind) + forceVec(j);

            stepper->addForce(ind, - forceVec(j));
        }
    }
}

void elasticBendingForce::computeJb()
{
    for (int i = 0; i < plate->bendingNum; i++)
    {
        arrayNum = plate->v_bendingElement[i].arrayNum;

        x_1 = plate->v_bendingElement[i].x_1;
        x_2 = plate->v_bendingElement[i].x_2;
        x_3 = plate->v_bendingElement[i].x_3;

        nBar = plate->v_bendingElement[i].nBar;

        e_1 = plate->v_bendingElement[i].e_1;
        e_2 = plate->v_bendingElement[i].e_2;

        t_1 = plate->v_bendingElement[i].t_1;
        t_2 = plate->v_bendingElement[i].t_2;

        norm_1 = plate->v_bendingElement[i].norm_1;
        norm_2 = plate->v_bendingElement[i].norm_2;

        voroni = plate->v_bendingElement[i].voroniLength;
        
        gradN1 = ( Id3 - t_1 * t_1.transpose() ) / norm_1;
        gradN2 = ( Id3 - t_2 * t_2.transpose() ) / norm_2;

        // force
        dEde1 = gradN1 * (t_1 - t_2);
        dEde2 = gradN2 * (t_2 - t_1);

        hessionMatrix1 = ( ( gradN1.col(0) * t_1.transpose() + t_1 * gradN1.row(0) + t_1(0) * gradN1 ) ) / ( norm_1 );
        hessionMatrix2 = ( ( gradN1.col(1) * t_1.transpose() + t_1 * gradN1.row(1) + t_1(1) * gradN1 ) ) / ( norm_1 );
        hessionMatrix3 = ( ( gradN1.col(2) * t_1.transpose() + t_1 * gradN1.row(2) + t_1(2) * gradN1 ) ) / ( norm_1 );

        d2Ede12 = gradN1 * gradN1 + hessionMatrix1 * ( t_2(0) - t_1(0) - nBar(0) ) + hessionMatrix2 * ( t_2(1) - t_1(1) - nBar(1) ) + hessionMatrix3 * ( t_2(2) - t_1(2) - nBar(2) );

        hessionMatrix1 = ( ( gradN2.col(0) * t_2.transpose() + t_2 * gradN2.row(0) + t_2(0) * gradN2 ) ) / ( norm_2 );
        hessionMatrix2 = ( ( gradN2.col(1) * t_2.transpose() + t_2 * gradN2.row(1) + t_2(1) * gradN2 ) ) / ( norm_2 );
        hessionMatrix3 = ( ( gradN2.col(2) * t_2.transpose() + t_2 * gradN2.row(2) + t_2(2) * gradN2 ) ) / ( norm_2 );

        d2Ede22 = gradN2 * gradN2 + hessionMatrix1 * ( t_1(0) - t_2(0) + nBar(0) ) + hessionMatrix2 * ( t_1(1) - t_2(1) + nBar(1) ) + hessionMatrix3 * ( t_1(2) - t_2(2) + nBar(2) );


        d2Ede1de2 =  - gradN1 * gradN2; 
        d2Ede2de1 = d2Ede1de2.transpose();

        Jbb = MatrixXd::Zero(9, 9);

        Jbb.block(0,0,3,3) =   d2Ede12;
        Jbb.block(3,3,3,3) =   d2Ede12 - d2Ede1de2 - d2Ede2de1 + d2Ede22;
        Jbb.block(6,6,3,3) =                                     d2Ede22;

        Jbb.block(0,3,3,3) = - d2Ede12 + d2Ede1de2;
        Jbb.block(0,6,3,3) =           - d2Ede1de2;
        Jbb.block(3,6,3,3) =             d2Ede1de2 - d2Ede22;

        Jbb.block(3,0,3,3) = Jbb.block(0,3,3,3).transpose();
        Jbb.block(6,0,3,3) = Jbb.block(0,6,3,3).transpose();
        Jbb.block(6,3,3,3) = Jbb.block(3,6,3,3).transpose();

        Jbb = - Jbb * plate->EI / voroni;

        for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                ind1 = arrayNum(j);
                ind2 = arrayNum(k);

                stepper->addJacobian(ind1, ind2, -Jbb(j,k));
            }
        }
    }
}

void elasticBendingForce::setFirstJacobian()
{
    for (int i = 0; i < plate->bendingNum; i++)
    {
        arrayNum = plate->v_bendingElement[i].arrayNum;

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