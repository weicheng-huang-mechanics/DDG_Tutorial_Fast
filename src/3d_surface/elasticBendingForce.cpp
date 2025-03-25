#include "elasticBendingForce.h"

elasticBendingForce::elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;

	Id3<<1,0,0,
         0,1,0,
         0,0,1;

	double EI = plate->EI;

    arrayNum = VectorXi::Zero(12);

    forceVec = VectorXd::Zero(12);
    Jbb = MatrixXd::Zero(12, 12);

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
        x_4 = plate->v_bendingElement[i].x_4;

        nBar = plate->v_bendingElement[i].nBar;

        e_1 = x_3 - x_1;
        e_2 = x_2 - x_1;
        e_3 = x_2 - x_4;
        e_4 = x_3 - x_4;

        n_1 = e_1.cross(e_2);
        n_2 = e_3.cross(e_4);

        norm_1 = n_1 / n_1.norm();
        norm_2 = n_2 / n_2.norm();
        
        gradN1 = ( Id3 - norm_1 * norm_1.transpose() ) / n_1.norm();
        gradN2 = ( Id3 - norm_2 * norm_2.transpose() ) / n_2.norm();

        // force
        dEde1 =   1 * (norm_1 - norm_2 - nBar).transpose() * gradN1 * skemMatrix(e_2);
        dEde2 = - 1 * (norm_1 - norm_2 - nBar).transpose() * gradN1 * skemMatrix(e_1);
        dEde3 = - 1 * (norm_1 - norm_2 - nBar).transpose() * gradN2 * skemMatrix(e_4);
        dEde4 =   1 * (norm_1 - norm_2 - nBar).transpose() * gradN2 * skemMatrix(e_3);

        forceVec = VectorXd::Zero(12);

        forceVec.segment(0,3) = - dEde1 - dEde2;
        forceVec.segment(3,3) =   dEde3 + dEde2;
        forceVec.segment(6,3) =   dEde1 + dEde4;
        forceVec.segment(9,3) = - dEde4 - dEde3;

        // force scale by EI and delta L

        forceVec = - forceVec * plate->EI;

        for (int j = 0; j < 12; j++)
        {
            ind = arrayNum(j);

            TotalForceVec(ind) = TotalForceVec(ind) + forceVec(j);

            stepper->addForce(ind, -forceVec(j));
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
        x_4 = plate->v_bendingElement[i].x_4;

        nBar = plate->v_bendingElement[i].nBar;

        e_1 = x_3 - x_1;
        e_2 = x_2 - x_1;
        e_3 = x_2 - x_4;
        e_4 = x_3 - x_4;

        n_1 = e_1.cross(e_2);
        n_2 = e_3.cross(e_4);

        norm_1 = n_1 / n_1.norm();
        norm_2 = n_2 / n_2.norm();

        gradN1 = ( Id3 - norm_1 * norm_1.transpose() ) / n_1.norm();
        gradN2 = ( Id3 - norm_2 * norm_2.transpose() ) / n_2.norm();

        // jacobian
        d2Ede12 = - skemMatrix(e_2) * gradN1 * gradN1 * skemMatrix(e_2) 
        - skemMatrix(e_2) * ( (norm_1(0) - norm_2(0) - nBar(0) ) * hessionMatrix_1(n_1) + ( norm_1(1) - norm_2(1) - nBar(1) ) * hessionMatrix_2(n_1) + ( norm_1(2) - norm_2(2) - nBar(2) ) * hessionMatrix_3(n_1) ) * skemMatrix(e_2);

        d2Ede22 = - skemMatrix(e_1) * gradN1 * gradN1 * skemMatrix(e_1) 
        - skemMatrix(e_1) * ( (norm_1(0) - norm_2(0) - nBar(0) ) * hessionMatrix_1(n_1) + ( norm_1(1) - norm_2(1) - nBar(1) ) * hessionMatrix_2(n_1) + ( norm_1(2) - norm_2(2) - nBar(2) ) * hessionMatrix_3(n_1) ) * skemMatrix(e_1);

        d2Ede32 = - skemMatrix(e_4) * gradN2 * gradN2 * skemMatrix(e_4) 
        + skemMatrix(e_4) * ( (norm_1(0) - norm_2(0) - nBar(0) ) * hessionMatrix_1(n_2) + ( norm_1(1) - norm_2(1) - nBar(1) ) * hessionMatrix_2(n_2) + ( norm_1(2) - norm_2(2) - nBar(2) ) * hessionMatrix_3(n_2) ) * skemMatrix(e_4);

        d2Ede42 = - skemMatrix(e_3) * gradN2 * gradN2 * skemMatrix(e_3) 
        + skemMatrix(e_3) * ( (norm_1(0) - norm_2(0) - nBar(0) ) * hessionMatrix_1(n_2) + ( norm_1(1) - norm_2(1) - nBar(1) ) * hessionMatrix_2(n_2) + ( norm_1(2) - norm_2(2) - nBar(2) ) * hessionMatrix_3(n_2) ) * skemMatrix(e_3);

        d2Ede1de2 = skemMatrix(e_1) * gradN1 * gradN1 * skemMatrix(e_2) 
        + skemMatrix(e_1) * ( (norm_1(0) - norm_2(0) - nBar(0) ) * hessionMatrix_1(n_1) + ( norm_1(1) - norm_2(1) - nBar(1) ) * hessionMatrix_2(n_1) + ( norm_1(2) - norm_2(2) - nBar(2) ) * hessionMatrix_3(n_1) ) * skemMatrix(e_2);
        d2Ede2de1 = d2Ede1de2.transpose();

        d2Ede1de3 = skemMatrix(e_2) * gradN1 * gradN2 * skemMatrix(e_4);
        d2Ede3de1 = d2Ede1de3.transpose();

        d2Ede1de4 = - skemMatrix(e_2) * gradN1 * gradN2 * skemMatrix(e_3);
        d2Ede4de1 = d2Ede1de4.transpose();

        d2Ede2de3 = - skemMatrix(e_1) * gradN1 * gradN2 * skemMatrix(e_4);
        d2Ede3de2 = d2Ede2de3.transpose();

        d2Ede2de4 = skemMatrix(e_1) * gradN1 * gradN2 * skemMatrix(e_3);
        d2Ede4de2 = d2Ede2de4.transpose();

        d2Ede3de4 = skemMatrix(e_3) * gradN2 * gradN2 * skemMatrix(e_4) 
        - skemMatrix(e_3) * ( ( norm_1(0) - norm_2(0) - nBar(0) ) * hessionMatrix_1(n_2) + ( norm_1(1) - norm_2(1) - nBar(1) ) * hessionMatrix_2(n_2) + ( norm_1(2) - norm_2(2) - nBar(2) ) * hessionMatrix_3(n_2) ) * skemMatrix(e_4);
        d2Ede4de3 = d2Ede3de4.transpose();

        Jbb = MatrixXd::Zero(12, 12);

        Jbb.block(0,0,3,3) = d2Ede12 + d2Ede22 + d2Ede1de2 + d2Ede2de1;
        Jbb.block(3,3,3,3) = d2Ede22 + d2Ede32 + d2Ede2de3 + d2Ede3de2;
        Jbb.block(6,6,3,3) = d2Ede12 + d2Ede42 + d2Ede1de4 + d2Ede4de1;
        Jbb.block(9,9,3,3) = d2Ede32 + d2Ede42 + d2Ede3de4 + d2Ede4de3;

        Jbb.block(0,3,3,3) = - d2Ede22 - d2Ede1de2 - d2Ede2de3 - d2Ede1de3;
        Jbb.block(3,0,3,3) = Jbb.block(0,3,3,3).transpose();

        Jbb.block(0,6,3,3) = - d2Ede12 - d2Ede1de4 - d2Ede2de1 - d2Ede2de4;
        Jbb.block(6,0,3,3) = Jbb.block(0,6,3,3).transpose();

        Jbb.block(0,9,3,3) = d2Ede1de3 + d2Ede1de4 + d2Ede2de3 + d2Ede2de4;
        Jbb.block(9,0,3,3) = Jbb.block(0,9,3,3).transpose();

        Jbb.block(3,6,3,3) = d2Ede2de1 + d2Ede2de4 + d2Ede3de1 + d2Ede3de4;
        Jbb.block(6,3,3,3) = Jbb.block(3,6,3,3).transpose();

        Jbb.block(3,9,3,3) = - d2Ede32 - d2Ede3de4 - d2Ede2de4 - d2Ede2de3;
        Jbb.block(9,3,3,3) = Jbb.block(3,9,3,3).transpose();

        Jbb.block(6,9,3,3) = - d2Ede42 - d2Ede4de3 - d2Ede1de4 - d2Ede1de3;
        Jbb.block(9,6,3,3) = Jbb.block(6,9,3,3).transpose();

	
        Jbb = - Jbb * plate->EI;

        for (int j = 0; j < 12; j++)
        {
            for (int k = 0; k < 12; k++)
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

        for (int j = 0; j < 12; j++)
        {
            for (int k = 0; k < 12; k++)
            {
                ind1 = arrayNum(j);
                ind2 = arrayNum(k);

                stepper->addJacobian(ind1, ind2, 1);
            }
        }
    }
}

Matrix3d elasticBendingForce::skemMatrix(Vector3d a)
{
    Matrix3d b;

    b<<0,a(2),-a(1),
    -a(2),0,a(0),
    a(1),-a(0),0;

    return b;
}

Matrix3d elasticBendingForce::hessionMatrix_1(Vector3d a)
{
    double normA = a.norm();

    Vector3d t = a / normA;

    Matrix3d grad = ( Id3 - t * t.transpose() ) / normA;

    Matrix3d matrix_temp = Id3 - t * t.transpose();

    Matrix3d hession_1;

    hession_1 = - ( ( (grad.col(0) * t.transpose() + t * grad.row(0)) + t(0) * grad ) * normA + t(0) * matrix_temp +  ( t * matrix_temp.row(0) + matrix_temp.col(0) * t.transpose() ) ) / (2 * normA * normA);

    return hession_1;
}

Matrix3d elasticBendingForce::hessionMatrix_2(Vector3d a)
{
    double normA = a.norm();

    Vector3d t = a / normA;

    Matrix3d grad = ( Id3 - t * t.transpose() ) / normA;

    Matrix3d matrix_temp = Id3 - t * t.transpose();

    Matrix3d hession_2;

    hession_2 = - ( ( (grad.col(1) * t.transpose() + t * grad.row(1)) + t(1) * grad ) * normA + t(1) * matrix_temp +  ( t * matrix_temp.row(1) + matrix_temp.col(1) * t.transpose() ) ) / (2 * normA * normA);

    return hession_2;
}

Matrix3d elasticBendingForce::hessionMatrix_3(Vector3d a)
{
    double normA = a.norm();

    Vector3d t = a / normA;

    Matrix3d grad = ( Id3 - t * t.transpose() ) / normA;

    Matrix3d matrix_temp = Id3 - t * t.transpose();

    Matrix3d hession_3;

    hession_3 = - ( ( (grad.col(2) * t.transpose() + t * grad.row(2)) + t(2) * grad ) * normA + t(2) * matrix_temp +  ( t * matrix_temp.row(2) + matrix_temp.col(2) * t.transpose() ) ) / (2 * normA * normA);

    return hession_3;
}
