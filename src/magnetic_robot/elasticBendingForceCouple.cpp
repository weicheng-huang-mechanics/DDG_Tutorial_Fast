#include "elasticBendingForceCouple.h"

elasticBendingForceCouple::elasticBendingForceCouple(GlobalStepper &m_stepper, double m_EI, 
    int m_rodDOF, int m_plateDOF)
{
	Globalstepper = &m_stepper;

	Id3<<1,0,0,
         0,1,0,
         0,0,1;

	EI = m_EI;

    indexVec = VectorXi::Zero(15);

    forceVec = VectorXd::Zero(15);
    Jbb = MatrixXd::Zero(15, 15);

    rodDOF = m_rodDOF;
    plateDOF = m_plateDOF;
}

elasticBendingForceCouple::~elasticBendingForceCouple()
{
	;
}

void elasticBendingForceCouple::computeFb(int m_plate_nv_1, int m_plate_nv_2, int m_plate_nv_3, 
        Vector3d m_plate_v1, Vector3d m_plate_v2, Vector3d m_plate_v3,
        int m_rodIndex, Vector3d m_rod_v1, Vector3d m_rod_v2, Vector3d m_nBar)
{
    indexVec(0) = 3 * m_plate_nv_1 + 0;
    indexVec(1) = 3 * m_plate_nv_1 + 1;
    indexVec(2) = 3 * m_plate_nv_1 + 2;

    indexVec(3) = 3 * m_plate_nv_2 + 0;
    indexVec(4) = 3 * m_plate_nv_2 + 1;
    indexVec(5) = 3 * m_plate_nv_2 + 2;

    indexVec(6) = 3 * m_plate_nv_3 + 0;
    indexVec(7) = 3 * m_plate_nv_3 + 1;
    indexVec(8) = 3 * m_plate_nv_3 + 2;

    indexVec(9)  = plateDOF + m_rodIndex * rodDOF + 0;
    indexVec(10) = plateDOF + m_rodIndex * rodDOF + 1;
    indexVec(11) = plateDOF + m_rodIndex * rodDOF + 2;

    indexVec(12) = plateDOF + m_rodIndex * rodDOF + 4;
    indexVec(13) = plateDOF + m_rodIndex * rodDOF + 5;
    indexVec(14) = plateDOF + m_rodIndex * rodDOF + 6;

    nBar = m_nBar;

    x_1 = m_plate_v1;
    x_2 = m_plate_v2;
    x_3 = m_plate_v3;

    e_1 = x_3 - x_1;
    e_2 = x_2 - x_1;

    n_1 = e_1.cross(e_2);

    norm_1 = n_1 / n_1.norm();

    gradN1 = ( Id3 - norm_1 * norm_1.transpose() ) / n_1.norm();

    x_4 = m_rod_v1;
    x_5 = m_rod_v2;

    n_2 = x_5 - x_4;

    norm_2 = n_2 / n_2.norm();

    gradN2 = ( Id3 - norm_2 * norm_2.transpose() ) / n_2.norm();

    dEde1 =   (norm_1 - norm_2 - nBar).transpose() * gradN1 * skemMatrix(e_2);
    dEde2 = - (norm_1 - norm_2 - nBar).transpose() * gradN1 * skemMatrix(e_1);

    dEde3 = - (norm_1 - norm_2 - nBar).transpose() * gradN2;

    forceVec.segment(0,3)  = - dEde1 - dEde2;
    forceVec.segment(3,3)  =   dEde2;
    forceVec.segment(6,3)  =   dEde1;
    forceVec.segment(9,3)  = - dEde3;
    forceVec.segment(12,3) =   dEde3;

    //cout << " error " << (norm_1 - norm_2 - nBar).norm() << " " << forceVec.norm() << endl;

    forceVec = - forceVec * EI;

    for (int j = 0; j < 15; j++)
    {
        ind = indexVec(j);

        Globalstepper->addForce(ind, - forceVec(j));
    }
}

void elasticBendingForceCouple::computeJb(int m_plate_nv_1, int m_plate_nv_2, int m_plate_nv_3, 
        Vector3d m_plate_v1, Vector3d m_plate_v2, Vector3d m_plate_v3,
        int m_rodIndex, Vector3d m_rod_v1, Vector3d m_rod_v2, Vector3d m_nBar)
{
    indexVec(0) = 3 * m_plate_nv_1 + 0;
    indexVec(1) = 3 * m_plate_nv_1 + 1;
    indexVec(2) = 3 * m_plate_nv_1 + 2;

    indexVec(3) = 3 * m_plate_nv_2 + 0;
    indexVec(4) = 3 * m_plate_nv_2 + 1;
    indexVec(5) = 3 * m_plate_nv_2 + 2;

    indexVec(6) = 3 * m_plate_nv_3 + 0;
    indexVec(7) = 3 * m_plate_nv_3 + 1;
    indexVec(8) = 3 * m_plate_nv_3 + 2;

    indexVec(9)  = plateDOF + m_rodIndex * rodDOF + 0;
    indexVec(10) = plateDOF + m_rodIndex * rodDOF + 1;
    indexVec(11) = plateDOF + m_rodIndex * rodDOF + 2;

    indexVec(12) = plateDOF + m_rodIndex * rodDOF + 4;
    indexVec(13) = plateDOF + m_rodIndex * rodDOF + 5;
    indexVec(14) = plateDOF + m_rodIndex * rodDOF + 6;

    nBar = m_nBar;

    x_1 = m_plate_v1;
    x_2 = m_plate_v2;
    x_3 = m_plate_v3;

    e_1 = x_3 - x_1;
    e_2 = x_2 - x_1;

    n_1 = e_1.cross(e_2);

    norm_1 = n_1 / n_1.norm();

    gradN1 = ( Id3 - norm_1 * norm_1.transpose() ) / n_1.norm();

    x_4 = m_rod_v1;
    x_5 = m_rod_v2;

    n_2 = x_5 - x_4;

    norm_2 = n_2 / n_2.norm();

    gradN2 = ( Id3 - norm_2 * norm_2.transpose() ) / n_2.norm();

    hession1 = ( ( gradN1.col(0) * norm_1.transpose() + norm_1 * gradN1.row(0) + norm_1(0) * gradN1 ) ) / ( n_1.norm() );
    hession2 = ( ( gradN1.col(1) * norm_1.transpose() + norm_1 * gradN1.row(1) + norm_1(1) * gradN1 ) ) / ( n_1.norm() );
    hession3 = ( ( gradN1.col(2) * norm_1.transpose() + norm_1 * gradN1.row(2) + norm_1(2) * gradN1 ) ) / ( n_1.norm() );

    hession = ( - gradN1 * gradN1 + hession1 * ( norm_1(0) - norm_2(0) - nBar(0) ) + hession2 * ( norm_1(1) - norm_2(1) - nBar(1) ) + hession3 * ( norm_1(2) - norm_2(2) - nBar(2) ));

    d2Ede12 = skemMatrix(e_2) * hession * skemMatrix(e_2);

    d2Ede22 = skemMatrix(e_1) * hession * skemMatrix(e_1);

    d2Ede1de2 = - skemMatrix(e_1) * hession * skemMatrix(e_2);
    d2Ede2de1 = d2Ede1de2.transpose();

    d2Ede1de3 = skemMatrix(e_2) * gradN1 * gradN2;
    d2Ede3de1 = d2Ede1de3.transpose();

    d2Ede2de3 = - skemMatrix(e_1) * gradN1 * gradN2;
    d2Ede3de2 = d2Ede2de3.transpose();

    hession1 = ( ( gradN2.col(0) * norm_2.transpose() + norm_2 * gradN2.row(0) + norm_2(0) * gradN2 ) ) / ( n_2.norm() );
    hession2 = ( ( gradN2.col(1) * norm_2.transpose() + norm_2 * gradN2.row(1) + norm_2(1) * gradN2 ) ) / ( n_2.norm() );
    hession3 = ( ( gradN2.col(2) * norm_2.transpose() + norm_2 * gradN2.row(2) + norm_2(2) * gradN2 ) ) / ( n_2.norm() );

    d2Ede32 = gradN2 * gradN2 + hession1 * ( norm_1(0) - norm_2(0) - nBar(0) ) + hession2 * ( norm_1(1) - norm_2(1) - nBar(1) ) + hession3 * ( norm_1(2) - norm_2(2) - nBar(2) );

    d2Ede32 = d2Ede32;

    Jbb = MatrixXd::Zero(15, 15);

    Jbb.block(0,0,3,3)   = d2Ede12 + d2Ede22 + d2Ede1de2 + d2Ede2de1;
    Jbb.block(3,3,3,3)   = d2Ede22;
    Jbb.block(6,6,3,3)   = d2Ede12;
    Jbb.block(9,9,3,3)   = d2Ede32;
    Jbb.block(12,12,3,3) = d2Ede32;



    Jbb.block(0,3,3,3) =  - d2Ede22 - d2Ede1de2;
    Jbb.block(3,0,3,3) = Jbb.block(0,3,3,3).transpose();

    Jbb.block(0,6,3,3) =  - d2Ede12 - d2Ede2de1;
    Jbb.block(6,0,3,3) = Jbb.block(0,6,3,3).transpose();

    Jbb.block(0,9,3,3) = d2Ede1de3 + d2Ede2de3;
    Jbb.block(9,0,3,3) = Jbb.block(0,9,3,3).transpose();

    Jbb.block(0,12,3,3) = - d2Ede1de3 - d2Ede2de3;
    Jbb.block(12,0,3,3) = Jbb.block(0,12,3,3).transpose();



    Jbb.block(3,6,3,3) = d2Ede2de1;
    Jbb.block(6,3,3,3) = Jbb.block(3,6,3,3).transpose();

    Jbb.block(3,9,3,3) = - d2Ede2de3;
    Jbb.block(9,3,3,3) = Jbb.block(3,9,3,3).transpose();

    Jbb.block(3,12,3,3) = d2Ede2de3;
    Jbb.block(12,3,3,3) = Jbb.block(3,12,3,3).transpose();



    Jbb.block(6,9,3,3) = - d2Ede1de3;
    Jbb.block(9,6,3,3) = Jbb.block(6,9,3,3).transpose();

    Jbb.block(6,12,3,3) = d2Ede1de3;
    Jbb.block(12,6,3,3) = Jbb.block(6,12,3,3).transpose();


    Jbb.block(9,12,3,3) = - d2Ede32;
    Jbb.block(12,9,3,3) = Jbb.block(9,12,3,3).transpose();

    Jbb = - Jbb * EI;

    for (int j = 0; j < 15; j++)
    {
        for (int k = 0; k < 15; k++)
        {
            ind1 = indexVec(j);
            ind2 = indexVec(k);

            Globalstepper->addJacobian(ind1, ind2, -Jbb(j,k));
        }
    }
}

void elasticBendingForceCouple::setFirstJacobian(int m_plate_nv_1, int m_plate_nv_2, int m_plate_nv_3, int m_rodIndex)
{
    indexVec(0) = 3 * m_plate_nv_1 + 0;
    indexVec(1) = 3 * m_plate_nv_1 + 1;
    indexVec(2) = 3 * m_plate_nv_1 + 2;

    indexVec(3) = 3 * m_plate_nv_2 + 0;
    indexVec(4) = 3 * m_plate_nv_2 + 1;
    indexVec(5) = 3 * m_plate_nv_2 + 2;

    indexVec(6) = 3 * m_plate_nv_3 + 0;
    indexVec(7) = 3 * m_plate_nv_3 + 1;
    indexVec(8) = 3 * m_plate_nv_3 + 2;

    indexVec(9)  = plateDOF + m_rodIndex * rodDOF + 0;
    indexVec(10) = plateDOF + m_rodIndex * rodDOF + 1;
    indexVec(11) = plateDOF + m_rodIndex * rodDOF + 2;

    indexVec(12) = plateDOF + m_rodIndex * rodDOF + 4;
    indexVec(13) = plateDOF + m_rodIndex * rodDOF + 5;
    indexVec(14) = plateDOF + m_rodIndex * rodDOF + 6;

    for (int j=0;j<15;j++)
    {
        for (int k=0;k<15;k++)
        {
            ind1 = indexVec(j);
            ind2 = indexVec(k);
            Globalstepper->addJacobian(ind1, ind2, 1);
        }
    }
}

Matrix3d elasticBendingForceCouple::skemMatrix(Vector3d a)
{
    Matrix3d b;

    b<<0,a(2),-a(1),
    -a(2),0,a(0),
    a(1),-a(0),0;

    return b;
}