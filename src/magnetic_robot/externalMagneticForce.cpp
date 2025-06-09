#include "externalMagneticForce.h"

externalMagneticForce::externalMagneticForce(elasticRod &m_rod, timeStepper &m_stepper, 
	Vector3d m_bAVector, Vector3d m_bRVector, double m_muZero, double m_omega)
{
	rod = &m_rod;
	stepper = &m_stepper;

	baVector_ref = m_bAVector;
	brVector_ref = m_bRVector;

	muZero = m_muZero;
	omega = m_omega;

	Id3<<1,0,0,
         0,1,0,
         0,0,1;

    force.setZero(7, 1);
    jacob.setZero(7, 7);
}

externalMagneticForce::~externalMagneticForce()
{
	;
}

void externalMagneticForce::computeFm(double m_currentTime, double m_rol, double m_row, int m_cut)
{
	for (int i = 0; i < rod->ne; i++)
	{
		m1_current = rod->m1_old.row(i);
		m2_current = rod->m2_old.row(i);
		m3_current = rod->tangent_old.row(i);

		m1_start = rod->m1_initial.row(i);
		m2_start = rod->m2_initial.row(i);
		m3_start = rod->tangent_initial.row(i);

		x1 = rod->getVertexOld(i); 
		x2 = rod->getVertexOld(i+1); 

		edge = (x2 - x1).norm();

		brVector.setZero(3, 1);
		
		brVector(1) = brVector_ref(2) * sin(m_rol * 4 * M_PI);
		brVector(2) = brVector_ref(2) * cos(m_rol * 4 * M_PI);
		

		baVector.setZero(3, 1);
		baVector(1) = baVector_ref(1) * cos(omega * m_currentTime);
		baVector(2) = - baVector_ref(1) * sin(omega * m_currentTime);

		gradientBa.setZero(3, 3);

		Mag = (m1_start.dot(brVector) * m1_current + m2_start.dot(brVector) * m2_current + m3_start.dot(brVector) * m3_current);

		dm3de = ( Id3 - m3_current * m3_current.transpose() ) / edge;
		dm1de = - ( m3_current * m1_current.transpose() ) / edge;
		dm2de = - ( m3_current * m2_current.transpose() ) / edge;

		dm1dtheta =  m2_current;
		dm2dtheta = -m1_current;
		dm3dtheta.setZero(3,1);

		dMde = m1_start.dot(brVector) * dm1de + m2_start.dot(brVector) * dm2de + m3_start.dot(brVector) * dm3de;
		dMdtheta = m1_start.dot(brVector) * dm1dtheta + m2_start.dot(brVector) * dm2dtheta + m3_start.dot(brVector) * dm3dtheta;

		dEde = dMde.transpose() * baVector;
		dEdtheta = dMdtheta.dot(baVector);

		force.setZero(7, 1);

		force.segment(0, 3) = - dEde + (gradientBa * Mag) / 2;
		force(3) = dEdtheta;
		force.segment(4, 3) = dEde + (gradientBa * Mag) / 2;

		force = - edge * ( force * rod->crossSectionalArea / muZero);

		for (int k = 0; k < 7; k++)
		{
			int ind = 4 * i + k;
			stepper->addForce(ind, -force[k]); // subtracting elastic force
		}

	}
}

void externalMagneticForce::computeJm()
{
	;
}