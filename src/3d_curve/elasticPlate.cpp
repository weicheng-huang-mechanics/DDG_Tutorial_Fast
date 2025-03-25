#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_radius, 
		double m_Possion, double m_dt)
{
	YoungM = m_YoungM;
	density = m_density;
	radius = m_radius;
	Possion = m_Possion;
	dt = m_dt;

	EA = YoungM * M_PI * radius * radius;
	EI = YoungM * M_PI * radius * radius * radius * radius / 4;
	GJ = YoungM/(2.0*(1.0+Possion)) * M_PI * radius * radius * radius * radius / 2;

	crossSectionalArea = M_PI * radius * radius;

	setupGeometry();

	ndof = 3 * nv + ne;
	x = VectorXd::Zero(ndof);
	x_initial = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	for (int i = 0; i < nv; i++)
	{
		x(3 * i + 0) = v_nodes[i](0);
		x(3 * i + 1) = v_nodes[i](1);
		x(3 * i + 2) = v_nodes[i](2);
	}
	x0 = x;

	x_initial = x;

	computeEdge();
	computeBending();

	setupMass();

	//set up constraint map
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{
	massArray = VectorXd::Zero(ndof);

	double deltaMass;

	int index1;
	int index2;

	for (int i = 0; i < edgeNum; i++)
	{
		deltaMass = M_PI * radius * radius * density * v_edgeElement[i].refLength / 2;

		index1 = v_edgeElement[i].nv_1;
		index2 = v_edgeElement[i].nv_2;

		massArray(3 * index1 + 0) = massArray(3 * index1 + 0) + deltaMass;
		massArray(3 * index1 + 1) = massArray(3 * index1 + 1) + deltaMass;
		massArray(3 * index1 + 2) = massArray(3 * index1 + 2) + deltaMass;
	
		massArray(3 * index2 + 0) = massArray(3 * index2 + 0) + deltaMass;
		massArray(3 * index2 + 1) = massArray(3 * index2 + 1) + deltaMass;
		massArray(3 * index2 + 2) = massArray(3 * index2 + 2) + deltaMass;

		massArray(3 * nv + i) = deltaMass * radius * radius / 2;
	}

}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector3d position, int k)
{
	isConstrained[3 * k + 0] = 1;
	isConstrained[3 * k + 1] = 1;
	isConstrained[3 * k + 2] = 1;
	
	// Store in the constrained dof vector
	x(3 * k + 0) = position(0);
	x(3 * k + 1) = position(1);
	x(3 * k + 2) = position(2);
}


void elasticPlate::setConstraint(double position, int k)
{
	isConstrained[k] = 1;

	// Store in the constrained dof vector
	x(k) = position;
}

void elasticPlate::setThetaBoundaryCondition(double position, int k)
{
	isConstrained[3 * nv + k] = 1;

	x(3 * nv + k) = position;
}

void elasticPlate::setOneVertexBoundaryCondition(double position, int i, int k)
{
	isConstrained[3 * i + k] = 1;
	
	// Store in the constrained dof vector
	x(3 * i + k) = position;
}

Vector3d elasticPlate::getVertex(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x(3 * i + 0);
	xCurrent(1) = x(3 * i + 1);
	xCurrent(2) = x(3 * i + 2);
	
	return xCurrent;
}

double elasticPlate::getTheta(int i)
{
	double thetaCurrent;

	thetaCurrent = x(3 * nv + i);
	
	return thetaCurrent;
}

Vector3d elasticPlate::getVertexStart(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x_initial(3 * i + 0);
	xCurrent(1) = x_initial(3 * i + 1);
	xCurrent(2) = x_initial(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVertexOld(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x0(3 * i + 0);
	xCurrent(1) = x0(3 * i + 1);
	xCurrent(2) = x0(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVelocity(int i)
{
	Vector3d uCurrent;

	uCurrent(0) = ( x(3 * i + 0) - x0(3 * i + 0) ) / dt;
	uCurrent(1) = ( x(3 * i + 1) - x0(3 * i + 1) ) / dt;
	uCurrent(2) = ( x(3 * i + 2) - x0(3 * i + 2) ) / dt;
	
	return uCurrent;
}

void elasticPlate::updateTimeStep()
{
	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;

	for (int i = 0; i < bendingNum; i++)
	{
		v_bendingElement[i].t_1_old = v_bendingElement[i].t_1;
		v_bendingElement[i].d_11_old = v_bendingElement[i].d_11;
		v_bendingElement[i].d_12_old = v_bendingElement[i].d_12;

		v_bendingElement[i].t_2_old = v_bendingElement[i].t_2;
		v_bendingElement[i].d_21_old = v_bendingElement[i].d_21;
		v_bendingElement[i].d_22_old = v_bendingElement[i].d_22;

		v_bendingElement[i].refTwist_old = v_bendingElement[i].refTwist;
	}
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
	}
}

void elasticPlate::prepareForIteration()
{
	updateEdgePair();
	updateBendingPair();
}

void elasticPlate::setupGeometry()
{
	v_nodes.clear();
    edge.clear();
    bending.clear();

    ifstream inFile1;
	inFile1.open("inputdata/nodesInput.txt");
	double a, b, c;
	nv = 0;
	while(inFile1 >> a >> b >> c)
	{
		Vector3d xCurrent;

		xCurrent(0) = a;
		xCurrent(1) = b;
		xCurrent(2) = c;

		nv = nv + 1;

		v_nodes.push_back(xCurrent);
	}
	inFile1.close();

	ifstream inFile2;
	inFile2.open("inputdata/edgeInput.txt");
	int d, e;
	ne = 0;
	while(inFile2 >> d >> e)
	{
		Vector2i edgeCurrent;

    	edgeCurrent(0) = d;
    	edgeCurrent(1) = e;

    	ne = ne + 1;

    	edge.push_back(edgeCurrent);
	}
	inFile2.close();

	ifstream inFile3;
	inFile3.open("inputdata/bendingInput.txt");
	int f, g;
	while(inFile3 >> f >> g)
	{
		Vector2i bendingCurrent;

    	bendingCurrent(0) = f;
    	bendingCurrent(1) = g; 

    	bending.push_back(bendingCurrent);
	}
	inFile3.close();

	ifstream inFile4;
	inFile4.open("inputdata/constraint.txt");
	int aa1;
	while(inFile4 >> aa1)
	{
    	constraint.push_back(aa1);
	}
	inFile4.close();
}

void elasticPlate::computeEdge()
{
	edgeNum = 0;
	v_edgeElement.clear();

	for (int i = 0; i < edge.size(); i++)
	{
		Vector2i edgeCurrent = edge[i];

		edgeElement m_edgeElement;

		m_edgeElement.nv_1 = edgeCurrent(0);
		m_edgeElement.nv_2 = edgeCurrent(1);

		m_edgeElement.x_1 = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2 = getVertex(m_edgeElement.nv_2);

		m_edgeElement.refLength = (m_edgeElement.x_2- m_edgeElement.x_1).norm();
		m_edgeElement.edgeLength = m_edgeElement.refLength;

		// compute frame
		Vector3d t_temp;
		Vector3d d1Tmp;

		m_edgeElement.t1 = (m_edgeElement.x_2 - m_edgeElement.x_1) / m_edgeElement.refLength;
		t_temp(0) = 0.0;
		t_temp(1) = 0.0;
		t_temp(2) = 1.0;
		d1Tmp = m_edgeElement.t1.cross(t_temp);
		if (fabs(d1Tmp.norm()) < 1.0e-6)
		{
			t_temp(0) = 0.0;
			t_temp(1) = 1.0;
			t_temp(2) = 0.0;
			d1Tmp = m_edgeElement.t1.cross(t_temp);
		}
		m_edgeElement.d1 = d1Tmp;
		m_edgeElement.d2 = m_edgeElement.t1.cross(d1Tmp);

		// material frame
		m_edgeElement.theta = 0.0;


		m_edgeElement.arrayNum = VectorXi::Zero(6);

		m_edgeElement.arrayNum(0) = 3 * m_edgeElement.nv_1 + 0;
		m_edgeElement.arrayNum(1) = 3 * m_edgeElement.nv_1 + 1;
		m_edgeElement.arrayNum(2) = 3 * m_edgeElement.nv_1 + 2;
		
		m_edgeElement.arrayNum(3) = 3 * m_edgeElement.nv_2 + 0;
		m_edgeElement.arrayNum(4) = 3 * m_edgeElement.nv_2 + 1;
		m_edgeElement.arrayNum(5) = 3 * m_edgeElement.nv_2 + 2;
		
		v_edgeElement.push_back(m_edgeElement);

		edgeNum = edgeNum + 1;
	}

	//cout << "total edge number " << edgeNum << endl;
}

void elasticPlate::computeBending()
{
	bendingNum = 0;
	v_bendingElement.clear();

	for (int i = 0; i < bending.size(); i++)
	{
		Vector2i bendingCurrent = bending[i];

		bendingElement m_bendingElement;

		int ne_1 = bendingCurrent[0];
		m_bendingElement.ne_1 = ne_1;
		int ne_2 = bendingCurrent[1];
		m_bendingElement.ne_2 = ne_2;

		edgeElement edge_1 = v_edgeElement[ne_1];
		edgeElement edge_2 = v_edgeElement[ne_2];

		int nv_1 = edge_1.nv_1;
		int nv_2 = edge_1.nv_2;
		int nv_3 = edge_2.nv_1;
		int nv_4 = edge_2.nv_2;

		if (nv_1 == nv_3)
		{
			m_bendingElement.nv_1 = nv_2;
			m_bendingElement.nv_2 = nv_1;
			m_bendingElement.nv_3 = nv_4;
		}

		if (nv_1 == nv_4)
		{
			m_bendingElement.nv_1 = nv_2;
			m_bendingElement.nv_2 = nv_1;
			m_bendingElement.nv_3 = nv_3;
		}

		if (nv_2 == nv_3)
		{
			m_bendingElement.nv_1 = nv_1;
			m_bendingElement.nv_2 = nv_2;
			m_bendingElement.nv_3 = nv_4;
		}

		if (nv_2 == nv_4)
		{
			m_bendingElement.nv_1 = nv_1;
			m_bendingElement.nv_2 = nv_2;
			m_bendingElement.nv_3 = nv_3;
		}

		m_bendingElement.x_1 = getVertex(m_bendingElement.nv_1);
		m_bendingElement.x_2 = getVertex(m_bendingElement.nv_2);
		m_bendingElement.x_3 = getVertex(m_bendingElement.nv_3);

		m_bendingElement.e_1 = m_bendingElement.x_2 - m_bendingElement.x_1;
		m_bendingElement.e_2 = m_bendingElement.x_3 - m_bendingElement.x_2;

		m_bendingElement.norm_1 =  m_bendingElement.e_1.norm();
		m_bendingElement.norm_2 =  m_bendingElement.e_2.norm();

		m_bendingElement.voroniLength = (m_bendingElement.norm_1 + m_bendingElement.norm_2) / 2;

		Vector3d t_temp;
		Vector3d d1Tmp;

		// first frame
		m_bendingElement.t_1 = m_bendingElement.e_1 / m_bendingElement.norm_1;
		t_temp(0) = 0.0;
		t_temp(1) = 0.0;
		t_temp(2) = 1.0;
		d1Tmp = m_bendingElement.t_1.cross(t_temp);
		if (fabs(d1Tmp.norm()) < 1.0e-6)
		{
			t_temp(0) = 0.0;
			t_temp(1) = 1.0;
			t_temp(2) = 0.0;
			d1Tmp = m_bendingElement.t_1.cross(t_temp);
		}
		m_bendingElement.d_11 = d1Tmp;
		m_bendingElement.d_12 = m_bendingElement.t_1.cross(d1Tmp);

		// second frame
		m_bendingElement.t_2 = m_bendingElement.e_2 / m_bendingElement.norm_2;
		parallelTansport(m_bendingElement.d_11, m_bendingElement.t_1, m_bendingElement.t_2, m_bendingElement.d_21);
		m_bendingElement.d_22 = m_bendingElement.t_2.cross(m_bendingElement.d_21);


		m_bendingElement.t_1_old = m_bendingElement.t_1;
		m_bendingElement.d_11_old = m_bendingElement.d_11;
		m_bendingElement.d_12_old = m_bendingElement.d_12;

		m_bendingElement.t_2_old = m_bendingElement.t_2;
		m_bendingElement.d_21_old = m_bendingElement.d_21;
		m_bendingElement.d_22_old = m_bendingElement.d_22;

		m_bendingElement.refTwist = 0.0;
		m_bendingElement.refTwist_old = 0.0;

		m_bendingElement.theta_1 = 0.0;
		m_bendingElement.theta_2 = 0.0;

		double cs, ss;

		// material frame
		cs = cos(m_bendingElement.theta_1);
        ss = sin(m_bendingElement.theta_1);
        m_bendingElement.m_11 =  cs * m_bendingElement.d_11 + ss * m_bendingElement.d_12;
        m_bendingElement.m_12 = -ss * m_bendingElement.d_11 + cs * m_bendingElement.d_12;

        cs = cos(m_bendingElement.theta_2);
        ss = sin(m_bendingElement.theta_2);
        m_bendingElement.m_21 =  cs * m_bendingElement.d_21 + ss * m_bendingElement.d_22;
        m_bendingElement.m_22 = -ss * m_bendingElement.d_21 + cs * m_bendingElement.d_22;

        // curvature
		m_bendingElement.kb = 2.0 * m_bendingElement.t_1.cross(m_bendingElement.t_2) / (1.0 + m_bendingElement.t_1.dot(m_bendingElement.t_2));
		m_bendingElement.kappa(0) = 0.5 * m_bendingElement.kb.dot(m_bendingElement.m_12 + m_bendingElement.m_22);
		m_bendingElement.kappa(1) =-0.5 * m_bendingElement.kb.dot(m_bendingElement.m_11 + m_bendingElement.m_21);

		m_bendingElement.kappaBar = m_bendingElement.kappa;

		/*

		cout << " t1: " << m_bendingElement.t_1.transpose() << endl;
		cout << " d11: " << m_bendingElement.d_11.transpose() << endl;
		cout << " d12: " << m_bendingElement.d_12.transpose() << endl;

		cout << " t2: " << m_bendingElement.t_2.transpose() << endl;
		cout << " d21: " << m_bendingElement.d_21.transpose() << endl;
		cout << " d22: " << m_bendingElement.d_22.transpose() << endl;

		cout << m_bendingElement.kappa.transpose() << endl;

		*/

		if ( edge_1.t1.dot(m_bendingElement.t_1) > 0 )
		{
			m_bendingElement.sign_1 = 1;
		}
		else
		{
			m_bendingElement.sign_1 = -1;
		}

		if ( edge_2.t1.dot(m_bendingElement.t_2) > 0 )
		{
			m_bendingElement.sign_2 = 1;
		}
		else
		{
			m_bendingElement.sign_2 = -1;
		}

		//cout << m_bendingElement.sign_1 << " " << m_bendingElement.sign_2 << endl;

		// dof
		m_bendingElement.arrayNum = VectorXi::Zero(11);

		m_bendingElement.arrayNum(0) = 3 * m_bendingElement.nv_1 + 0;
		m_bendingElement.arrayNum(1) = 3 * m_bendingElement.nv_1 + 1;
		m_bendingElement.arrayNum(2) = 3 * m_bendingElement.nv_1 + 2;

		m_bendingElement.arrayNum(3) = 3 * nv + ne_1;
		
		m_bendingElement.arrayNum(4) = 3 * m_bendingElement.nv_2 + 0;
		m_bendingElement.arrayNum(5) = 3 * m_bendingElement.nv_2 + 1;
		m_bendingElement.arrayNum(6) = 3 * m_bendingElement.nv_2 + 2;

		m_bendingElement.arrayNum(7) = 3 * nv + ne_2;
		
		m_bendingElement.arrayNum(8) = 3 * m_bendingElement.nv_3 + 0;
		m_bendingElement.arrayNum(9) = 3 * m_bendingElement.nv_3 + 1;
		m_bendingElement.arrayNum(10) = 3 * m_bendingElement.nv_3 + 2;

		v_bendingElement.push_back(m_bendingElement);

		bendingNum = bendingNum + 1;
	}

	//cout << "total bend number " << bendingNum << endl;
}

void elasticPlate::updateEdgePair()
{
	for (int i = 0; i < edgeNum; i++)
	{
		v_edgeElement[i].x_1 = getVertex(v_edgeElement[i].nv_1);
		v_edgeElement[i].x_2 = getVertex(v_edgeElement[i].nv_2);
		v_edgeElement[i].edgeLength = (v_edgeElement[i].x_1 - v_edgeElement[i].x_2).norm();
	}
}

void elasticPlate::updateBendingPair()
{
	for (int i = 0; i < bendingNum; i++)
	{
		v_bendingElement[i].x_1 = getVertex(v_bendingElement[i].nv_1);
		v_bendingElement[i].x_2 = getVertex(v_bendingElement[i].nv_2);
		v_bendingElement[i].x_3 = getVertex(v_bendingElement[i].nv_3);

		v_bendingElement[i].theta_1 = getTheta(v_bendingElement[i].ne_1);
		v_bendingElement[i].theta_2 = getTheta(v_bendingElement[i].ne_2);

		if ( v_bendingElement[i].sign_1 < 0 )
        {
            v_bendingElement[i].theta_1 = - v_bendingElement[i].theta_1;
        }

        if ( v_bendingElement[i].sign_2 < 0 )
        {
            v_bendingElement[i].theta_2 = - v_bendingElement[i].theta_2;
        }

		v_bendingElement[i].e_1 = v_bendingElement[i].x_2 - v_bendingElement[i].x_1;
		v_bendingElement[i].e_2 = v_bendingElement[i].x_3 - v_bendingElement[i].x_2;

		v_bendingElement[i].norm_1 =  v_bendingElement[i].e_1.norm();
		v_bendingElement[i].norm_2 =  v_bendingElement[i].e_2.norm();

		v_bendingElement[i].t_1 = v_bendingElement[i].e_1 / v_bendingElement[i].norm_1;
		v_bendingElement[i].t_2 = v_bendingElement[i].e_2 / v_bendingElement[i].norm_2;

		// computeTimeParallel
		Vector3d d1, d2, t1;
		Vector3d d1_old, d2_old, t1_old;

		t1_old = v_bendingElement[i].t_1_old;
		d1_old = v_bendingElement[i].d_11_old;
		d2_old = v_bendingElement[i].d_12_old;
		t1 = v_bendingElement[i].t_1;
		parallelTansport(d1_old, t1_old, t1, d1);
		d2 = t1.cross(d1);
		v_bendingElement[i].d_11 = d1;
		v_bendingElement[i].d_12 = d2;

		t1_old = v_bendingElement[i].t_2_old;
		d1_old = v_bendingElement[i].d_21_old;
		d2_old = v_bendingElement[i].d_22_old;
		t1 = v_bendingElement[i].t_2;
		parallelTansport(d1_old, t1_old, t1, d1);
		d2 = t1.cross(d1);
		v_bendingElement[i].d_21 = d1;
		v_bendingElement[i].d_22 = d2;

		// refTwist
		Vector3d u_e, u_f, t_e, t_f, u_temp;

		u_e = v_bendingElement[i].d_11;
		u_f = v_bendingElement[i].d_21;
		t_e = v_bendingElement[i].t_1;
		t_f = v_bendingElement[i].t_2;
		parallelTansport(u_e, t_e, t_f, u_temp);
		rotateAxisAngle(u_temp, t_f, v_bendingElement[i].refTwist_old);
		double sgnAngle = signedAngle(u_temp, u_f, t_f);
        v_bendingElement[i].refTwist = v_bendingElement[i].refTwist_old + sgnAngle;

        // material director
        double cs, ss;

        cs = cos(v_bendingElement[i].theta_1);
        ss = sin(v_bendingElement[i].theta_1);
        v_bendingElement[i].m_11 =  cs * v_bendingElement[i].d_11 + ss * v_bendingElement[i].d_12;
        v_bendingElement[i].m_12 = -ss * v_bendingElement[i].d_11 + cs * v_bendingElement[i].d_12;

        cs = cos(v_bendingElement[i].theta_2);
        ss = sin(v_bendingElement[i].theta_2);
        v_bendingElement[i].m_21 =  cs * v_bendingElement[i].d_21 + ss * v_bendingElement[i].d_22;
        v_bendingElement[i].m_22 = -ss * v_bendingElement[i].d_21 + cs * v_bendingElement[i].d_22;

        // curvature
		v_bendingElement[i].kb = 2.0 * v_bendingElement[i].t_1.cross(v_bendingElement[i].t_2) / (1.0 + v_bendingElement[i].t_1.dot(v_bendingElement[i].t_2));
		v_bendingElement[i].kappa(0) = 0.5 * v_bendingElement[i].kb.dot(v_bendingElement[i].m_12 + v_bendingElement[i].m_22);
		v_bendingElement[i].kappa(1) =-0.5 * v_bendingElement[i].kb.dot(v_bendingElement[i].m_11 + v_bendingElement[i].m_21);
	}
}

void elasticPlate::parallelTansport(const Vector3d &d1_1,const Vector3d &t1, const Vector3d &t2, Vector3d &d1_2)
{
	Vector3d b;
	Vector3d n1,n2;

	b=t1.cross(t2);
	
	if(b.norm()==0)
		d1_2=d1_1;
	else
	{
		b = b / b.norm();
		b = b - b.dot(t1) * t1;
		b = b / b.norm();
		b = b - b.dot(t1) * t2;
		b = b / b.norm();
		
		n1=t1.cross(b);
		n2=t2.cross(b);
		d1_2=d1_1.dot(t1)*t2+d1_1.dot(n1)*n2+d1_1.dot(b)*b;
		d1_2=d1_2-d1_2.dot(t2)*t2;
		d1_2=d1_2/d1_2.norm();
	}
}

double elasticPlate::signedAngle(const Vector3d &u, const Vector3d &v, const Vector3d &n)
{
	//Compute the angle between two vectors
	Vector3d w=u.cross(v);

	double angle=atan2(w.norm(),u.dot(v));
	if (n.dot(w)<0)
		return -angle;
	else 
		return angle;
}

void elasticPlate::rotateAxisAngle(Vector3d &v,const Vector3d &z,const double &theta)
{
	//Compute the vector when it rotates along another vector into certain angle
	if (theta!=0) // if theta=0, v = v
	{
		double cs,ss;
		cs=cos(theta);
		ss=sin(theta);
		v=cs*v+ss*z.cross(v)+z.dot(v)*(1.0-cs)*z;
	}
}