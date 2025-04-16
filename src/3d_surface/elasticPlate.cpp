#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_thickness, double m_Possion, double m_dt)
{
	YoungM = m_YoungM;
	density = m_density;
	thickness = m_thickness;
	Possion = m_Possion;
	dt = m_dt;

	EA = 2 * YoungM * thickness * sqrt(3) / 4;
	EI = 2 * YoungM * thickness * thickness * thickness / ( 12 * sqrt(3) );

	setupGeometry();

	ndof = 3 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	for (int i = 0; i < nv; i++)
	{
		x(3 * i + 0) = v_nodes[i](0);
		x(3 * i + 1) = v_nodes[i](1);
		x(3 * i + 2) = v_nodes[i](2);
	}
	x0 = x;

	// set up edge element
	readInputEdge();

	// set up triangular element
	readInputTriangular();

	// set up bending element
	computeBendingPair();

	// update all pairs
	updateEdgePair();
	updateBendingPair();

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
	boundaryIndex = VectorXi::Zero(nv);

	massArray = VectorXd::Zero(ndof);

	for (int i = 0; i < triangularNum; i++)
	{
		int index1 = v_triangularElement[i].nv_1;
		int index2 = v_triangularElement[i].nv_2;
		int index3 = v_triangularElement[i].nv_3;

		massArray(3 * index1 + 0) = massArray(3 * index1 + 0) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index1 + 1) = massArray(3 * index1 + 1) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index1 + 2) = massArray(3 * index1 + 2) + v_triangularElement[i].area * thickness * density / 3;
		
		massArray(3 * index2 + 0) = massArray(3 * index2 + 0) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index2 + 1) = massArray(3 * index2 + 1) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index2 + 2) = massArray(3 * index2 + 2) + v_triangularElement[i].area * thickness * density / 3;

		massArray(3 * index3 + 0) = massArray(3 * index3 + 0) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index3 + 1) = massArray(3 * index3 + 1) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index3 + 2) = massArray(3 * index3 + 2) + v_triangularElement[i].area * thickness * density / 3;
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

	uCurrent(0) = u(3 * i + 0);
	uCurrent(1) = u(3 * i + 1);
	uCurrent(2) = u(3 * i + 2);

	return uCurrent;
}

void elasticPlate::readInputEdge()
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

		m_edgeElement.refLength = (m_edgeElement.x_1 - m_edgeElement.x_2).norm();
		m_edgeElement.edgeLength = m_edgeElement.refLength;

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
}

void elasticPlate::readInputTriangular()
{
	triangularNum = 0;
	v_triangularElement.clear();

	for (int i = 0; i < triangular.size(); i++)
	{
		Vector3i triangularCurrent = triangular[i];

		triangularElement m_triangularElement;

		m_triangularElement.nv_1 = triangularCurrent(0);
		m_triangularElement.nv_2 = triangularCurrent(1);
		m_triangularElement.nv_3 = triangularCurrent(2);

		m_triangularElement.x_1 = getVertex(m_triangularElement.nv_1);
		m_triangularElement.x_2 = getVertex(m_triangularElement.nv_2);
		m_triangularElement.x_3 = getVertex(m_triangularElement.nv_3);

		m_triangularElement.v_c = (m_triangularElement.x_1 + m_triangularElement.x_2 + m_triangularElement.x_3) / 3;

		Vector3d e_1 = m_triangularElement.x_2 - m_triangularElement.x_1;
		Vector3d e_2 = m_triangularElement.x_3 - m_triangularElement.x_1;

		m_triangularElement.area = 0.5 * ( e_1.cross(e_2) ).norm();

		m_triangularElement.a1 = m_triangularElement.x_3 - m_triangularElement.x_2;
		m_triangularElement.a2 = m_triangularElement.x_1 - m_triangularElement.x_3;
		m_triangularElement.a3 = m_triangularElement.x_2 - m_triangularElement.x_1;

		m_triangularElement.a1 = m_triangularElement.a1 / m_triangularElement.a1.norm();
		m_triangularElement.a2 = m_triangularElement.a2 / m_triangularElement.a2.norm();
		m_triangularElement.a3 = m_triangularElement.a3 / m_triangularElement.a3.norm();

		m_triangularElement.norm = m_triangularElement.a1.cross(m_triangularElement.a2) / ( m_triangularElement.a1.cross(m_triangularElement.a2).norm() );

		m_triangularElement.abar(0, 0) = m_triangularElement.a1(0);
		m_triangularElement.abar(1, 0) = m_triangularElement.a1(1);
		m_triangularElement.abar(2, 0) = m_triangularElement.a1(2);

		m_triangularElement.abar(0, 1) = m_triangularElement.a2(0);
		m_triangularElement.abar(1, 1) = m_triangularElement.a2(1);
		m_triangularElement.abar(2, 1) = m_triangularElement.a2(2);

		m_triangularElement.abar(0, 2) = m_triangularElement.norm(0);
		m_triangularElement.abar(1, 2) = m_triangularElement.norm(1);
		m_triangularElement.abar(2, 2) = m_triangularElement.norm(2);

		m_triangularElement.abarinv = m_triangularElement.abar.inverse();


		v_triangularElement.push_back(m_triangularElement);

		triangularNum = triangularNum + 1;		
	}
}

void elasticPlate::computeBendingPair()
{
	bendingNum = 0;

	v_bendingElement.clear();

	for (int i = 0; i < triangularNum; i++)
	{
		triangularElement m_triangularElement_1 = v_triangularElement[i];

		Vector3i triangularVertex_1;
		triangularVertex_1(0) = m_triangularElement_1.nv_1;
		triangularVertex_1(1) = m_triangularElement_1.nv_2;
		triangularVertex_1(2) = m_triangularElement_1.nv_3;

		for (int j = i+1; j < triangularNum; j++)
		{
			triangularElement m_triangularElement_2 = v_triangularElement[j];

			Vector3i triangularVertex_2;
			triangularVertex_2(0) = m_triangularElement_2.nv_1;
			triangularVertex_2(1) = m_triangularElement_2.nv_2;
			triangularVertex_2(2) = m_triangularElement_2.nv_3;

			int connectNodes = 0;
			Vector2i connectPair;

			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					if (triangularVertex_1(k) == triangularVertex_2(l))
					{
						connectPair(connectNodes) = triangularVertex_1(k);

						connectNodes = connectNodes + 1;
					}
				}
			}

			if (connectNodes == 2)
			{
				bendingElement m_bendingElement;

				m_bendingElement.triangular_1 = i;
				m_bendingElement.triangular_2 = j;

				m_bendingElement.nv_2 = connectPair(0);
				m_bendingElement.nv_3 = connectPair(1);

				for (int k = 0; k < 3; k++)
				{
					if ( triangularVertex_1(k) != connectPair(0) && triangularVertex_1(k) != connectPair(1) )
					{
						m_bendingElement.nv_1 = triangularVertex_1(k);
					}
				}

				for (int k = 0; k < 3; k++)
				{
					if ( triangularVertex_2(k) != connectPair(0) && triangularVertex_2(k) != connectPair(1) )
					{
						m_bendingElement.nv_4 = triangularVertex_2(k);
					}
				}

				m_bendingElement.arrayNum = VectorXi::Zero(12);

				m_bendingElement.arrayNum(0)  = 3 * m_bendingElement.nv_1 + 0;
				m_bendingElement.arrayNum(1)  = 3 * m_bendingElement.nv_1 + 1;
				m_bendingElement.arrayNum(2)  = 3 * m_bendingElement.nv_1 + 2;

				m_bendingElement.arrayNum(3)  = 3 * m_bendingElement.nv_2 + 0;
				m_bendingElement.arrayNum(4)  = 3 * m_bendingElement.nv_2 + 1;
				m_bendingElement.arrayNum(5)  = 3 * m_bendingElement.nv_2 + 2;

				m_bendingElement.arrayNum(6)  = 3 * m_bendingElement.nv_3 + 0;
				m_bendingElement.arrayNum(7)  = 3 * m_bendingElement.nv_3 + 1;
				m_bendingElement.arrayNum(8)  = 3 * m_bendingElement.nv_3 + 2;

				m_bendingElement.arrayNum(9)  = 3 * m_bendingElement.nv_4 + 0;
				m_bendingElement.arrayNum(10) = 3 * m_bendingElement.nv_4 + 1;
				m_bendingElement.arrayNum(11) = 3 * m_bendingElement.nv_4 + 2;

				m_bendingElement.x_1 = getVertex(m_bendingElement.nv_1);
				m_bendingElement.x_2 = getVertex(m_bendingElement.nv_2);
				m_bendingElement.x_3 = getVertex(m_bendingElement.nv_3);
				m_bendingElement.x_4 = getVertex(m_bendingElement.nv_4);

				m_bendingElement.e_1 = m_bendingElement.x_3 - m_bendingElement.x_1;
				m_bendingElement.e_2 = m_bendingElement.x_2 - m_bendingElement.x_1;
				m_bendingElement.e_3 = m_bendingElement.x_2 - m_bendingElement.x_4;
				m_bendingElement.e_4 = m_bendingElement.x_3 - m_bendingElement.x_4;

				m_bendingElement.n_1 = (m_bendingElement.e_1).cross(m_bendingElement.e_2);
				m_bendingElement.n_2 = (m_bendingElement.e_3).cross(m_bendingElement.e_4);

				m_bendingElement.norm_1 = m_bendingElement.n_1 / (m_bendingElement.n_1).norm();
				m_bendingElement.norm_2 = m_bendingElement.n_2 / (m_bendingElement.n_2).norm();

				//m_bendingElement.nBar = (m_bendingElement.norm_1 - m_bendingElement.norm_2).norm();

				m_bendingElement.nBar(0) = 0.0;
				m_bendingElement.nBar(1) = 0.0;
				m_bendingElement.nBar(2) = 0.0;

				v_bendingElement.push_back(m_bendingElement);

				bendingNum = bendingNum + 1;
			}
		}
	}
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
		v_bendingElement[i].x_4 = getVertex(v_bendingElement[i].nv_4);

		v_bendingElement[i].e_1 = v_bendingElement[i].x_3 - v_bendingElement[i].x_1;
		v_bendingElement[i].e_2 = v_bendingElement[i].x_2 - v_bendingElement[i].x_1;
		v_bendingElement[i].e_3 = v_bendingElement[i].x_2 - v_bendingElement[i].x_4;
		v_bendingElement[i].e_4 = v_bendingElement[i].x_3 - v_bendingElement[i].x_4;

		v_bendingElement[i].n_1 = (v_bendingElement[i].e_1).cross(v_bendingElement[i].e_2);
		v_bendingElement[i].n_2 = (v_bendingElement[i].e_3).cross(v_bendingElement[i].e_4);

		v_bendingElement[i].norm_1 = v_bendingElement[i].n_1 / (v_bendingElement[i].n_1).norm();
		v_bendingElement[i].norm_2 = v_bendingElement[i].n_2 / (v_bendingElement[i].n_2).norm();
	}
}

void elasticPlate::prepareForIteration()
{
	updateEdgePair();
	updateBendingPair();
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;
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

void elasticPlate::setupGeometry()
{
	nv = 0;
	edgeNum = 0;
	triangularNum = 0;

	ifstream inFile1;
	inFile1.open("inputdata/3d_surface/inputdata/nodesInput.txt");
	v_nodes.clear();
	double a, b, c;
	while(inFile1 >> a >> b >> c)
	{
		Vector3d xCurrent;

		xCurrent(0) = a;
		xCurrent(1) = b;
		xCurrent(2) = c;

		v_nodes.push_back(xCurrent);
	}
	nv = v_nodes.size();
	inFile1.close();


	ifstream inFile2;
	inFile2.open("inputdata/3d_surface/inputdata/constraint.txt");
	constraint.clear();
	int aa1;
	while(inFile2 >> aa1)
	{
    	constraint.push_back(aa1);
	}
	inFile2.close();


	ifstream inFile3;
	inFile3.open("inputdata/3d_surface/inputdata/edgeInput.txt");
	v_edgeElement.clear();
	int d, e;
	while(inFile3 >> d >> e)
	{
		Vector2i edgeCurrent;

    	edgeCurrent(0) = d;
    	edgeCurrent(1) = e;

    	edge.push_back(edgeCurrent);
	}
	edgeNum = edge.size();
	inFile3.close();


	ifstream inFile4;
	inFile4.open("inputdata/3d_surface/inputdata/triangleInput.txt");
	v_triangularElement.clear();
	int aa, bb, cc;
	while(inFile4 >> aa >> bb >> cc)
	{
		aa = aa;
		bb = bb;
		cc = cc;

		Vector3i triangularCurrent;
		triangularCurrent(0) = aa;
		triangularCurrent(1) = bb;
		triangularCurrent(2) = cc;
		triangular.push_back(triangularCurrent);
	}
	triangularNum = triangular.size();
	inFile4.close();
}