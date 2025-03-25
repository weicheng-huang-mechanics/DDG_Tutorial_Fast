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

	crossSectionalArea = M_PI * radius * radius;

	setupGeometry();

	ndof = 3 * nv;
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
}

void elasticPlate::setupGeometry()
{
	v_nodes.clear();
    edge.clear();
    
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

void elasticPlate::updateEdgePair()
{
	for (int i = 0; i < edgeNum; i++)
	{
		v_edgeElement[i].x_1 = getVertex(v_edgeElement[i].nv_1);
		v_edgeElement[i].x_2 = getVertex(v_edgeElement[i].nv_2);
	}
}