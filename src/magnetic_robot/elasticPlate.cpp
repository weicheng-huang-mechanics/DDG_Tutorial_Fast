#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_thickness, double m_Possion, double m_dt,
	double m_length, double m_width, double m_deltaEdge, double m_startHeight)
{
	YoungM = m_YoungM;
	density = m_density;
	thickness = m_thickness;
	Possion = m_Possion;
	dt = m_dt;

	startHeight = m_startHeight;

	length = m_length;
    width = m_width;
    deltaEdge = m_deltaEdge;
    height = deltaEdge / 2 * sqrt(3);

	EA = 2 * YoungM * thickness * sqrt(3) / 4;
	EI = 2 * YoungM * thickness * thickness * thickness / ( 12 * sqrt(3) );

	setupGeometry();

	// set up nodes
	readInputNodes();

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

	deltaArea = deltaEdge * deltaEdge * sqrt(3) / 4;
	deltaLength = deltaEdge;

	for (int i = 0; i < nv; i++)
	{
		if (boundaryIndex(i) == 0)
		{
			massArray(3 * i + 0) = 2 * deltaArea * thickness * density;
			massArray(3 * i + 1) = 2 * deltaArea * thickness * density;
			massArray(3 * i + 2) = 2 * deltaArea * thickness * density;
		}
		else
		{
			massArray(3 * i + 0) = deltaArea * thickness * density;
			massArray(3 * i + 1) = deltaArea * thickness * density;
			massArray(3 * i + 2) = deltaArea * thickness * density;
		}
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

void elasticPlate::setVertexOneBoundaryCondition(double position, int i, int k)
{
	isConstrained[3 * i + k] = 1;
	
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

	return uCurrent;
}


Vector3d elasticPlate::getVelocityOld(int i)
{
	Vector3d uCurrent;

	uCurrent(0) = u(3 * i + 0);
	uCurrent(1) = u(3 * i + 1);
	uCurrent(2) = u(3 * i + 2);

	return uCurrent;
}

void elasticPlate::readInputNodes()
{
	nv = 0;
	v_nodes.clear();

	for (int i = 0; i < nodes.size(); i++)
	{
		Vector3d xCurrent = nodes[i];

		v_nodes.push_back(xCurrent);

		nv = nv + 1;
	}
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

		Vector3d e_1 = m_triangularElement.x_2 - m_triangularElement.x_1;
		Vector3d e_2 = m_triangularElement.x_3 - m_triangularElement.x_1;

		m_triangularElement.area = 0.5 * ( e_1.cross(e_2) ).norm();

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

				m_bendingElement.nBar = m_bendingElement.norm_1 - m_bendingElement.norm_2;

			//	m_bendingElement.nBar(0) = 0.0;
			//	m_bendingElement.nBar(1) = 0.0;
			//	m_bendingElement.nBar(2) = 0.0;

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
	xNodeNum = length / deltaEdge;
    yNodeNum = width / height;

    if (xNodeNum % 2 == 1)
    {
    	xNodeNum = xNodeNum + 2;
    }

    if (xNodeNum % 2 == 0)
    {
    	xNodeNum = xNodeNum + 1;
    } 

    if (yNodeNum % 2 == 1)
    {
    	yNodeNum = yNodeNum + 2;
    }

    if (yNodeNum % 2 == 0)
    {
    	yNodeNum = yNodeNum + 1;
    }

    ratio = (yNodeNum * height) / (xNodeNum * deltaEdge);

    nodes.clear();
    startnode.clear();
    edge.clear();
    triangular.clear();

    boundaryArray_1.clear();
    boundaryArray_2.clear();

    // set up nodes
    temp = 0;

    for (int i = 0; i < yNodeNum; i++)
    {
        // odd line
        if ( i % 2 == 0)
        {
            for (int j = 0; j < xNodeNum; j++)
            {
                if (j == 0)
                {
                    startnode.push_back(temp);
                }

                Vector3d xCurrent;

                xCurrent(0) = j * deltaEdge;
                xCurrent(1) = i * height;
                xCurrent(2) = startHeight;

                nodes.push_back(xCurrent);
                temp = temp + 1;
            }
        }

        // even line
        if ( i % 2 == 1)
        {
            for (int j = 0; j < xNodeNum-1; j++)
            {
                if (j == 0)
                {
                    startnode.push_back(temp);
                }

                Vector3d xCurrent;

                xCurrent(0) = j * deltaEdge + deltaEdge / 2;
                xCurrent(1) = i * height;
                xCurrent(2) = startHeight;

                nodes.push_back(xCurrent);
                temp = temp + 1;
            }
        }
    }

    // set up edge
    for (int i = 0; i < yNodeNum; i++)
    {
        // odd line
        if ( i % 2 == 0)
        {
            for (int j = 0; j < xNodeNum-1; j++)
            {
                Vector2i edgeCurrent;

                edgeCurrent(0) = startnode[i] + j;
                edgeCurrent(1) = startnode[i] + j + 1;

                edge.push_back(edgeCurrent);
            }
        }

        if ( i % 2 == 1)
        {
            for (int j = 0; j < xNodeNum-2; j++)
            {
                Vector2i edgeCurrent;

                edgeCurrent(0) = startnode[i] + j;
                edgeCurrent(1) = startnode[i] + j + 1;

                edge.push_back(edgeCurrent);
            }
        }
    }

    for (int i = 0; i < yNodeNum-1; i++)
    {
        if (i % 2 == 0)
        {
            for (int j = 0; j < xNodeNum; j++)
            {
                if (j == 0)
                {
                    Vector2i edgeCurrent;
                    edgeCurrent(0) = startnode[i];
                    edgeCurrent(1) = startnode[i+1];

                    edge.push_back(edgeCurrent);
                }

                if (j > 0 && j < xNodeNum - 1)
                {
                    Vector2i edgeCurrent_1;
                    edgeCurrent_1(0) = startnode[i] + j;
                    edgeCurrent_1(1) = startnode[i+1] + j - 1;

                    edge.push_back(edgeCurrent_1);

                    Vector2i edgeCurrent_2;
                    edgeCurrent_2(0) = startnode[i] + j;
                    edgeCurrent_2(1) = startnode[i+1] + j;

                    edge.push_back(edgeCurrent_2);
                }

                if (j == xNodeNum - 1)
                {
                    Vector2i edgeCurrent;
                    edgeCurrent(0) = startnode[i] + j;
                    edgeCurrent(1) = startnode[i+1] + j - 1;

                    edge.push_back(edgeCurrent);
                }
            }
        }

        if ( i % 2 == 1)
        {
            for (int j = 0; j < xNodeNum - 1; j++)
            {
                Vector2i edgeCurrent_1;
                edgeCurrent_1(0) = startnode[i] + j;
                edgeCurrent_1(1) = startnode[i+1] + j;

                edge.push_back(edgeCurrent_1);

                Vector2i edgeCurrent_2;
                edgeCurrent_2(0) = startnode[i] + j;
                edgeCurrent_2(1) = startnode[i+1] + j + 1;

                edge.push_back(edgeCurrent_2);
            }
        }
    }

    /**********/
    
    edgeOutput.clear();
    edgeOutput = edge;

    for (int i = 0; i < yNodeNum-1; i++)
    {
        if (i % 2 == 0)
        {
        	Vector2i edgeCurrent;
        	edgeCurrent(0) = startnode[i];
        	edgeCurrent(1) = startnode[i+2];

        	edgeOutput.push_back(edgeCurrent);
        }
    }

    for (int i = 0; i < yNodeNum-2; i++)
    {
        if (i % 2 == 1)
        {
        	Vector2i edgeCurrent;
        	edgeCurrent(0) = startnode[i]-1;
        	edgeCurrent(1) = startnode[i+2]-1;

        	edgeOutput.push_back(edgeCurrent);
        }
    }

    Vector2i edgeCurrent;
    edgeCurrent(0) = startnode[startnode.size()-2]-1;
    edgeCurrent(1) = nv-1;

    edgeOutput.push_back(edgeCurrent);

    /**********/

    // set up triangular
    for (int i = 0; i < yNodeNum; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < xNodeNum-1; j++)
            {
                Vector3i triangularCurrent;

                triangularCurrent(0) = startnode[i] + j;
                triangularCurrent(1) = startnode[i] + j + 1;
                triangularCurrent(2) = startnode[i+1] + j;

                triangular.push_back(triangularCurrent);
            }
        }

        if (i == yNodeNum-1)
        {
            if (i % 2 == 0)
            {
                for (int j = 0; j < xNodeNum-1; j++)
                {
                    Vector3i triangularCurrent;

                    triangularCurrent(0) = startnode[i] + j;
                    triangularCurrent(1) = startnode[i] + j + 1;
                    triangularCurrent(2) = startnode[i-1] + j;

                    triangular.push_back(triangularCurrent);
                }
            }

            if (i % 2 == 1)
            {
                for (int j = 0; j < xNodeNum-2; j++)
                {
                    Vector3i triangularCurrent;

                    triangularCurrent(0) = startnode[i] + j;
                    triangularCurrent(1) = startnode[i] + j + 1;
                    triangularCurrent(2) = startnode[i-1] + j + 1;

                    triangular.push_back(triangularCurrent);
                }
            }
        }

        if (i > 0 && i < yNodeNum-1)
        {
            if (i % 2 == 0)
            {
                for (int j = 0; j < xNodeNum-1; j++)
                {
                    Vector3i triangularCurrent_1;
                    triangularCurrent_1(0) = startnode[i] + j;
                    triangularCurrent_1(1) = startnode[i] + j + 1;
                    triangularCurrent_1(2) = startnode[i-1] + j;
                    triangular.push_back(triangularCurrent_1);

                    Vector3i triangularCurrent_2;
                    triangularCurrent_2(0) = startnode[i] + j;
                    triangularCurrent_2(1) = startnode[i] + j + 1;
                    triangularCurrent_2(2) = startnode[i+1] + j;
                    triangular.push_back(triangularCurrent_2);
                }
            }

            if (i % 2 == 1)
            {
                for (int j = 0; j < xNodeNum-2; j++)
                {
                    Vector3i triangularCurrent_1;
                    triangularCurrent_1(0) = startnode[i] + j;
                    triangularCurrent_1(1) = startnode[i] + j + 1;
                    triangularCurrent_1(2) = startnode[i-1] + j + 1;
                    triangular.push_back(triangularCurrent_1);

                    Vector3i triangularCurrent_2;
                    triangularCurrent_2(0) = startnode[i] + j;
                    triangularCurrent_2(1) = startnode[i] + j + 1;
                    triangularCurrent_2(2) = startnode[i+1] + j + 1;
                    triangular.push_back(triangularCurrent_2);
                }
            }
        }
    }


    /**********/

    triangularOutput.clear();
    triangularOutput = triangular;

    for (int i = 0; i < yNodeNum-1; i++)
    {
    	if (i % 2 == 0)
    	{
    		Vector3i triangularCurrent;
    		triangularCurrent(0) = startnode[i];
    		triangularCurrent(1) = startnode[i+1];
    		triangularCurrent(2) = startnode[i+2];
    		triangularOutput.push_back(triangularCurrent);
    	}
    }

    for (int i = 0; i < yNodeNum-3; i++)
    {
    	if (i % 2 == 0)
    	{
    		Vector3i triangularCurrent;
    		triangularCurrent(0) = startnode[i+1]-1;
    		triangularCurrent(1) = startnode[i+2]-1;
    		triangularCurrent(2) = startnode[i+3]-1;
    		triangularOutput.push_back(triangularCurrent);
    	}
    }

    Vector3i triangularCurrent;
    triangularCurrent(0) = startnode[startnode.size()-2]-1;
    triangularCurrent(1) = startnode[startnode.size()-1]-1;
    triangularCurrent(2) = nv-1;
    triangularOutput.push_back(triangularCurrent);

    /**********/

    // set up first b.c array
    for (int i = 0; i < startnode[2]; i++)
    {
        boundaryArray_1.push_back(i);
    }

    // set up second b.c array
    for (int i = startnode[startnode.size() - 2]; i < nodes.size(); i++)
    {
        boundaryArray_2.push_back(i);
    }

    // check
        
    /*
    cout << "triangular:" << endl;
    for (int i = 0; i < triangular.size(); i++)
    {
        Vector3i triangularCurrent = triangular[i];

        cout << triangularCurrent(0) << " " << triangularCurrent(1) << " " << triangularCurrent(2) << endl;
    }
    
    cout << "edge:" << endl;
    for (int i = 0; i < edge.size(); i++)
    {
        Vector2i edgeCurrent = edge[i];

        cout << edgeCurrent(0) << " " << edgeCurrent(1) << endl;
    }

    cout << "startnode:" << endl;
    for (int i = 0; i < startnode.size(); i++)
    {
        cout << startnode[i] << endl;
    }

    cout << "nodes:" << endl;
    for (int i = 0; i < nodes.size(); i++)
    {
        Vector3d xCurrent = nodes[i];

        cout << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
    }

    cout << "boundaryArray 1:" << endl;
    for (int i = 0; i < boundaryArray_1.size(); i++)
    {
        cout << boundaryArray_1[i] << endl;
    }

    cout << "boundaryArray 2:" << endl;
    for (int i = 0; i < boundaryArray_2.size(); i++)
    {
        cout << boundaryArray_2[i] << endl;
    }

    */
}