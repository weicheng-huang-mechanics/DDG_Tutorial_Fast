#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");			
	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungRod = m_inputData.GetScalarOpt("YoungRod");
	YoungPlate = m_inputData.GetScalarOpt("YoungPlate");
	densityRod = m_inputData.GetScalarOpt("densityRod");
	densityPlate = m_inputData.GetScalarOpt("densityPlate");
	thickness = m_inputData.GetScalarOpt("thickness");
	Possion = m_inputData.GetScalarOpt("Possion");
	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");

	length = m_inputData.GetScalarOpt("length");
	width = m_inputData.GetScalarOpt("width");
	deltaLength = m_inputData.GetScalarOpt("deltaLength");

	xNum = m_inputData.GetIntOpt("xNum");			        // number of x rods
	yNum = m_inputData.GetIntOpt("yNum");			        // number of x rods
	rodLength = m_inputData.GetScalarOpt("rodLength");
	rodRadius = m_inputData.GetScalarOpt("rodRadius");
	xLength = m_inputData.GetScalarOpt("xLength");
	yLength = m_inputData.GetScalarOpt("yLength");
	nv = m_inputData.GetIntOpt("nv");

	cutIndex = m_inputData.GetIntOpt("cutIndex");

	baVector = m_inputData.GetVecOpt("baVector");
	brVector = m_inputData.GetVecOpt("brVector");
	muZero = m_inputData.GetScalarOpt("muZero");
	omega = m_inputData.GetScalarOpt("omega");

	stiffness = m_inputData.GetScalarOpt("stiffness");
	dBar = m_inputData.GetScalarOpt("dBar");

	mu = m_inputData.GetScalarOpt("mu");
	epsilonV = m_inputData.GetScalarOpt("epsilonV");

	shearRod = YoungRod/(2.0*(1.0+Possion));					// shear modulus

	deltaLen = rodLength / (nv - 1);

	numRods = xNum * yNum;

	deltaX = xLength / (xNum - 1);
	deltaY = yLength / (yNum - 1);

	rodDOF = 4 * nv - 1;

	if (render) 
	{
		cout << " deltaTime " << deltaTime << endl;
		cout << " totalTime " << totalTime << endl;
		cout << " YoungRod " << YoungRod << endl;
		cout << " YoungPlate " << YoungPlate << endl;
		cout << " densityRod " << densityRod << endl;
		cout << " densityPlate " << densityPlate << endl;
		cout << " Possion " << Possion << endl;
		cout << " thickness " << thickness << endl;
		cout << " length " << length << endl;
		cout << " width " << width << endl;
		cout << " deltaLength " << deltaLength << endl;
		cout << " xNum " << xNum << endl;
		cout << " yNum " << yNum << endl;
		cout << " rodLength " << rodLength << endl;
		cout << " xLength " << xLength << endl;
		cout << " yLength " << yLength << endl;
		cout << " nv " << nv << endl;
		cout << " baVector " << baVector.transpose() << endl;
		cout << " brVector " << brVector.transpose() << endl;
		cout << " muZero " << muZero << endl;
		cout << " omega " << omega << endl;
		cout << " stiffness " << stiffness << endl;
		cout << " dBar " << dBar << endl;
		cout << " mu " << mu << endl;
		cout << " epsilonV " << epsilonV << endl;
	}
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDiscretePlate";
    name << "_mu_" << mu;
    name << "_omega_" << omega;
    name << "_cutIndex_" << cutIndex;
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if ( timeStep % 10 != 0)
	{
		//return;
	}

	/*

	int plateIndex = 0;
	for (int i = 0; i < plate->nv; i++)
	{
		Vector3d xCurrent = plate->getVertex(i);

		outfile << currentTime << " " << plateIndex << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) - deltaLen << endl;
	}

	for (int n = 0; n < numRods; n++)
	{
		for (int i = 0; i < rodsVector[n]->nv; i++)
		{
			Vector3d xCurrent = rodsVector[n]->getVertex(i);
			outfile << currentTime << " " << n+1 << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
		}
	}
	*/

	//if (timeStep == Nstep)
	{
		Vector3d xCurrent1 = plate->getVertex(5);
		Vector3d xCurrent2 = plate->getVertex(257);

		//Vector3d uMid = plate->getVelocityOld(131);

		outfile << currentTime << " " << xCurrent1(0) << " " << xCurrent1(1) << " " << xCurrent2(0) << " " << xCurrent2(1) << endl;
	}

	//Vector3d xCurrent = plate->getVertex(130);

	//outfile << currentTime << " " << xCurrent(1) + 0.2 - 0.104 << endl;

	/*

	//if (timeStep == Nstep)
	{
		int plateIndex = 0;
		for (int i = 0; i < plate->nv; i++)
		{
			Vector3d xCurrent = plate->getVertex(i);

			outfile << currentTime << " " << plateIndex << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) - deltaLen << endl;
		}

		ofstream out1("datafiles/edge.txt");
		for (int i = 0; i < plate->edgeOutput.size(); i++)
		{
			Vector2i edgeCurrent = plate->edgeOutput[i];
			out1 << edgeCurrent(0) << " " << edgeCurrent(1) << endl;
		}
        out1.close();

        ofstream out2("datafiles/triangular.txt");
		for (int i = 0; i < plate->triangularOutput.size(); i++)
		{
			Vector3i triangularCurrent = plate->triangularOutput[i];
			out2 << triangularCurrent(0) << " " << triangularCurrent(1) << " " << triangularCurrent(2) << endl;
		}
        out2.close();

		ofstream out3("datafiles/rod.txt");
		for (int n = 0; n < numRods; n++)
		{
			for (int i = 0; i < rodsVector[n]->nv; i++)
			{
				Vector3d xCurrent = rodsVector[n]->getVertex(i);
				outfile << currentTime << " " << n+1 << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
			}
		}
		out3.close();
	}

	if (timeStep == Nstep)
	{
		ofstream out4("datafiles/coupleElement.txt");
		for (int i = 0; i < v_coupleBending.size(); i++)
		{
			out4 << v_coupleBending[i].rodIndex << " ";
			out4 << v_coupleBending[i].triangularIndex << " ";
			out4 << v_coupleBending[i].alpha << " ";
			out4 << v_coupleBending[i].beta << " ";
			out4 << v_coupleBending[i].gamma << endl;
		}
        out4.close();
	}

	*/
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungPlate, densityPlate, thickness, Possion, deltaTime, length, width, deltaLength, dBar + rodLength+deltaLen);

	plateBoundaryCondition();

	plate->setup();

	stepperPlate = new timeStepperPlate(*plate);

	// set up force
	m_inertialForcePlate = new inertialForcePlate(*plate, *stepperPlate);
	m_gravityForcePlate = new externalGravityForcePlate(*plate, *stepperPlate, gVector);
	m_stretchForcePlate = new elasticStretchingForcePlate(*plate, *stepperPlate);
	m_bendingForcePlate = new elasticBendingForcePlate(*plate, *stepperPlate);
	m_dampingForcePlate = new dampingForcePlate(*plate, *stepperPlate, viscosity);
	m_externalContactForcePlate = new externalContactForcePlate(*plate, *stepperPlate, stiffness, dBar, mu, epsilonV);

	plate->updateTimeStep();

	rolIndex = VectorXd::Zero(numRods);
	rowIndex = VectorXd::Zero(numRods);

	for (int i=0; i < numRods; i++)
	{
		// Set up geometry
		rodGeometry(i);	

		// Create the rod 
		rodsVector.push_back( new elasticRod(vertices, vertices, densityRod, rodRadius, deltaTime,
			YoungRod, shearRod, localLength, theta) );

		// Set up boundary condition
		rodBoundaryCondition(i);
	
		// setup the rod so that all the relevant variables are populated
		rodsVector[i]->setup();
		// End of rod setup

		// set up the time stepper
		stepperVector.push_back( new timeStepper(*rodsVector[i]) );

		// declare the forces
		v_stretchForce.push_back( new elasticStretchingForce( *rodsVector[i], *stepperVector[i]) );
		v_bendingForce.push_back( new elasticBendingForce( *rodsVector[i], *stepperVector[i]) );
		v_twistingForce.push_back( new elasticTwistingForce( *rodsVector[i], *stepperVector[i]) );
		v_inertialForce.push_back( new inertialForce( *rodsVector[i], *stepperVector[i]) );
		v_gravityForce.push_back( new externalGravityForce( *rodsVector[i], *stepperVector[i], gVector) );
		v_externalMagneticForce.push_back(new externalMagneticForce(*rodsVector[i], *stepperVector[i], baVector, brVector, muZero, omega) );
		v_dampingForce.push_back(new dampingForce(*rodsVector[i], *stepperVector[i], viscosity));
		v_externalContactForce.push_back(new externalContactForce(*rodsVector[i], *stepperVector[i], stiffness, dBar, mu, epsilonV));

		// Allocate every thing to prepare for the first iteration
		rodsVector[i]->updateTimeStep();
	}

	v_coupleElement.clear();
	v_coupleBending.clear();

	for (int n = 0; n < numRods; n++)
	{
		Vector3d xStart = rodsVector[n]->getVertex(0);

		Vector3d xCurrent;

		double minDistance = 1000.0;
		double localDistance;

		int minIndex;

		for (int i = 0; i < plate->v_triangularElement.size(); i++)
		{
			Vector3d x1 = plate->v_triangularElement[i].x_1;
			Vector3d x2 = plate->v_triangularElement[i].x_2;
			Vector3d x3 = plate->v_triangularElement[i].x_3;

			localDistance = (xStart - x1).norm() + (xStart - x2).norm() + (xStart - x3).norm();

			if (localDistance < minDistance)
			{
				minDistance = localDistance;
				minIndex = i;
			}
		}

		coupleElement m_coupleElement;

		m_coupleElement.rodIndex = n;
		m_coupleElement.plateIndex = plate->v_triangularElement[minIndex].nv_1;
		xCurrent = plate->getVertex(m_coupleElement.plateIndex);
		m_coupleElement.edge = (xStart - xCurrent).norm();
		v_coupleElement.push_back(m_coupleElement);

		m_coupleElement.rodIndex = n;
		m_coupleElement.plateIndex = plate->v_triangularElement[minIndex].nv_2;
		xCurrent = plate->getVertex(m_coupleElement.plateIndex);
		m_coupleElement.edge = (xStart - xCurrent).norm();
		v_coupleElement.push_back(m_coupleElement);

		m_coupleElement.rodIndex = n;
		m_coupleElement.plateIndex = plate->v_triangularElement[minIndex].nv_3;
		xCurrent = plate->getVertex(m_coupleElement.plateIndex);
		m_coupleElement.edge = (xStart - xCurrent).norm();
		v_coupleElement.push_back(m_coupleElement);

		coupleBending m_coupleBending;
		m_coupleBending.rodIndex = n;
		m_coupleBending.triangularIndex = minIndex;
		v_coupleBending.push_back(m_coupleBending);
	}

	for (int i = 0; i < v_coupleBending.size(); i++)
	{
		int rodIndex = v_coupleBending[i].rodIndex;
		int triangularIndex = v_coupleBending[i].triangularIndex;

		int nv1 = plate->v_triangularElement[triangularIndex].nv_1;
		int nv2 = plate->v_triangularElement[triangularIndex].nv_2;
		int nv3 = plate->v_triangularElement[triangularIndex].nv_3;

		Vector3d x1 = plate->getVertex(nv1);
		Vector3d x2 = plate->getVertex(nv2);
		Vector3d x3 = plate->getVertex(nv3);

		Vector3d e_1 = x3 - x1;
		Vector3d e_2 = x2 - x1;

		Vector3d n_1 = e_1.cross(e_2);

		Vector3d x4 = rodsVector[rodIndex]->getVertex(0);
		Vector3d x5 = rodsVector[rodIndex]->getVertex(1);

		Vector3d n_2 = x5 - x4;

		v_coupleBending[i].nbar = n_1 / n_1.norm() - n_2 / n_2.norm();

		if ( (v_coupleBending[i].nbar).norm() > 1.5 )
		{
			x2 = plate->getVertex(nv3);
			x3 = plate->getVertex(nv2);

			e_1 = x3 - x1;
			e_2 = x2 - x1;

			n_1 = e_1.cross(e_2);

			v_coupleBending[i].nbar = n_1 / n_1.norm() - n_2 / n_2.norm();

			v_coupleBending[i].ifInverse = 1;
		}
		else
		{
			v_coupleBending[i].ifInverse = 0;
		}
	}

	for (int i = 0; i < v_coupleBending.size(); i++)
	{
		int rodIndex = v_coupleBending[i].rodIndex;
		int triangularIndex = v_coupleBending[i].triangularIndex;

		int nv1 = plate->v_triangularElement[triangularIndex].nv_1;
		int nv2 = plate->v_triangularElement[triangularIndex].nv_2;
		int nv3 = plate->v_triangularElement[triangularIndex].nv_3;

		Vector3d x1 = plate->getVertex(nv1);
		Vector3d x2 = plate->getVertex(nv2);
		Vector3d x3 = plate->getVertex(nv3);

		Matrix3d matrixA;
		matrixA.col(0) = x1;
		matrixA.col(1) = x2;
		matrixA.col(2) = x3;

		Vector3d vectorB;

		vectorB = rodsVector[rodIndex]->getVertex(0);

		vectorB(2) = plate->startHeight;

		Vector3d coeff = matrixA.colPivHouseholderQr().solve(vectorB);

		v_coupleBending[i].alpha = coeff(0);
		v_coupleBending[i].beta = coeff(1);
		v_coupleBending[i].gamma = coeff(2);

	}

	total_dof_mapped = plate->ndof + numRods * rodDOF;

	isConstrained_global = VectorXi::Zero(total_dof_mapped);

	setGlobalBoundary();

	m_GlobalStepper = new GlobalStepper(uncons_dof_mapped, isConstrained_global, fullToUnconsMap_global);

	m_elasticStretchingForceCouple = new elasticStretchingForceCouple(*m_GlobalStepper, 1000 * YoungRod * rodRadius, rodDOF, plate->ndof);
	m_elasticBendingForceCouple = new elasticBendingForceCouple(*m_GlobalStepper, 1000 * YoungRod * rodRadius * rodRadius * rodRadius, rodDOF, plate->ndof);

	// set up first jacobian
	m_inertialForcePlate->setFirstJacobian();
	m_stretchForcePlate->setFirstJacobian();
	m_bendingForcePlate->setFirstJacobian();
	m_dampingForcePlate->setFirstJacobian();

	for (int i = 0; i < plate->ndof; i++)
	{
		for (int j = 0; j < plate->ndof; j++)
		{
			double jacobian_element = stepperPlate->JacobianMatrix(i, j);

			if (jacobian_element != 0)
			{
				m_GlobalStepper->addJacobian(i, j, jacobian_element);
			}
		}
	}

	for (int n = 0; n < numRods; n++)
	{
		v_inertialForce[n]->setFirstJacobian();
		v_stretchForce[n]->setFirstJacobian();
		v_bendingForce[n]->setFirstJacobian();
		v_twistingForce[n]->setFirstJacobian();

		for (int i = 0; i < rodsVector[n]->ndof; i++)
		{
			for (int j = 0; j < rodsVector[n]->ndof; j++)
			{
				int local_dof1 = plate->ndof + n * rodDOF + i;
				int local_dof2 = plate->ndof + n * rodDOF + j;

				double jacobian_element = stepperVector[n]->JacobianMatrix(i, j);

				if (jacobian_element != 0)
				{
					m_GlobalStepper->addJacobian(local_dof1, local_dof2, jacobian_element);
				}
			}
		}
	}

	for (int i = 0; i < v_coupleElement.size(); i++)
	{
		int rodIndex = v_coupleElement[i].rodIndex;
		int plateIndex = v_coupleElement[i].plateIndex;

		m_elasticStretchingForceCouple->setFirstJacobian(rodIndex, plateIndex);
	}

	for (int i = 0; i < v_coupleBending.size(); i++)
	{
		int rodIndex = v_coupleBending[i].rodIndex;
		int triangularIndex = v_coupleBending[i].triangularIndex;

		int nv1 = plate->v_triangularElement[triangularIndex].nv_1;
		int nv2 = plate->v_triangularElement[triangularIndex].nv_2;
		int nv3 = plate->v_triangularElement[triangularIndex].nv_3;

		m_elasticBendingForceCouple->setFirstJacobian(nv1, nv2, nv3, rodIndex);
	}

	m_GlobalStepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;
}

// Setup geometry
void world::rodGeometry(int i)
{
	vertices = MatrixXd::Zero(nv, 3);

	Vector3d xMid = plate->getVertex(plate->nv / 2);

	int rolNum = floor( i / xNum);
	int colNum = i - rolNum * xNum;

	rolIndex(i) = rolNum * 1.0;

	rolNum = rolNum - (yNum - 1) / 2;
	colNum = colNum - (xNum - 1) / 2;

	rowIndex(i) = colNum * 1.0;

	localLength = rodLength;

	for (int k = 0; k < nv; k++)
	{
		vertices(k, 0) = colNum * deltaX + xMid(0);
		vertices(k, 1) = rolNum * deltaY + xMid(1) - 0.001;
		vertices(k, 2) = dBar + rodLength + deltaLen - deltaLen * (k+1);
	}

    // initial theta should be zeros
    theta = VectorXd::Zero(nv - 1);

    cout << rolNum << endl;

    if (rolNum == 100)
    {
    	if (colNum == -2 || colNum == 2 || colNum == 0)
    	{
    		vertices = MatrixXd::Zero(2, 3);

    		localLength = 2 * deltaLen;

    		for (int k = 0; k < 2; k++)
			{
				vertices(k, 0) = colNum * deltaX + xMid(0);
				vertices(k, 1) = rolNum * deltaY + xMid(1) - 0.001;
				vertices(k, 2) = dBar + rodLength + deltaLen - deltaLen * (k+1);
			}

			theta = VectorXd::Zero(1);
    	}
    }
}

void world::rodBoundaryCondition(int n)
{
	//rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(0), 0);
	//rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(1), 1);

	Vector3d xCurrent1 = plate->getVertex(5);
	Vector3d xCurrent2 = plate->getVertex(257);

	Vector3d tangent = (xCurrent2 - xCurrent1) / (xCurrent2 - xCurrent1).norm();

	tangent(2) = 0.0;

	tangent = tangent / tangent.norm();

	Vector3d nHead;

	nHead(0) = 0.0;
	nHead(1) = 1.0;
	nHead(2) = 0.0;

	double cosTheta = tangent.dot(nHead) / (tangent.norm() * nHead.norm());

	double turnTheta = acos(cosTheta);

	if (n == 0)
	{
		//cout << turnTheta << endl;
	}

	//cout << turnTheta << endl;

	for (int i = 0; i < rodsVector[n]->ne; i++)
	{
		rodsVector[n]->setThetaBoundaryCondition(- turnTheta, i);
	}

	//rodsVector[n]->setThetaBoundaryCondition(turnTheta, 0);

	for (int i = 0; i < rodsVector[n]->nv; i++)
	{
		Vector3d xCurrent = rodsVector[n]->getVertex(i);

		//rodsVector[n]->setOneVertexBoundaryCondition(xCurrent(0), i, 0);
	}
}

void world::plateBoundaryCondition()
{
	for (int i = 0; i < plate->boundaryArray_1.size(); i++)
	{
		int ind = plate->boundaryArray_1[i];
		
		//plate->setVertexBoundaryCondition(plate->getVertex(ind), ind);
	}


	for (int i = 0; i < plate->nv; i++)
	{
		Vector3d xCurrent = plate->getVertex(i);

		//plate->setVertexOneBoundaryCondition(xCurrent(0), i, 0);
	}
}

void world::computeForce()
{
	plate->prepareForIteration();

	stepperPlate->setZero();

	m_inertialForcePlate->computeFi();
	m_gravityForcePlate->computeFg();
	m_stretchForcePlate->computeFs();
	m_bendingForcePlate->computeFb();
	m_dampingForcePlate->computeFd();
	m_externalContactForcePlate->computeFc();

	for (int i = 0; i < plate->ndof; i++)
	{
		m_GlobalStepper->addForce(i, stepperPlate->ForceVector(i));
	}

	for (int n = 0; n < numRods; n++)
	{
		rodsVector[n]->prepareForIteration();

		stepperVector[n]->setZero();

		// Compute the forces
		v_inertialForce[n]->computeFi();
		v_stretchForce[n]->computeFs();	
		v_bendingForce[n]->computeFb();
		v_twistingForce[n]->computeFt();
		v_gravityForce[n]->computeFg();
		v_dampingForce[n]->computeFd();
		v_externalContactForce[n]->computeFc();
		v_externalMagneticForce[n]->computeFm(currentTime, rolIndex(n)/(yNum-1), rowIndex(n), cutIndex);

		for (int i = 0; i < rodsVector[n]->ndof; i++)
		{
			int local_dof = plate->ndof + n * rodDOF + i;

			if ( getIfConstrainGlobal(local_dof) == 0 )
			{
				m_GlobalStepper->addForce(local_dof, stepperVector[n]->ForceVector(i));
			}
		}
	}

	for (int i = 0; i < v_coupleElement.size(); i++)
	{
		int rodIndex = v_coupleElement[i].rodIndex;
		int plateIndex = v_coupleElement[i].plateIndex;

		Vector3d rodNode = rodsVector[rodIndex]->getVertex(0);
		Vector3d plateNode = plate->getVertex(plateIndex);

		double refLength = v_coupleElement[i].edge;

		m_elasticStretchingForceCouple->computeFs(rodIndex, plateIndex, rodNode, plateNode, refLength);
	}

	for (int i = 0; i < v_coupleBending.size(); i++)
	{
		int rodIndex = v_coupleBending[i].rodIndex;
		int triangularIndex = v_coupleBending[i].triangularIndex;

		int nv1;
		int nv2;
		int nv3;

		if (v_coupleBending[i].ifInverse == 0)
		{
			nv1 = plate->v_triangularElement[triangularIndex].nv_1;
			nv2 = plate->v_triangularElement[triangularIndex].nv_2;
			nv3 = plate->v_triangularElement[triangularIndex].nv_3;
		}

		if (v_coupleBending[i].ifInverse == 1)
		{
			nv1 = plate->v_triangularElement[triangularIndex].nv_1;
			nv2 = plate->v_triangularElement[triangularIndex].nv_3;
			nv3 = plate->v_triangularElement[triangularIndex].nv_2;
		}

		Vector3d x1 = plate->getVertex(nv1);
		Vector3d x2 = plate->getVertex(nv2);
		Vector3d x3 = plate->getVertex(nv3);

		Vector3d x4 = rodsVector[rodIndex]->getVertex(0);
		Vector3d x5 = rodsVector[rodIndex]->getVertex(1);

		Vector3d nBar = v_coupleBending[i].nbar;

		m_elasticBendingForceCouple->computeFb(nv1, nv2, nv3, x1, x2, x3, rodIndex, x4, x5, nBar);
	}
}


void world::computeJacobian()
{
	m_inertialForcePlate->computeJi();
	m_gravityForcePlate->computeJg();
	m_stretchForcePlate->computeJs();
	m_bendingForcePlate->computeJb();
	m_dampingForcePlate->computeJd();
	m_externalContactForcePlate->computeJc();

	for (int i = 0; i < plate->ndof; i++)
	{
		for (int j = 0; j < plate->ndof; j++)
		{
			m_GlobalStepper->addJacobian(i, j, stepperPlate->JacobianMatrix(i, j));
		}
	}

	for(int n = 0; n < numRods; n++)
	{
		v_inertialForce[n]->computeJi();
		v_stretchForce[n]->computeJs();
		v_bendingForce[n]->computeJb();
		v_twistingForce[n]->computeJt();
		v_gravityForce[n]->computeJg();
		v_dampingForce[n]->computeJd();
		v_externalContactForce[n]->computeJc();

		for (int i = 0; i < rodsVector[n]->ndof; i++)
		{
			for (int j = 0; j < rodsVector[n]->ndof; j++)
			{
				int local_dof1 = plate->ndof + n * rodDOF + i;
				int local_dof2 = plate->ndof + n * rodDOF + j;

				m_GlobalStepper->addJacobian(local_dof1, local_dof2, stepperVector[n]->JacobianMatrix(i, j) );
			}
		}
	}

	for (int i = 0; i < v_coupleElement.size(); i++)
	{
		int rodIndex = v_coupleElement[i].rodIndex;
		int plateIndex = v_coupleElement[i].plateIndex;

		Vector3d rodNode = rodsVector[rodIndex]->getVertex(0);
		Vector3d plateNode = plate->getVertex(plateIndex);

		double refLength = v_coupleElement[i].edge;

		m_elasticStretchingForceCouple->computeJs(rodIndex, plateIndex, rodNode, plateNode, refLength);
	}

	for (int i = 0; i < v_coupleBending.size(); i++)
	{
		int rodIndex = v_coupleBending[i].rodIndex;
		int triangularIndex = v_coupleBending[i].triangularIndex;

		int nv1;
		int nv2;
		int nv3;

		if (v_coupleBending[i].ifInverse == 0)
		{
			nv1 = plate->v_triangularElement[triangularIndex].nv_1;
			nv2 = plate->v_triangularElement[triangularIndex].nv_2;
			nv3 = plate->v_triangularElement[triangularIndex].nv_3;
		}

		if (v_coupleBending[i].ifInverse == 1)
		{
			nv1 = plate->v_triangularElement[triangularIndex].nv_1;
			nv2 = plate->v_triangularElement[triangularIndex].nv_3;
			nv3 = plate->v_triangularElement[triangularIndex].nv_2;
		}

		Vector3d x1 = plate->getVertex(nv1);
		Vector3d x2 = plate->getVertex(nv2);
		Vector3d x3 = plate->getVertex(nv3);

		Vector3d x4 = rodsVector[rodIndex]->getVertex(0);
		Vector3d x5 = rodsVector[rodIndex]->getVertex(1);

		Vector3d nBar = v_coupleBending[i].nbar;

		m_elasticBendingForceCouple->computeJb(nv1, nv2, nv3, x1, x2, x3, rodIndex, x4, x5, nBar);
	}
}

void world::updataMotion()
{
	local_vec_motion = VectorXd::Zero(plate->uncons);

	for (int i = 0; i < plate->ndof; i++)
	{
		if ( getIfConstrainGlobal(i) == 0 )
		{
			local_vec_motion(plate->fullToUnconsMap[i]) = m_GlobalStepper->GlobalMotionVec(fullToUnconsMap_global[i]);
		}
	}

	plate->updateNewtonMethod(local_vec_motion);

	for (int n = 0; n < numRods; n++)
	{
		local_vec_motion = VectorXd::Zero(rodsVector[n]->uncons);

		for (int i = 0; i < rodsVector[n]->ndof; i++)
		{
			int local_dof = plate->ndof + n * rodDOF + i;

			if ( getIfConstrainGlobal(local_dof) == 0 )
			{
				local_vec_motion(rodsVector[n]->fullToUnconsMap[i]) = m_GlobalStepper->GlobalMotionVec(fullToUnconsMap_global[local_dof]);
			}
		}

		rodsVector[n]->updateNewtonMethod(local_vec_motion);
	}
}

void world::updateTimeStep()
{
	//cout << "m1: " << rodsVector[0]->m1.row(0).transpose() << endl;
	//cout << "m2: " << rodsVector[0]->m2.row(0).transpose() << endl;
	//cout << "ta: " << rodsVector[0]->tangent.row(0).transpose() << endl;

	//cout << "m1 " << (rodsVector[0]->m1.row(rodsVector[0]->ne-1)).transpose() << endl;

	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	for (int n = 0; n < numRods; n++)
	{
		rodBoundaryCondition(n);
		rodsVector[n]->updateGuess();
	}

	plate->updateGuess(); // x = x0 + u * dt

	while (solved == false)
	{
		m_GlobalStepper->setZero();

		computeForce();

		normf = m_GlobalStepper->GlobalForceVec.norm();

		if (iter == 0)
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		if (solved == false)
		{
			computeJacobian();

			m_GlobalStepper->integrator();

    		updataMotion();

    		iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}

		if (render)
		{
			cout << "time: " << currentTime << endl;
		}
	}

	plate->updateTimeStep();

	for(int n = 0; n < numRods; n++)
	{
		rodsVector[n]->updateTimeStep();
	}

	currentTime += deltaTime;
		
	timeStep++;
}

int world::simulationRunning()
{
	if (timeStep < Nstep) 
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

Vector3d world::getScaledCoordinatePlate(int i, int j)
{
	Vector3d xCurrent;
	
	if (j == 0)
	{
		xCurrent = plate->v_edgeElement[i].x_1 * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = plate->v_edgeElement[i].x_2 * scaleRendering;
	}

	return xCurrent;
}

int world::numStretchingPair()
{
	return plate->edgeNum;
}

int world::numPoints(int n)
{
	return rodsVector[n]->nv;
}

double world::getScaledCoordinate(int n, int i)
{
	return rodsVector[n]->x[i] * scaleRendering;
}

Vector3d world::getScaledCoordinatem1(int n)
{
	return rodsVector[n]->m1.row(rodsVector[n]->ne-1);
}

Vector3d world::getScaledCoordinateCouple(int n, int i)
{
	if (i == 0)
	{
		int index = v_coupleElement[n].rodIndex;

		return rodsVector[index]->getVertex(0) * scaleRendering;
	}

	if (i == 1)
	{
		int index = v_coupleElement[n].plateIndex;

		return plate->getVertex(index) * scaleRendering;
	}
}

int world::getnumRods()
{
	return numRods;
}

Vector3d world::getMidPoint()
{
	Vector3d xMid = plate->getVertex(plate->nv / 2);

	return xMid * scaleRendering;
}

int world::getIfConstrainGlobal(int k)
{
	return isConstrained_global[k];
}

void world::setGlobalBoundary()
{
	cons_dof_mapped = 0;

	for (int i = 0; i < plate->ndof; i++)
	{
		if (plate->getIfConstrained(i) == 1)
		{
			cons_dof_mapped = cons_dof_mapped + 1;
			isConstrained_global[i] = 1;
		}
	}

	for (int i = 0; i < numRods; i++)
	{
		for (int j = 0; j < rodsVector[i]->ndof; j++)
		{
			int localIndex = plate->ndof + i * rodDOF + j;

			if (rodsVector[i]->getIfConstrained(j) == 1)
			{
				cons_dof_mapped = cons_dof_mapped + 1;
				isConstrained_global[localIndex] = 1;
			}
		}
	}

	uncons_dof_mapped = total_dof_mapped - cons_dof_mapped;

	unconstrainedMap_global = VectorXi::Zero(uncons_dof_mapped); // maps xUncons to x
	fullToUnconsMap_global = VectorXi::Zero(total_dof_mapped);

	int c = 0;
	for (int i = 0; i < total_dof_mapped; i++)
	{
		if (isConstrained_global[i] == 0)
		{
			unconstrainedMap_global[c] = i;
			fullToUnconsMap_global[i] = c;
			c++;
		}
	}
}

