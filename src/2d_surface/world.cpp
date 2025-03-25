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
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");

	thickness = m_inputData.GetScalarOpt("thickness");
	Possion = m_inputData.GetScalarOpt("Possion");

	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");
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

    name << "datafiles/simDER";
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

	if ( timeStep % 10 != 0)
	{
		return;
	}

	//if (timeStep == Nstep)
	{
		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);

			outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << endl;
		}
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, thickness, Possion, deltaTime);

	plateBoundaryCondition();

	plate->setup();

	stepper = new timeStepper(*plate);

	// set up force
	m_inertialForce = new inertialForce(*plate, *stepper);
	m_gravityForce = new externalGravityForce(*plate, *stepper, gVector);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosity);
	m_elasticStretchingForce = new elasticStretchingForce(*plate, *stepper);
	m_elasticBendingForce = new elasticBendingForce(*plate, *stepper);
	
	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();
	m_elasticStretchingForce->setFirstJacobian();
	m_elasticBendingForce->setFirstJacobian();

	stepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;
}

void world::plateBoundaryCondition()
{
	for (int i = 0; i < plate->constraint.size(); i++)
	{
		double pos = plate->x[plate->constraint[i]];
		
		plate->setConstraint(pos, plate->constraint[i]);
	}
}

void world::updateTimeStep()
{
	Vector2d x1 = plate->getVertex(0);
	x1(1) = x1(1) - deltaTime * 0.05;
	plate->setVertexBoundaryCondition(x1, 0);

	Vector2d x2 = plate->getVertex(1);
	x2(1) = x2(1) - deltaTime * 0.05;
	plate->setVertexBoundaryCondition(x2, 1);


	bool goodSolved = false;

	while (goodSolved == false)
	{
		// Start with a trial solution for our solution x
		plate->updateGuess(); // x = x0 + u * dt

		updateEachStep();

		goodSolved = true;
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
}

void world::updateEachStep()
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		m_dampingForce->computeFd();
		m_gravityForce->computeFg();
		
		m_elasticStretchingForce->computeFs();
		m_elasticBendingForce->computeFb();
		
		normf = stepper->GlobalForceVec.norm();

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

		//cout << normf << " ";

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			m_gravityForce->computeJg();
			m_dampingForce->computeJd();

			m_elasticStretchingForce->computeJs();
			m_elasticBendingForce->computeJb();
			
			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			timeStep = Nstep;
			break;
		}
	}

	if (render)
	{
		cout << "iter " << iter << endl;
	}
}

int world::simulationRunning()
{

	if (timeStep < Nstep)
	{
		return 1;
	}
	else
	{
		return - 1;
	}
}

Vector2d world::getScaledCoordinate(int i, int j)
{
	Vector2d xCurrent;
	
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