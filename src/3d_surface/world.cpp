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

    name << "datafiles/simDiscretePlate";
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
		return;
	}

	//if (timeStep == Nstep)
	{
		
		for (int i = 0; i < plate->nv; i++)
		{
			Vector3d xCurrent = plate->getVertex(i);

			outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
		}
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
	m_stretchForce = new elasticStretchingForce(*plate, *stepper);
	m_bendingForce = new elasticBendingForce(*plate, *stepper);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosity);
	
	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_stretchForce->setFirstJacobian();
	m_bendingForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();

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

	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	// Start with a trial solution for our solution x
	plate->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		m_gravityForce->computeFg();
		m_stretchForce->computeFs();
		m_bendingForce->computeFb();
		m_dampingForce->computeFd();

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

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			m_gravityForce->computeJg();
			m_stretchForce->computeJs();
			m_bendingForce->computeJb();
			m_dampingForce->computeJd();

			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " iter=" << iter << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		timeStep = Nstep; // we are exiting
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
		return -1;
	}
}

Vector3d world::getScaledCoordinate(int i, int j)
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