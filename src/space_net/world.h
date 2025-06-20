#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include input file and option
#include "setInput.h"

// include elastic Plate class
#include "elasticPlate.h"

// include time stepper
#include "timeStepper.h"

// include force
#include "inertialForce.h"
#include "externalGravityForce.h"
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "dampingForce.h"
#include "externalContactForce.h"
#include "externalContactForceBody.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	void setPlateStepper();

	void updateTimeStep();

	int simulationRunning();

	int numStretchingPair();
	Vector3d getScaledCoordinate(int i, int j);
	Vector3d getScaledCoordinatePoint();

	int getBoundaryEdge();

	double getBoxSize();
		
private:

	// physical parameters
	bool render;
	bool saveData;
	double deltaTime;
	double totalTime;
	double YoungM;
	double density;
	double Possion;
	double stol;
	double forceTol;
	double scaleRendering;
	int maxIter;
	Vector3d gVector;
	double viscosity;
	double rodRadius;
	double deltaLength;
	double boxSize;

	double dBarBody;
	double stiffnessBody;

	double stiffness;
    double dBar; 

	int Nstep;
	int timeStep;

	double characteristicForce;

	double currentTime;

	void plateBoundaryCondition();

	// Plate
	elasticPlate *plate;

	// stepper
	timeStepper *stepper;

	// force
	inertialForce *m_inertialForce;
	externalGravityForce *m_gravityForce;
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	dampingForce *m_dampingForce;
	externalContactForce *m_externalContactForce;
	externalContactForceBody *m_externalContactForceBody;

	void updateEachStep();
};

#endif
