#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include input file and option
#include "setInput.h"

// include elastic Plate class
#include "elasticPlate.h"

// include time stepper
#include "timeStepperPlate.h"

// include force
#include "inertialForcePlate.h"
#include "externalGravityForcePlate.h"
#include "elasticStretchingForcePlate.h"
#include "elasticBendingForcePlate.h"
#include "dampingForcePlate.h"
#include "externalContactForcePlate.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"
#include "externalMagneticForce.h"
#include "dampingForce.h"
#include "timeStepper.h"
#include "externalContactForce.h"

#include "GlobalStepper.h"

#include "elasticStretchingForceCouple.h"
#include "elasticBendingForceCouple.h"

struct coupleElement
{
	int rodIndex;
	int plateIndex;

	double edge;
};

struct coupleBending
{
	int rodIndex;
	int triangularIndex;

	Vector3d nbar;

	int ifInverse;

	double alpha;
	double beta;
	double gamma;
};

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
	Vector3d getScaledCoordinatePlate(int i, int j);

	int numPoints(int n);
	double getScaledCoordinate(int n, int i);
	int getnumRods();

	Vector3d getMidPoint();

	Vector3d getScaledCoordinateCouple(int n, int i);

	Vector3d getScaledCoordinatem1(int n);
		
private:

	std::vector<coupleElement> v_coupleElement;
	std::vector<coupleBending> v_coupleBending;

	// physical parameters
	bool render;
	bool saveData;
	double deltaTime;
	double totalTime;
	double YoungRod;
	double YoungPlate;
	double densityRod;
	double densityPlate;
	double thickness;
	double Possion;
	double stol;
	double forceTol;
	double scaleRendering;
	int maxIter;
	Vector3d gVector;
	double viscosity;

	double length;
	double width;
	double deltaLength;

	double stiffness;
	double dBar;

	double shearRod;
	int xNum;
	int yNum;
	double rodLength;
	double rodRadius;
	double xLength;
	double yLength;
	int nv;
	double deltaX;
	double deltaY;
	int numRods;
	double deltaLen;

	Vector3d baVector;
	Vector3d brVector;
	double muZero;
	double omega;

	double mu;
	double epsilonV;

	// Geometry
	MatrixXd vertices;
	VectorXd theta;

	int rodDOF;

	int Nstep;
	int timeStep;

	int cutIndex;

	double characteristicForce;

	double currentTime;

	void plateBoundaryCondition();

	// Plate
	elasticPlate *plate;

	// stepper
	timeStepperPlate *stepperPlate;

	// force
	inertialForcePlate *m_inertialForcePlate;
	externalGravityForcePlate *m_gravityForcePlate;
	elasticStretchingForcePlate *m_stretchForcePlate;
	elasticBendingForcePlate *m_bendingForcePlate;
	dampingForcePlate *m_dampingForcePlate;
	externalContactForcePlate *m_externalContactForcePlate;

	// Rod
	std::vector<elasticRod*> rodsVector;
	
	// set up the time stepper
	std::vector<timeStepper*> stepperVector;
	
	// declare the forces
	std::vector<elasticStretchingForce*> v_stretchForce;
	std::vector<elasticBendingForce*> v_bendingForce;
	std::vector<elasticTwistingForce*> v_twistingForce;
	std::vector<inertialForce*> v_inertialForce;
	std::vector<externalGravityForce*> v_gravityForce;
	std::vector<externalMagneticForce*> v_externalMagneticForce;
	std::vector<dampingForce*> v_dampingForce;
	std::vector<externalContactForce*> v_externalContactForce;

	GlobalStepper* m_GlobalStepper;
	elasticStretchingForceCouple* m_elasticStretchingForceCouple;
	elasticBendingForceCouple* m_elasticBendingForceCouple;

	void rodGeometry(int i);
	void rodBoundaryCondition(int i);

	int xTemp;
	int yTemp;

	int total_dof_mapped;
	int uncons_dof_mapped;
	int cons_dof_mapped;

	VectorXi isConstrained_global;
	VectorXi unconstrainedMap_global;
	VectorXi fullToUnconsMap_global;

	void setGlobalBoundary();

	int getIfConstrainGlobal(int k);

	void computeForce();
	void computeJacobian();
	void updataMotion();

	VectorXd local_vec_motion;

	VectorXd rolIndex;
	VectorXd rowIndex;


	double localLength;
};

#endif
