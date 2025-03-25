#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include <fstream>

struct edgeElement
{
	int nv_1;
	int nv_2;

	Vector3d x_1;
	Vector3d x_2;

	double refLength;
	double edgeLength;

	VectorXi arrayNum;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_radius, 
		double m_Possion, double m_dt);
	~elasticPlate();

	double YoungM;
	double radius;
	double Possion;
	double dt;
	double density;

	Vector3d getVertex(int i);
	Vector3d getVertexOld(int i);
	Vector3d getVelocity(int i);
	Vector3d getVertexStart(int i);
	double getTheta(int i);

	VectorXd x;
	VectorXd x0;
	VectorXd u;
	VectorXd x_initial;

	std::vector<Vector3d> v_nodes;
    std::vector<Vector2i> edge;
    std::vector<int> constraint;

	std::vector<edgeElement> v_edgeElement;
	
	int temp;

	double crossSectionalArea;

	int nv;
	int ne;
	int edgeNum;

	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector3d position, int k);
	void setOneVertexBoundaryCondition(double position, int i, int k);
	void setConstraint(double position, int k);

	void computeEdge();

	// boundary conditions
	int* isConstrained;
	int getIfConstrained(int k);
	int* unconstrainedMap;
	int* fullToUnconsMap;
	void setup();
	void setupMap();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
	void prepareForIteration();

	VectorXd massArray;
	void setupMass();
	
	double EA;
	
	void updateEdgePair();
	
	private:
};

#endif
