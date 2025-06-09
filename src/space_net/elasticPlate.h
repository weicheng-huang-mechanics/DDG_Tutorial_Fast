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

struct bendingElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	VectorXi arrayNum;

	Vector3d x_1;
	Vector3d x_2;
	Vector3d x_3;

	Vector3d e_1;
	Vector3d e_2;

	double norm_1;
	double norm_2;

	Vector3d t_1;
	Vector3d t_2;

	Vector3d nBar;

	double voroniLength;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_radius, 
		double m_Possion, double m_dt, double m_deltaLength);
	~elasticPlate();

	double YoungM;
	double radius;
	double Possion;
	double dt;
	double density;

	double deltaEdge;

	Vector3d getVertex(int i);
	Vector3d getVertexOld(int i);
	Vector3d getVelocity(int i);

	VectorXd x;
	VectorXd x0;
	VectorXd u;

	std::vector<Vector3d> v_nodes;
    std::vector<Vector2i> edge;
    std::vector<Vector3i> bending;
    std::vector<int> constraint;

    std::vector<VectorXi> boundaryRod;
    std::vector<VectorXi> ringIndex;

    MatrixXi boundaryRodMatrix;
    MatrixXi ringIndexMatrix;

	std::vector<edgeElement> v_edgeElement;
	std::vector<bendingElement> v_bendingElement;

	int temp;

	int nv;
	int edgeNum;
	int bendingNum;

	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector3d position, int k);

	void computeEdge();
	void computeBending();

	void updateEdgePair();
	void updateBendingPair();

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
	VectorXi boundaryIndex;

	double EA;
	double EI;

	VectorXi ifContact;

	int edgeStart;
	int edgeEnd;

	int bendingStart;
	int bendingEnd;

	VectorXi connerIndex;

	MatrixXd boundaryStart;

	private:
};

#endif
