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

	Vector3d t1;
	Vector3d d1;
	Vector3d d2;

	double theta;

	VectorXi arrayNum;
};

struct bendingElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	int ne_1;
	int ne_2;

	VectorXi arrayNum;

	Vector3d x_1;
	Vector3d x_2;
	Vector3d x_3;

	Vector3d e_1;
	Vector3d e_2;

	double theta_1;
	double theta_2;

	double norm_1;
	double norm_2;

	Vector3d t_1;
	Vector3d d_11;
	Vector3d d_12;
	Vector3d t_1_old;
	Vector3d d_11_old;
	Vector3d d_12_old;
	Vector3d m_11;
	Vector3d m_12;

	Vector3d t_2;
	Vector3d d_21;
	Vector3d d_22;
	Vector3d t_2_old;
	Vector3d d_21_old;
	Vector3d d_22_old;
	Vector3d m_21;
	Vector3d m_22;

	double voroniLength;

	Vector3d kb;
	Vector2d kappa;
	Vector2d kappaBar;

	double refTwist;
	double refTwist_old;

	int sign_1;
	int sign_2;
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
    std::vector<Vector2i> bending;
    std::vector<int> constraint;

	std::vector<edgeElement> v_edgeElement;
	std::vector<bendingElement> v_bendingElement;
	
	int temp;

	double crossSectionalArea;

	int nv;
	int ne;
	int edgeNum;
	int bendingNum;

	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector3d position, int k);
	void setOneVertexBoundaryCondition(double position, int i, int k);
	void setThetaBoundaryCondition(double position, int k);
	void setConstraint(double position, int k);

	void computeEdge();
	void computeBending();

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
	double GJ;

	void updateEdgePair();
	void updateBendingPair();

	void parallelTansport(const Vector3d &d1_1,const Vector3d &t1, const Vector3d &t2, Vector3d &d1_2);
	double signedAngle(const Vector3d &u, const Vector3d &v, const Vector3d &n);
	void rotateAxisAngle(Vector3d &v,const Vector3d &z,const double &theta);

	private:
};

#endif
