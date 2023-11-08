
#include <iostream>
#include <vector>

#include<eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues> 

#include "lqr.h"

using namespace Eigen;
using namespace std;


int main()
{
// Here we construct a test case. It is a mass-spring-damper system. 

//# masses, spring and damper constants
double m1=2  ; double m2=2   ; double k1=100  ; double k2=200 ; double d1=1  ; double d2=5; 
//# define the continuous-time system matrices and initial condition

    Matrix <double,4,4> Ac {{0, 1, 0, 0},
                            {-(k1+k2)/m1 ,  -(d1+d2)/m1 , k2/m1 , d2/m1},
                            {0 , 0 ,  0 , 1},
                            {k2/m2,  d2/m2, k2/m2, d2/m2}};
    Matrix <double,4,1> Bc {{0},{0},{0},{1/m2}}; 
    Matrix <double,1,4> Cc {{1,0,0,0}};
    
    // initial condition for simulation of the closed loop system
    Matrix <double,4,1> x0 {{5},{-3},{10},{-1}};

    // total number of simulation steps
    int simulationTimeSteps=50;
    // discretization time step 
    double h=0.1;
    // extract the number of rows and columns
    unsigned int n=Ac.rows();
    unsigned int m=Bc.cols();
    
    // construct the weigthing matrices
    // state weighting matrix
    MatrixXd weightMatrixQ;
    weightMatrixQ= 100*MatrixXd::Identity(n,n);

    // input weighting matrix
    MatrixXd weightMatrixR;
    weightMatrixR= 0.01*MatrixXd::Identity(m,m);
    
    // mixed state-input weighting matrix
    MatrixXd weightMatrixS;
    weightMatrixS.resize(n, m); 
    weightMatrixS.setZero();
    int maxIteration=5000;
    double toleranceConvergence=1e-8;
    
    // construct the LQR controller object
    LQRController lqr(Ac, Bc, weightMatrixQ, weightMatrixR,weightMatrixS);
    // compute the solution
    lqr.ComputeSolution(maxIteration,toleranceConvergence);
    // simulate the system
    lqr.SimulateSystem(x0,simulationTimeSteps,h);
    // save the data
    //void SaveData(string KFile,string AclFile,string solutionRiccatiFile,string simulatedStateTrajectoryFile) const;
    lqr.SaveData("computedK.csv","computedAcl.csv","computedSolutionRiccati.csv","computedSimulatedStateTrajectory.csv");

return 0;
}