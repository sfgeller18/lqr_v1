#ifndef LQRCONTROLLER_H
#define LQRCONTROLLER_H

#include<string>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;


class LQRController{
    public:
        LQRController();
        
        // Input arguments"
        // Ainput, Binput - A and B system matrices 
        // Qinput - state weighting matrix
        // Rinput - control input weighting matrix
        // Sinput - state-input mixed weighting matrix
                
        LQRController(MatrixXd Ainput, MatrixXd Binput, 
                       MatrixXd Qinput, MatrixXd Rinput,MatrixXd Sinput);
        
        /*
        NOTE: The initial guess of the solution should be a stabilizing matrix
        That is, the matrix A-B*K0 is a stable matrix, where K0 is the initial gain computed
        on the basis of the initial guess X0 of the solution of the Riccati equation.
        This function will compute X0.2
        */

        void ComputeInitialGuess(MatrixXd &initialGuessMatrix);

        void SolveLyapunovEquation(const MatrixXd &Am, const MatrixXd &  Rm, MatrixXd &SolutionMatrix, double &residual);
        
        void ComputeSolution(int numberIterations, double tolerance); //Newton Method Solution

        void SimulateSystem(MatrixXd x0,int simulationTimeSteps, double h);

        void SaveData(const string& KFile,const string& AclFile,const string& solutionRiccatiFile,const string& simulatedStateTrajectoryFile) const;

    private:

	    MatrixXd A,B; // system matrices
        MatrixXd Q,R,S;   // weighting matrices
	   
        MatrixXd K; // LQR control gain matrix 
        MatrixXd Acl; // LQR closed-loop matrix Acl=A-B*K
        MatrixXd solutionRiccati; // solution of the Riccati equation
        MatrixXd In; // identity matrix for computations
        MatrixXd simulatedStateTrajectory; // simulated state trajectory of the closed-loop state-space model

        unsigned int n; // state dimension 
        unsigned int m; // input dimension
};

#endif