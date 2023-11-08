
#include <iostream>
#include<tuple>
#include<string>
#include<fstream>
#include<vector>
#include<cerrno>
#include <cstring>
#include <unistd.h>
#include<eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <eigen3/Eigen/Eigenvalues> 

#include "lqr.h"


LQRController::LQRController(MatrixXd Ainput, MatrixXd Binput, 
                       MatrixXd Qinput, MatrixXd Rinput,MatrixXd Sinput)
{
    A=Ainput; B=Binput; Q=Qinput; R=Rinput; S=Sinput;
    n=A.rows(); m=B.cols();
    if (A.rows() != A.cols()) {
    throw std::invalid_argument("Matrix A must be square (n x n) for LQR control.");
}

if (B.rows() != n) {
    throw std::invalid_argument("Matrix B must have n rows for LQR control.");
}

if (Q.rows() != n || Q.cols() != n) {
    throw std::invalid_argument("Matrix Q must be (n x n) for LQR control.");
}

if (R.rows() != m || R.cols() != m) {
    throw std::invalid_argument("Matrix R must be (m x m) for LQR control.");
}

if (S.rows() != n || S.cols() != m) {
    throw std::invalid_argument("Matrix S must be (n x m) for LQR control.");
}

    K.resize(m,n); K.setZero();

    Acl.resize(n,n); Acl.setZero();

    solutionRiccati.resize(n,n);
    solutionRiccati.setZero();

    In= MatrixXd::Identity(n,n);
    
    
}

void LQRController::ComputeInitialGuess(MatrixXd &initialGuessMatrix)
{
    
    // here, we have to introduce new matrices in order to fit everything 
    // in the computational framework explained in the book
    MatrixXd Anew;
    MatrixXd Bnew;
    Anew=A-B*(R.inverse())*(S.transpose());
    Bnew=B*(R.inverse())*(B.transpose());
    
    // compute the eigenvalues of the open-loop matrix
    EigenSolver<MatrixXd> eigenValueSolver(Anew);
    // used to store the current eigenvalue
    complex<double> lambda; 
    // beta parameter from the book
    double bParam=0;
    // we store the real parts of the eigenvalues in this vector
    vector <double> realPartEigenValues;
    // extract the real parts of the eigenvalues and store them in the vector
    for (int i=0; i<n; i++)
    {
        lambda= eigenValueSolver.eigenvalues()[i];
        realPartEigenValues.push_back(lambda.real());
    }
    // compute the bParam
    bParam=(*min_element(realPartEigenValues.begin(), realPartEigenValues.end()));
    //cout<<bParam;
    bParam=(-1)*bParam;
    // be carefull about selecting this constant 0.02, you might need to change it if the convergence 
    // is not fast, 
    bParam=max(bParam,0.0)+0.02;
    // Abar is A+beta*I matrix from the equation 3.14 in the book
    MatrixXd Abar;
    // note that we need to transpose, since our Lyapunov solver is written like this
    // A^{T}X+XA=RHS
    // in the book, it is 
    // AX+XA^{T}=RHS 
    Abar=(bParam*In+Anew).transpose();
    // right-hand side of the equation 3.14 from the book
    MatrixXd RHS;
    RHS=2*Bnew*(Bnew.transpose());
    // solve the Lyapunov equation
    MatrixXd solutionLyapunov;
    solutionLyapunov.resize(n, n); 
    solutionLyapunov.setZero();
    // this residual is filled by the solution residual computed by the Lyapunov solver
    double residualLyap=10e10;
    // call the Lyapunov solver
    SolveLyapunovEquation(Abar, RHS, solutionLyapunov, residualLyap);
    // compute the initial guess matrix by using the equation from the book
    initialGuessMatrix=(Bnew.transpose())*(solutionLyapunov.inverse());
    // ensure that this solution is symmetric - this is a pure heuristic...
    initialGuessMatrix=0.5*(initialGuessMatrix+initialGuessMatrix.transpose());
    //cout<<endl<<"Initial guess of the solution:"<<endl;
    //cout<<endl<<initialGuessMatrix<<endl;
}


void LQRController::SolveLyapunovEquation(const MatrixXd &Am, const MatrixXd &  Rm, MatrixXd &SolutionMatrix, double &residual)
{
    MatrixXd LHSMatrix;
    LHSMatrix=kroneckerProduct(In,Am.transpose())+kroneckerProduct(Am.transpose(),In);
    MatrixXd solutionLyapunovVector;
    MatrixXd RHSVector;
    RHSVector=Rm.reshaped(n*n,1);

    //solutionLyapunovVector=LHSMatrix.colPivHouseholderQr().solve(RHSVector); (if known not to be singular)
     solutionLyapunovVector=LHSMatrix.completeOrthogonalDecomposition().solve(RHSVector);
    SolutionMatrix=solutionLyapunovVector.reshaped(n,n);

    MatrixXd residualMatrix;
    residualMatrix=(Am.transpose())*SolutionMatrix+SolutionMatrix*Am-Rm;
    residual=residualMatrix.squaredNorm();
}



void LQRController::ComputeSolution(int maxNumberIterations, double tolerance)
{
    MatrixXd initialGuess;
    initialGuess.resize(n, n);
    initialGuess.setZero();

    ComputeInitialGuess(initialGuess);

    MatrixXd initialK;
    MatrixXd initialAcl;
    initialK = (R.inverse()) * ((B.transpose()) * initialGuess + S.transpose());
    initialAcl = A - B * initialK;

    solutionRiccati = initialGuess;
    MatrixXd solutionLyapunov;
    solutionLyapunov.resize(n, n);
    solutionLyapunov.setZero();

    MatrixXd RHS;
    MatrixXd tmpMatrix;
    MatrixXd updateMatrix;
    double residualLyap = 10e10;
    double errorConvergence = 10e10;
    int currentIteration = 0;

    while (currentIteration <= maxNumberIterations && errorConvergence >= tolerance)
    {
        tmpMatrix = (B.transpose()) * solutionRiccati + S.transpose();
        K = (R.inverse()) * tmpMatrix;
        Acl = A - B * K;
        RHS = Q + (K.transpose()) * R * K - S * K - (K.transpose()) * (S.transpose());
        RHS = -RHS;
        SolveLyapunovEquation(Acl, RHS, solutionLyapunov, residualLyap);

        updateMatrix = solutionLyapunov - solutionRiccati;
        errorConvergence = (updateMatrix.cwiseAbs().colwise().sum().maxCoeff()) / (solutionRiccati.cwiseAbs().colwise().sum().maxCoeff());
        solutionRiccati = solutionLyapunov;
        currentIteration = currentIteration + 1;
    }

    if (errorConvergence < tolerance)
    {
        cout << endl << "Solution computed within prescribed error tolerances!" << endl;
        cout << "Converged error: " << errorConvergence << endl;
        cout << "Number of iterations:" << currentIteration << endl;
    }

    if (currentIteration > maxNumberIterations)
    {
        cout << endl << "Maximum number of maximum iterations exceeded!" << endl;
        cout << "However, the current error is not below the error tolerance." << endl;
        cout << "Consider increasing the maximum number of iterations or decreasing the tolerances." << endl;
        cout << "Current error:" << errorConvergence << endl;
    }

    EigenSolver<MatrixXd> eigenValueSolver2(Acl);
    cout << endl << "The eigenvalues of the final closed-loop matrix: " << endl << eigenValueSolver2.eigenvalues() << endl;
    cout << endl << "Computed K matrix: " << K << endl;
}

void LQRController::SimulateSystem(MatrixXd x0, int simulationTimeSteps, double h)
{
    simulatedStateTrajectory.resize(n, simulationTimeSteps);
    simulatedStateTrajectory.setZero();
    simulatedStateTrajectory.col(0) = x0;

    MatrixXd AclDiscrete;
    AclDiscrete = (In - h * Acl).inverse();
    simulatedStateTrajectory.col(0) = x0;

    for (int i = 0; i < simulationTimeSteps - 1; i++)
    {
        simulatedStateTrajectory.col(i + 1) = AclDiscrete * simulatedStateTrajectory.col(i);
    }

Eigen::VectorXd columnNumbers = Eigen::VectorXd::LinSpaced(simulatedStateTrajectory.cols(), 1, simulatedStateTrajectory.cols());

// Add the column numbers as the first row
Eigen::MatrixXd dataWithColumnNumbers(simulatedStateTrajectory.rows() + 1, simulatedStateTrajectory.cols());
dataWithColumnNumbers << columnNumbers.transpose(), simulatedStateTrajectory;
simulatedStateTrajectory = dataWithColumnNumbers;
}

void LQRController::SaveData(const string& KFile, const string& AclFile, const string& solutionRiccatiFile, const string& simulatedStateTrajectoryFile) const
{
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    ofstream file1("/home/sfgeller18/projects/lqr/output/" + KFile);
    if (file1.is_open())
    {
        file1 << K.format(CSVFormat);
        file1.close();
        cout << "Data saved to file: " << KFile << endl;
    }
    else
    {
        cerr << "Error opening file " << KFile << ": " << strerror(errno) << endl;
    }

    ofstream file2("/home/sfgeller18/projects/lqr/output/" + AclFile);
    if (file2.is_open())
    {
        file2 << Acl.format(CSVFormat);
        file2.close();
        cout << "Data saved to file: " << AclFile << endl;
    }
    else
    {
        cerr << "Error opening file " << AclFile << ": " << strerror(errno) << endl;
    }

    ofstream file3("/home/sfgeller18/projects/lqr/output/" + solutionRiccatiFile);
    if (file3.is_open())
    {
        file3 << solutionRiccati.format(CSVFormat);
        file3.close();
        cout << "Data saved to file: " << solutionRiccatiFile << endl;
    }
    else
    {
        cerr << "Error opening file " << solutionRiccatiFile << ": " << strerror(errno) << endl;
    }

    ofstream file4("/home/sfgeller18/projects/lqr/output/" + simulatedStateTrajectoryFile);
    if (file4.is_open())
    {
        file4 << (simulatedStateTrajectory.transpose()).format(CSVFormat);
        file4.close();
        cout << "Data saved to file: " << simulatedStateTrajectoryFile << endl;
    }
    else
    {
        cerr << "Error opening file " << simulatedStateTrajectoryFile << ": " << strerror(errno) << endl;
    }
}


