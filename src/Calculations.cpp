/**
 *  \file Calculations.cpp
 *
 *  \brief  Contains non-member calculation functions
 *
 *   This library contains functions to perform miscellaneous calculations.
 *   Common calculations, such as EOM evaluation or primary location, should
 *   be included in the tpat_model class and its derivatives.
 *   
 *   \author Andrew Cox
 *  \version May 25, 2016
 *  \copyright GNU GPL v3.0
 */
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Astrohelion.
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Calculations.hpp"

#include "AsciiOutput.hpp"
#include "BaseArcset.hpp"
#include "BodyData.hpp"
#include "Common.hpp"
#include "CorrectionEngine.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "Nodeset_bc4bp.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SimEngine.hpp"
#include "SysData_2bp.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Traj.hpp"
#include "Traj_2bp.hpp"
#include "Traj_bc4bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Utilities.hpp"

#include <cspice/SpiceUsr.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstring>
#include <string>
#include <typeinfo>


namespace astrohelion{

//-----------------------------------------------------
//      General Utility Functions
//-----------------------------------------------------
/**
 * @addtogroup util
 * \{
 */


/**
 *  \brief convert a date to ephemeris time
 *  
 *  \param pDate a string representing the date. The string can be formatted in one
 *  of two ways. First, a Gregorian-style date: 'YYYY/MM/DD HH:II:SS' (UTC, 24-hour clock);
 *  The time 'HH:II:SS' can be ommited, and the function will assume the time is 0:0:00
 *  Second, a Julian date (UTC) can be input with the format 'jd #' where '#' represents
 *  the Julian date.
 *
 *  @return the J2000 ephemeris time, or number of seconds after Jan 1, 2000 at 0:0:00 UTC.
 *  \throws Exception if the SPICE kernels cannot be loaded: the kernel names and
 *  filepaths are located in the settings XML file
 */
double dateToEphemerisTime(const char *pDate){
    // Convert the date to ephemeris time
    double et = 0;
    str2et_c(pDate, &et);

    return et;
}//==========================================

/**
 *  \brief Convert a gregorian date to Julian Date
 * 
 *  \param yr Year
 *  \param mo Month (Jan = 1, ..., Dec = 12)
 *  \param d Day of mongth
 *  \param h Hour (24-hr clock)
 *  \param m minute
 *  \param s second
 *  @return Julian Date (days)
 */
double gregorianToJulian(double yr, double mo, double d, double h, double m, double s){
    return 367.0*yr - floor(7.0*(yr + floor((mo + 9.0)/12.0))/4.0) +
        floor(275.0*mo/9.0) + d + 1721013.5 + ((5.0/60.0 + m)/60.0 + h + s/3600.0)/24.0;
}//==========================================

/**
 *  \brief Determine the Greenwich Sidereal Time (i.e., angle) at the specified date
 *  @details Input date must be in UT1 time (always within +/- 0.9 seconds of UTC thanks
 *  to leapseconds)
 *  
 *  \param yr Year
 *  \param mo Month (Jan = 1, ..., Dec = 12)
 *  \param d Day of mongth
 *  \param h Hour (24-hr clock)
 *  \param m minute
 *  \param s second
 *  @return The Greenwich Sidereal Time (GST) in radians at the specified date
 */
double dateToGST(double yr, double mo, double d, double h, double m, double s){

    // Julian date at midnight
    double JD_midnight = gregorianToJulian(yr, mo, d, 0, 0, 0);

    double time_ut1 = (JD_midnight - 2451545.0)/36525.0; // number of Julian Centuries ellapsed since J2000

    // GST angle at midnight, degrees
    double GST_T0 = 100.4606184 + 36000.77005361*time_ut1 + 0.00038793*time_ut1*time_ut1 - (2.6e-8)*pow(time_ut1,3);
    
    // Add rotation past midnight using Earth's average rotation rate
    GST_T0 += 0.250684477337*(h*60.0 + m + s/60.0);

    // Scale back to a reasonable angle value
    GST_T0 -= std::floor(GST_T0/360.0)*360.0;
    
    // GST angle at midnight + earth rotation (deg/min) times minutes past midnight
    return GST_T0*PI/180.0;
}//==========================================

/**
 *  \brief use least squares to predict new values of variables in a continuation process
 *
 *  This function uses a 2nd-order polynomial fit to predict a set of
 *  dependent variables using past relationships between an independent
 *  variable and the dependent variables. The number of points considered
 *  in those past relationships can be adjusted to vary the "stiffness"
 *  of the fit.
 *
 *  Ocasionally the 2nd-order fit does not work well and we encounter a 
 *  singular (or very near singular) matrix. In this case, the algorithm
 *  will apply linear regression, which generally solves the problem and
 *  results in a safe inversion.
 *
 *  For all these inputs, the "state vector" must take the following form:
 *      [x, y, z, xdot, ydot, zdot, Period, Jacobi Constant]
 *
 *  \param indVarIx the index of the state variable to use as the indpendent variable
 *  \param nextInd The next value of the independent variable
 *  \param depVars a vector specifying the indices of the states that will be
 *  dependent variables. The algorithm will predict fugure values for these
 *  variables based on how they have changed with the independent variable.
 *  \param varHistory a vector representing an n x 8 matrix which contains
 *  information about previous states, period, and JC. n should be at least 
 *  3. If it is larger than MemSize, only the first set of MemSize rows will
 *  be used.
 *
 *  @return an 8-element vector with predictions for the dependent variables.
 *  If a particular variable has not been predicted, its place will be kept with
 *  a NAN value.
 *
 *  An example may make things more clear:
 *
 *  Say I am continuing a family and am using x as the natural parameter in the
 *  continuation. indVarIx would be 0 to represent x. We input the value of x
 *  for the next orbit in the family (nextInd) and specify which variables (from
 *  the 8-element "state") we would like to have predicted by least-squares.
 *  \throws Exception if the <tt>varHistory</tt> vector contains fewer than 
 *  8 elements or if the <tt>depVars</tt> vector has no data
 */
std::vector<double> familyCont_LS(int indVarIx, double nextInd, std::vector<int> depVars, std::vector<double> varHistory){
    const int STATE_SIZE = 8;
    const double EPS = 1e-14;

    if(varHistory.size() < STATE_SIZE)
        throw Exception("Calculations::familyCont_LS: Not enough data to create A matrix\n");

    if(depVars.size() == 0)
        throw Exception("Calculations::familyCont_LS: Not enough data to create B matrix\n");

    // Form A and B matrices
    std::vector<double> A_data;
    std::vector<double> B_data;

    for(unsigned int n = 0; n < varHistory.size()/STATE_SIZE; n++){
        double d = varHistory[n*STATE_SIZE + indVarIx];
        A_data.push_back(d*d);  // ind. var^2
        A_data.push_back(d);    // ind. var^1
        A_data.push_back(1);    // ind. var^0

        for(unsigned int p = 0; p < depVars.size(); p++){
            // vector of dependent variables
            B_data.push_back(varHistory[n*STATE_SIZE + depVars[p]]);
        }
    }

    MatrixXRd A = Eigen::Map<MatrixXRd>(&(A_data[0]), varHistory.size()/STATE_SIZE, 3);
    MatrixXRd B = Eigen::Map<MatrixXRd>(&(B_data[0]), varHistory.size()/STATE_SIZE, depVars.size());

    // Generate coefficient matrix; these are coefficients for second-order
    // polynomials in the new independent variable
    MatrixXRd G(A.cols(), A.cols());
    G.noalias() = A.transpose()*A;
    
    Eigen::JacobiSVD<MatrixXRd> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(1e-14);
    Eigen::VectorXd S = svd.singularValues();
    double smallestVal = S.minCoeff();
    
    Eigen::RowVectorXd P(depVars.size());

    if(smallestVal > EPS){
        // Use 2nd-order polynomial fit; solve GC = A.transpose()*B for C
        MatrixXRd C = G.fullPivLu().solve(A.transpose()*B);
        
        double indMatData[] = {nextInd*nextInd, nextInd, 1};
        Eigen::RowVector3d indMat = Eigen::Map<Eigen::RowVector3d>(indMatData, 1, 3);

        P.noalias() = indMat*C;
    }else{
        // User 1st-order polynomial fit
        std::vector<double> A_lin_data;
        for(unsigned int n = 0; n < varHistory.size()/STATE_SIZE; n++){
            A_lin_data.push_back(varHistory[n*STATE_SIZE + indVarIx]);
            A_lin_data.push_back(1);
        }

        MatrixXRd A_lin = Eigen::Map<MatrixXRd>(&(A_lin_data[0]), varHistory.size()/STATE_SIZE, 2);
        Eigen::Matrix2d G;
        G.noalias() = A_lin.transpose()*A_lin;

        MatrixXRd C = G.fullPivLu().solve(A_lin.transpose()*B);

        Eigen::RowVector2d indMat(nextInd, 1);

        P.noalias() = indMat*C;
    }

    // Insert NAN for states that have not been predicted
    std::vector<double> predicted;
    predicted.assign(STATE_SIZE, NAN);
    for(unsigned int i = 0; i < depVars.size(); i++)
        predicted[depVars[i]] = P(i);

    return predicted;
}//====================================================

/**
 *  \brief construct a matrix to mirror a 6-d state over the specified plane or axis
 *  \param mirrorType describes how to mirror a 6-d state
 *  @return a 6x6 matrix that will mirror a 6-d state over the specified plane or axis
 */
MatrixXRd getMirrorMat(Mirror_tp mirrorType){
    switch(mirrorType){
        case Mirror_tp::MIRROR_XZ:
        {
            double data[] = {1, 0, 0, 0, 0, 0,
                            0, -1, 0, 0, 0, 0,
                            0, 0, 1, 0, 0, 0,
                            0, 0, 0, -1, 0, 0,
                            0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, -1};
            return Eigen::Map<MatrixXRd>(data, 6, 6);
        }
        case Mirror_tp::MIRROR_X_AX_H:
        case Mirror_tp::MIRROR_X_AX_V:
        {
            double data[] = {1, 0, 0, 0, 0, 0,
                            0, -1, 0, 0, 0, 0,
                            0, 0, -1, 0, 0, 0,
                            0, 0, 0, -1, 0, 0,
                            0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, 1};
            return Eigen::Map<MatrixXRd>(data, 6, 6);
        }
        default:
            astrohelion::printErr("Calculations::getMirrorMat: Mirror type is not implemented; returning identiy\n");
            return Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Identity();
    }
}//====================================================

/**
 *  \brief Sort eigenvalues
 *  \details This algorithm leverages three simple axioms to determine if the eigenvalues
 *  are in a consistent order.
 *  
 *      1. Pairs of eigenvalues should occur in sequence ( [0,1] or [2,3], etc.) and have
 *      a product equal to one
 *      
 *      2. Eigenvalues evolve continuously along the family, thus the changes between subsequent 
 *      eigenvalues should be small
 *      
 *      3. Eigenvectors evolve continuously along the family, thus the dot product between subsequent
 *      eigenvectors should be small
 * 
 *  Axiom #1 is not always applicable to trajectories but is always applicable to families of periodic
 *  orbits in the CR3BP
 *  
 *  \param eigVals Vector of all eigenvalues along a family of trajectories. It is assumed that there
 *  are six eigenvalues per trajectory.
 *  
 *  \param eigVecs Vector of 6x6 eigenvector matrices; the eigenvectors are columns of the matrix and their
 *  order is consistent with the order of the six eigenvalues associated with the same trajectory
 * 
 *  \return Rearranged indices that describe the correct order of the eigenvalues/vectors
 */
std::vector<unsigned int> sortEig(std::vector<cdouble> eigVals, std::vector<MatrixXRcd> eigVecs){
    std::vector<unsigned int> sortedIxs(eigVals.size(), 0);

    if(eigVals.size() == 0){
        astrohelion::printErr("Calculations::sortEig: No eigenvalues - easy to sort! :P\n");
        return sortedIxs;
    }

    if(eigVals.size() % 6 != 0){
        astrohelion::printErr("Calculations::sortEig: Must have 6n eigenvalues!\n");
        return sortedIxs;
    }

    // Generate all permutations of the indices 0 through 5
    std::vector<unsigned int> vals {0,1,2,3,4,5};
    std::vector<unsigned int> ixPerms = astrohelion::generatePerms<unsigned int>(vals);
    std::vector<float> cost(ixPerms.size()/6, 0);

    // Sort the first set of eigenvalues so that reciprocal pairs occur near one another
    unsigned int s = 0;
    for(unsigned int p = 0; p < ixPerms.size()/6; p++){
        for(unsigned int i = 0; i < 3; i++){
            // Penalize configures with pairs that are not reciprocal
            cost[p] += std::abs(1.0 - eigVals[s*6 + ixPerms[p*6 + 2*i]] * eigVals[s*6 + ixPerms[p*6 + 2*i + 1]]);
        }
    }

    // Find the minimum cost
    std::vector<float>::iterator minCostIt = std::min_element(cost.begin(), cost.end());
    unsigned int minCostIx = minCostIt - cost.begin();

    // printf("Minimum cost for member %u is %f on permutation %u\n", s, *minCostIt, minCostIx);

    std::vector<cdouble> prevSortVal(6,0);
    MatrixXRcd prevSortVec = MatrixXRcd::Zero(6,6);
    for(unsigned int i = 0; i < 6; i++){
        sortedIxs[s*6 + i] = ixPerms[minCostIx*6 + i];
        prevSortVal[i] = eigVals[s*6 + ixPerms[minCostIx*6 + i]];
        prevSortVec.col(i) = eigVecs[s].col(ixPerms[minCostIx*6+i]);
    }

    // printf("[%f%+fi, %f%+fi, %f%+fi, %f%+fi, %f%+fi, %f%+fi]\n",
    //     std::real(eigVals[sortedIxs[s*6+0]]), std::imag(eigVals[sortedIxs[s*6+0]]),
    //     std::real(eigVals[sortedIxs[s*6+1]]), std::imag(eigVals[sortedIxs[s*6+1]]),
    //     std::real(eigVals[sortedIxs[s*6+2]]), std::imag(eigVals[sortedIxs[s*6+2]]),
    //     std::real(eigVals[sortedIxs[s*6+3]]), std::imag(eigVals[sortedIxs[s*6+3]]),
    //     std::real(eigVals[sortedIxs[s*6+4]]), std::imag(eigVals[sortedIxs[s*6+4]]),
    //     std::real(eigVals[sortedIxs[s*6+5]]), std::imag(eigVals[sortedIxs[s*6+5]]));

    for(s = 1; s < eigVals.size()/6; s++){
        MatrixXRd dp_err = MatrixXRd::Zero(6,6);
        for(unsigned int i = 0; i < 6; i++){
            for(unsigned int j = i; j < 6; j++){
                dp_err(i,j) = std::abs(1.0 - std::abs(prevSortVec.col(i).dot(eigVecs[s].col(j))));
                dp_err(j,i) = dp_err(i,j);
            }
        }

        cost.clear();
        cost.assign(ixPerms.size()/6, 0);
        for(unsigned int p = 0; p < ixPerms.size()/6; p++){
            for(unsigned int i = 0; i < 6; i++){
                // Assign cost based on the dot product between this new arrangement
                // of eigenvectors and the previous arrangement
                cost[p] += dp_err(i, ixPerms[6*p + i]);

                // Assign cost based on the distance from the previous sorted eigenvalues
                cost[p] += std::abs(eigVals[6*s + ixPerms[6*p + i]] - prevSortVal[i]);

                // Assign cost based on the reciprocal nature of the eigenvalues
                // This test could be removed to sort eigenvalues that don't occur in pairs
                if(i%2 == 1)
                    cost[p] += std::abs(1.0 - eigVals[s*6 + ixPerms[p*6 + i-1]] * eigVals[s*6 + ixPerms[p*6 + i]]);
            }
        }

        // Find the minimum cost
        minCostIt = std::min_element(cost.begin(), cost.end());
        minCostIx = minCostIt - cost.begin();

        // printf("Minimum cost for member %u is %f on permutation %u\n", s, *minCostIt, minCostIx);
        for(unsigned int i = 0; i < 6; i++){
            sortedIxs[s*6 + i] = ixPerms[minCostIx*6 + i];
            prevSortVal[i] = eigVals[s*6 + sortedIxs[s*6 + i]];
            prevSortVec.col(i) = eigVecs[s].col(sortedIxs[s*6 + i]);
        }
    }

    return sortedIxs;
}//====================================================

/**
 *  \brief Compute manifolds of a CR3BP trajectory
 *  @details [long description]
 * 
 *  \param type The type of manifolds to generate
 *  \param pPerOrbit A periodic, CR3BP orbit. No checks are made to ensure periodicity,
 *  so this function also performs well for nearly periodic segments of quasi-periodic
 *  arcs. If an arc that is not approximately periodic is input, the behavior may be... strange.
 *  \param numMans The number of manifolds to generate
 *  \param tof Time-of-flight along each manifold arc after stepping off the specified orbit
 *  @return a vector of trajectory objects, one for each manifold arc.
 *  \throws Exception if the eigenvalues cannot be computed, or if only one
 *  stable or unstable eigenvalue is computed (likely because of an impropper monodromy matrix)
 */
std::vector<Traj_cr3bp> getManifolds(Manifold_tp type, const Traj_cr3bp *pPerOrbit, int numMans, double tof){
    // Get eigenvalues of monodromy matrix
    MatrixXRd mono = pPerOrbit->getSTMByIx(-1);

    Eigen::EigenSolver<MatrixXRd> eigensolver(mono);
    if(eigensolver.info() != Eigen::Success)
        throw Exception("Calculations::getManifolds: Could not compute eigenvalues of monodromy matrix");

    Eigen::VectorXcd vals = eigensolver.eigenvalues();
    MatrixXRcd eigVecs = eigensolver.eigenvectors();
    std::vector<cdouble> eigData(vals.data(), vals.data()+6);

    // Sort eigenvalues to put them in a "propper" order, get indices
    // to sort the eigenvectors to match
    std::vector<MatrixXRcd> vecvec {eigVecs};
    std::vector<unsigned int> sortedIx = sortEig(eigData, vecvec);
    std::vector<cdouble> sortedEig;
    for(unsigned int i = 0; i < 6; i++){
        sortedEig.push_back(eigData[sortedIx[i]]);
    }

    // Figure out which eigenvalues are the stable and unstable ones,
    // delete the rest
    std::vector<cdouble> nonCenterVals;
    std::vector<cdouble> nonCenterVecs;
    for(unsigned int c = 0; c < sortedEig.size(); c++){
        double realErr = std::real(sortedEig[c]) - 1.0;
        double imagErr = std::imag(sortedEig[c]);

        if(std::abs(realErr) > 1e-5 && std::abs(imagErr) < 1e-5){
            // Keep this eigenvalue/eigenvector pair
            nonCenterVals.push_back(sortedEig[c]);
            unsigned int vecIx = sortedIx[c];
            nonCenterVecs.insert(nonCenterVecs.end(),
                eigVecs.data()+vecIx*6, eigVecs.data()+(vecIx+1)*6); 
        }
    }

    std::vector<Traj_cr3bp> allManifolds;

    if(nonCenterVals.size() == 0){
        astrohelion::printWarn("Calculations::getManifolds: No stable/unstable eigenvalues were found\n");
        return allManifolds;
    }

    if(nonCenterVals.size() == 1){
        throw Exception("Calculations::getManifolds: Only found one stable/unstable eigenvalue");
    }

    if(nonCenterVals.size() > 2){
        astrohelion::printWarn("Calculations::getManifolds: Stable/Unstable subspace is larger than 2D. Only the first pair will be considered\n");
    }

    /** TODO: If and when I can flexibily integrate generic trajectories and nodesets,
     *  I would like to replace this with a nodeset to space the points evenly
     *  in time and/or arclength
     */
    // Get a bunch of points to use as starting guesses for the manifolds
    if(numMans > pPerOrbit->getNumNodes()){
        astrohelion::printWarn("Calculations::getManifolds: Requested too many manifolds... will return fewer\n");
        numMans = pPerOrbit->getNumNodes();
    }

    double stepSize = (static_cast<double>(pPerOrbit->getNumSegs()))/(static_cast<double>(numMans));
    std::vector<int> pointIx(numMans, 0);
    for(int i = 0; i < numMans; i++){
        pointIx[i] = floor(i*stepSize+0.5);
    }

    std::vector<double> realVecs = astrohelion::real(nonCenterVecs);
    MatrixXRd vecs = Eigen::Map<MatrixXRd>(&(realVecs[0]), nonCenterVecs.size()/6, 6);
    vecs.transposeInPlace();    // Transpose so eigenvectors are columns

    // NOW, copute the manifolds!
    SimEngine sim;
    double stepDist = 200;
    double charL = pPerOrbit->getSysData()->getCharL();
    for(int n = 0; n < numMans; n++){
        // Transform the eigenvectors to this updated time
        MatrixXRd newVecs = pPerOrbit->getSTMByIx(pointIx[n])*vecs;

        // Pick the direction from one of the transformed eigenvectors
        Eigen::VectorXd direction(6);
        for(unsigned int v = 0; v < 2; v++){
            Eigen::VectorXd eigVec = newVecs.col(v);
            double mag = sqrt(eigVec(0)*eigVec(0) + 
                eigVec(1)*eigVec(1) + eigVec(2)*eigVec(2));
            if(type == Manifold_tp::MAN_U_P || type == Manifold_tp::MAN_U_M){
                if(std::abs(nonCenterVals[v]) > 1){
                    direction = eigVec/mag;
                    sim.setRevTime(false);
                    break;
                }
            }else{
                if(std::abs(nonCenterVals[v]) < 1){
                    direction = eigVec/mag;
                    sim.setRevTime(true);
                    break;
                }
            }
        }

        // Make sure it is pointing in +x direction
        direction *= astrohelion::sign(direction(0));
        
        // Orient according to specified type
        if(type == Manifold_tp::MAN_U_M || type == Manifold_tp::MAN_S_M)
            direction *= -1;

        // Step away from the point on the arc in the direction of the eigenvector
        std::vector<double> state = pPerOrbit->getStateByIx(pointIx[n]);
        Eigen::VectorXd q0 = Eigen::Map<Eigen::VectorXd>(&(state[0]), 6, 1);
        q0 += stepDist/charL * direction;

        // Simulate for some time to generate a manifold arc
        Traj_cr3bp traj(static_cast<const SysData_cr3bp *>(pPerOrbit->getSysData()));
        sim.runSim(q0.data(), tof, &traj);
        allManifolds.push_back(traj);
    }

    return allManifolds;
}//====================================================

/**
 *  \brief Compute the stability index of a periodic orbit from a set of eigenvalues
 *  @details This algorithm assumes the orbit is periodic and that the eigenvalues 
 *  have been sorted (so they come in pairs).
 * 
 *  \param eigs A 6-element vector of eigenvalues associated with a periodic orbit
 *  @return the stability index, or NAN if no real, reciprocal eigenvalue pair is found
 *  \throws Exception if <tt>eigs</tt> does not have six elements
 */
double getStabilityIndex(std::vector<cdouble> eigs){
    if(eigs.size() != 6)
        throw Exception("Calculations::getStabilityIndex: Must input 6 eigenvalues!");

    double okErr = 1e-3;
    cdouble one(1,0);

    std::vector<EigValSet_tp> setTypes;
    setTypes.reserve(3);

    for(int set = 0; set < 3; set++){
        double sumImag = (std::abs(std::imag(eigs[set*2])) + std::abs(std::imag(eigs[set*2+1])))/2;
        double sumDistFromOne = (std::abs(eigs[set*2] - one) + std::abs(eigs[set*2+1] - one))/2;

        if(sumImag > okErr){
            setTypes[set] = EigValSet_tp::EIGSET_COMP_CONJ;
        }else{
            if(sumDistFromOne < okErr){
                setTypes[set] = EigValSet_tp::EIGSET_ONES;
            }else{
                setTypes[set] = EigValSet_tp::EIGSET_REAL_RECIP;
                return 0.5*std::real(eigs[set*2] + eigs[set*2 + 1]);
            }
        }
    }
    return NAN;
}//====================================================

/**
 *  \brief Compute the total delta-V along a corrected nodeset
 * 
 *  \param pIt pointer to an MultShootData object associated with a corrections process
 *  @return the total delta-V, units consistent with the nodeset's stored velocity states
 */
double getTotalDV(const MultShootData *pIt){
    double total = 0;
    for(unsigned int n = 0; n < pIt->deltaVs.size()/3; n++){
        total += sqrt(pIt->deltaVs[3*n + 0]*pIt->deltaVs[3*n + 0] +
            pIt->deltaVs[3*n + 1]*pIt->deltaVs[3*n + 1] + 
            pIt->deltaVs[3*n + 2]*pIt->deltaVs[3*n + 2]);
    }
    return total;
}//=====================================================

/**
 *  \brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  @details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  \param pNodeset A nodeset with some constraints
 */
bool finiteDiff_checkMultShoot(const Nodeset *pNodeset, Verbosity_tp verbosity, bool writeToFile){
    CorrectionEngine engine;  // Create engine with default settings
    return finiteDiff_checkMultShoot(pNodeset, engine, verbosity, writeToFile);
}//====================================================

/**
 *  \brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  @details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  \param pNodeset A nodeset with some constraints
 *  \param engine correction engine object configured with the appropriate settings (i.e.,
 *  equal arc time, etc.). Note that the maxIts, verbosity, and ignoreDiverge
 *  attributes of the engine will be overridden by this function.
 */
bool finiteDiff_checkMultShoot(const Nodeset *pNodeset, CorrectionEngine engine, Verbosity_tp verbosity, bool writeToFile){
    printVerb(verbosity >= Verbosity_tp::SOME_MSG, "Finite Diff: Checking DF matrix... ");
    // Create multiple shooter that will only do 1 iteration
    CorrectionEngine corrector(engine);
    corrector.setMaxIts(1);
    corrector.setVerbosity(Verbosity_tp::NO_MSG);
    corrector.setIgnoreDiverge(true);

    // Run multiple shooter to get X, FX, and DF
    MultShootData it = corrector.multShoot(pNodeset, NULL);
    Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(it.FX[0]), it.totalCons, 1);
    MatrixXRd DF = Eigen::Map<MatrixXRd>(&(it.DF[0]), it.totalCons, it.totalFree);
    MatrixXRd DFest = MatrixXRd::Zero(it.totalCons, it.totalFree);

    double pertSize = 1e-8;
    #pragma omp parallel for firstprivate(it, corrector)
    for(int i = 0; i < it.totalFree; i++){
        std::vector<double> pertX = it.X0;      // Copy unperturbed state vetor
        pertX[i] += pertSize;                   // add perturbation
        it.X = pertX;                           // Copy into iteration data
        MultShootData pertIt = corrector.multShoot(it, NULL);     // Correct perturbed state
        Eigen::VectorXd FX_up = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);

        // Do another process for opposite direction
        pertX = it.X0;
        pertX[i] -= pertSize;
        it.X = pertX;
        pertIt = corrector.multShoot(it, NULL);
        Eigen::VectorXd FX_down = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);

        // An iteration for twice the perturbation up
        pertX = it.X0;
        pertX[i] += 2*pertSize;
        it.X = pertX;
        pertIt = corrector.multShoot(it, NULL);
        Eigen::VectorXd FX_2up = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);

        // An iteration for twice the perturbation down
        pertX = it.X0;
        pertX[i] -= 2*pertSize;
        it.X = pertX;
        pertIt = corrector.multShoot(it, NULL);
        Eigen::VectorXd FX_2down = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);


        // Compute central difference
        Eigen::VectorXd col = (-1*FX_2up + 8*FX_up - 8*FX_down + FX_2down)/std::abs(12*pertSize);   // Five-point stencil
        // Eigen::VectorXd col = (FX_up - FX_down)/std::abs(2*pertSize);   // Central Difference
        // Eigen::VectorXd col = (FX_up - FX)/std::abs(pertSize);       // Forward difference
        DFest.block(0, i, it.totalCons, 1) = col;
    }


    MatrixXRd diff = DF - DFest;
    MatrixXRd DF_abs = DF.cwiseAbs();       // Get coefficient-wise absolute value
    MatrixXRd DFest_abs = DFest.cwiseAbs();

    diff = diff.cwiseAbs();                     // Get coefficient-wise aboslute value

    // // Divide each element by the magnitude of the DF element to get a relative difference magnitude
    // for(int r = 0; r < diff.rows(); r++){
    //     for(int c = 0; c < diff.cols(); c++){
    //         // If one of the elements is zero, let the difference just be the difference; no division
    //         if(DF_abs(r,c) > 1e-13 && DFest_abs(r,c) > 1e-13)   // consider 1e-13 to be zero
    //             diff(r,c) = diff(r,c)/DF_abs(r,c);

    //         // if(r == 50 && c == 98){
    //         //     printf("DF(50, 98) = %.4e\n", DF(r,c));
    //         //     printf("DFest(50, 98) = %.4e\n", DFest(r,c));
    //         // }
    //     }
    // }
    
    if(writeToFile){
        astrohelion::toCSV(DF, "FiniteDiff_DF.csv");
        astrohelion::toCSV(DFest, "FiniteDiff_DFest.csv");
        astrohelion::toCSV(diff, "FiniteDiff_Diff.csv"); 
    }
    // 

    Eigen::VectorXd rowMax = diff.rowwise().maxCoeff();
    Eigen::RowVectorXd colMax = diff.colwise().maxCoeff();

    // Make a map that links the row number (column in the Jacobian matrix)
    // of the free variable to the MSVarMap_Key that represents the MSVarMap_Obj 
    // that contains information about the free variable
    std::map<int, MSVarMap_Key> freeVarMap_rowNumToKey;
    for(auto& obj : it.freeVarMap){
        for(int r = 0; r < obj.second.nRows; r++){
            freeVarMap_rowNumToKey[obj.second.row0 + r] = obj.first;
        }
    }

    double rowMaxMax = rowMax.maxCoeff();
    double colMaxMax = colMax.maxCoeff();
    int errScalar = 10000;

    if(rowMaxMax < errScalar*pertSize && colMaxMax < errScalar*colMaxMax){
        if(verbosity >= Verbosity_tp::SOME_MSG)
            astrohelion::printColor(BOLDGREEN, "No significant errors!\n");

        return true;
    }else{
        if(verbosity >= Verbosity_tp::SOME_MSG){
            astrohelion::printColor(BOLDRED, "Significant errors!\n");
            printf("Maximum relative difference between computed DF and estimated DF\n");

            int conCount = 0;
            for(long r = 0; r < rowMax.size(); r++){
                if(r == 0 && it.totalCons > 0){
                    printf("Applies to %s %d: %s Constraint:\n", 
                        Constraint::getAppTypeStr(it.allCons[conCount].getAppType()),
                        it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
                }else if(conCount+1 < static_cast<int>(it.allCons.size()) && r >= it.conRows[conCount+1]){
                    conCount++;
                    printf("Applies to %s %d: %s Constraint:\n", 
                        Constraint::getAppTypeStr(it.allCons[conCount].getAppType()),
                        it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
                }
                astrohelion::printColor(rowMax[r] > errScalar*pertSize || std::isnan(rowMax[r]) ? RED : GREEN,
                    "  row %03zu: %.6e\n", r, rowMax[r]);
            }
            for(long c = 0; c < colMax.size(); c++){
                std::string parent = "unknown";
                std::string type = "unknown";
                MSVarMap_Key key;
                try{
                    key = freeVarMap_rowNumToKey.at(c);
                    MSVarMap_Obj obj = it.freeVarMap.at(key);
                    parent = MSVarMap_Obj::parent2str(obj.parent);
                    type = MSVarMap_Key::type2str(key.type);
                }catch(std::out_of_range &e){}

                astrohelion::printColor(colMax[c] > errScalar*pertSize || std::isnan(colMax[c]) ? RED : GREEN,
                    "Col %03zu: %s (%d)-owned %s: %.6e\n", c, parent.c_str(), key.id, type.c_str(), colMax[c]);
            }
        }

        return false;
    }
}//====================================================

/**
 *  \brief Use numerical integration to find a node on a trajectory at
 *  the specified time
 *  @details 
 * 
 *  \param traj A trajectory 
 *  \param t The time at which the node is located
 * 
 *  @return A node on the specified trajectory at time <tt>t</tt>
 */
Node interpPointAtTime(const Traj *traj, double t){
    // Find the time that is closest to t and before it (assuming allTimes progresses smoothly forwards in time)
    std::vector<double> allTimes = traj->getEpochs();
    size_t i = 0;
    while(i < allTimes.size() && allTimes[i] < t)
        i += 1;
    
    if(i > 0)
        i -= 1; // Rewind to the point before for forward integration
    
    // printf("Beginning at node %zu:\n", i);
    Node node = traj->getNodeByIx(i);
    // node.print();

    Traj *temp = new Traj(traj->getSysData());
    SimEngine sim;
    // sim.setVerbosity(Verbosity_tp::ALL_MSG);
    sim.setNumSteps(2);
    sim.setVarStepSize(false);
    sim.setRevTime(t - node.getEpoch() < 0);
    sim.runSim(node.getState(), node.getEpoch(), t - node.getEpoch(), temp);

    Node node2 = temp->getNodeByIx(-1);
    delete(temp);
    return node2;
}//====================================================

//-----------------------------------------------------
//      Orbit Determination Utility Functions
//-----------------------------------------------------


/**
 *  \brief Determine a set of spherical coordinates the describe 
 *  the position vector specified by x, y, and z
 *  @details [long description]
 * 
 *  \param x x-coordinate
 *  \param y y-coordinate
 *  \param z z-coordinate
 *  
 *  @return a vector containing {lat, long, R} in radians
 *  and the distance units of x, y, and z. Latitude takes a 
 *  value between -pi/2 and pi/2 while longitude takes a value
 *  between -pi and pi
 */
std::vector<double> getSpherical(double x, double y, double z){
    double R = sqrt(x*x + y*y + z*z);   // Distance from center of body to point
    double unit[] = {x/R, y/R, z/R};      // Unit vector pointing from center to point

    double lon = atan2(unit[1], unit[0]);           // longitude
    double lat = atan2(unit[2], sqrt(unit[0]*unit[0] + unit[1]*unit[1]));    // latitude

    std::vector<double> values {lat, lon, R};
    return values;
}//====================================================

/**
 *  \brief Convert inertial coordinates to local tangent coordinates
 * 
 *  \param inertPos Object position in inertial, cartesian (x,y,z) coordinates
 *  \param lat Latitude angle, radians, measured from the inertial xy plane;
 *  above the plane (z > 0) is a positive latitude.
 *  \param lon Longitude angle, radians, measured from reference meridian in a positive,
 *  right-handed rotation about inertial z
 *  \param theta_mer Meridian longitude, radians, measured from inertial x in a positive,
 *  right-handed rotation about inertial z. For example, Earth ground stations' longitudes
 *  are relative to the Greenwich Meridian
 *  @return The position of the object in local tangent, cartesian coordinates (East, North, Up)
 *  as viewed from the specified latitude and longitude
 */
std::vector<double> inert2LocalTangent(std::vector<double> inertPos, double lat, double lon, double theta_mer){
    double nu = lon + theta_mer;    // Rotation about z by longitude + greenwich longitude
    Matrix3Rd dcm;
    dcm << -sin(nu),    -cos(nu)*sin(lat),  cos(nu)*cos(lat),
            cos(nu),    -sin(nu)*sin(lat),  sin(nu)*cos(lat),
            0,              cos(lat),               sin(lat);

    Eigen::RowVector3d r_xyz(inertPos[0], inertPos[1], inertPos[2]);    // Position vector in inertial (x, y, z)
    Eigen::RowVector3d r_enz = r_xyz*dcm;                               // Position vector in local tangent (e_E, e_N, e_Z)

    std::vector<double> local {r_enz(0), r_enz(1), r_enz(2)};
    return local;
}//====================================================

/**
 *  \brief Convert local tangent coordinates to inertial coordinates
 * 
 *  \param localPos Object position in local tangent, cartesian (East, North, Up) coordinates
 *  \param lat Latitude angle, radians, measured from the inertial xy plane;
 *  above the plane (z > 0) is a positive latitude.
 *  \param lon Longitude angle, radians, measured from reference meridian in a positive,
 *  right-handed rotation about inertial z
 *  \param theta_mer Meridian longitude, radians, measured from inertial x in a positive,
 *  right-handed rotation about inertial z. For example, Earth ground stations' longitudes
 *  are relative to the Greenwich Meridian
 *  @return The position of the object in inertial, cartesian coordinates (x, y, z)
 *  as viewed from the specified latitude and longitude
 */
std::vector<double> localTangent2Inert(std::vector<double> localPos, double lat, double lon, double theta_mer){
    double nu = lon + theta_mer;    // Rotation about z by longitude + greenwich longitude
    Matrix3Rd dcm;
    dcm << -sin(nu),    -cos(nu)*sin(lat),  cos(nu)*cos(lat),
            cos(nu),    -sin(nu)*sin(lat),  sin(nu)*cos(lat),
            0,              cos(lat),               sin(lat);

    Eigen::Vector3d r_enz(localPos[0], localPos[1], localPos[2]);      // Position vector in local tangent (east, north, up)
    Eigen::Vector3d r_xyz = dcm*r_enz;                                 // Position vector in inertial (x, y, z)

    std::vector<double> inert {r_xyz(0), r_xyz(1), r_xyz(2)};
    return inert;
}//====================================================

/**
 *  \brief Get the local tangent coordinates of an object given its 
 *  azimuth, elevation, and range
 *  @details [long description]
 * 
 *  \param s range distance
 *  \param az azimuth, measured from North toward East, radians
 *  \param el elevation, measured from local horizontal, radians
 *  @return [r_E, r_N, r_Z] the position of the object in local
 *  tangent coordinates (East, North, Up), units that match the input range
 */
std::vector<double> azEl2LocalTangent(double s, double az, double el){
    std::vector<double> cartCoord {
        s*cos(el)*sin(az),      // East component
        s*cos(el)*cos(az),      // North component
        s*sin(el)               // Up component
    };

    return cartCoord;
}//====================================================

//-----------------------------------------------------
//      2BP Utility Functions
//-----------------------------------------------------

/**
 *  \ingroup 2bp
 *  \brief Compute the Keplarian elements at all steps/nodes
 *  of a 2-body trajectory or nodeset
 * 
 *  \param pSet Pointer to a nodeset or trajectory
 */
void r2bp_computeAllKepler(BaseArcset *pSet){
    if(pSet->getSysData()->getType() == SysData_tp::R2BP_SYS){
        const SysData_2bp *pSys = static_cast<const SysData_2bp*>(pSet->getSysData());
        
        for(int n = 0; n < pSet->getNumNodes(); n++){
            Node& node = pSet->getNodeRefByIx(n);
            r2bp_computeKepler(pSys, &node);
        }
    }else{
        printErr("r2bp_computeAllKepler: Arcset object is not associated with the 2BP");
    }
}//====================================================

/**
 *  \ingroup 2bp
 *  \brief Compute the Keplarian orbital elements at a 
 *  specific node
 * 
 *  \param pSys Pointer to the dynamical system the node exists in
 *  \param pNode Pointer to the node
 */
void r2bp_computeKepler(const SysData_2bp *pSys, Node *pNode){
    std::vector<double> q = pNode->getState();

    double r = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);     // radius, km
    double v = sqrt(q[3]*q[3] + q[4]*q[4] + q[5]*q[5]);     // speed, km/s

    double energy = 0.5*v*v - pSys->getMu()/r;               // two-body energy, km^2/s^2
    double a = -pSys->getMu()/(2*energy);                    // semi-majora axis, km

    double h_vec[] = {  q[1]*q[5] - q[2]*q[4],              // specific angular momentum vector h = r x v, km^2/s
                        q[2]*q[3] - q[0]*q[5],
                        q[0]*q[4] - q[1]*q[3] };

    double h = sqrt(h_vec[0]*h_vec[0] + h_vec[1]*h_vec[1] + h_vec[2]*h_vec[2]);     // specific angular momentum km^2/s
    double e = sqrt(1 + 2*energy*h*h/(pSys->getMu() * pSys->getMu()));                // eccentricity, nondimensional
    double p = a*(1 - e*e);                                 // semi-latus rectum, km
    double fpa = acos(h/(r*v));                             // flight path angle, rad
    double ta = acos((p-r)/(r*e));;                         // true anomaly, rad, NOTE: two options from arccos, symmetric about ta = 0

    double h_unit[] = {h_vec[0]/h, h_vec[1]/h, h_vec[2]/h};                 // unit vector, nondimensional
    double r_unit[] = {q[0]/r, q[1]/r, q[2]/r};                             // unit vector, nondimensional
    double theta_unit[] = { h_unit[1]*r_unit[2] - h_unit[2]*r_unit[1],      // unit vector, theta = h x r, nondimensional
                            h_unit[2]*r_unit[0] - h_unit[0]*r_unit[2],
                            h_unit[0]*r_unit[1] - h_unit[1]*r_unit[0]};

    double r_dot = q[3]*r_unit[0] + q[4]*r_unit[1] + q[5]*r_unit[2];        // range rate, km/s

    // if s/c is descending, true anomaly is between pi and 2*pi
    if(r_dot < 0){
        ta = 2*PI - ta;
    }

    double inc = acos(boundValue(h_unit[2], -1.0, 1.0));   // inclination, rad, choose positive
    
    // Right ascension of the ascending node, rad
    double raan = resolveAngle(asin(boundValue(h_unit[0]/sin(inc), -1.0, 1.0)),
        acos(boundValue(-h_unit[1]/sin(inc), -1.0, 1.0)));

    // Angular distance from the ascending node, rad
    double theta = resolveAngle(asin(boundValue(r_unit[2]/sin(inc), -1.0, 1.0)),
        acos(boundValue(theta_unit[2]/sin(inc), -1.0, 1.0)));

    pNode->setExtraParam("range", r);
    pNode->setExtraParam("speed", v);
    pNode->setExtraParam("energy", energy);
    pNode->setExtraParam("sma", a);
    pNode->setExtraParam("angMom", h);
    pNode->setExtraParam("ecc", e);
    pNode->setExtraParam("fpa", fpa);
    pNode->setExtraParam("ta", ta);
    pNode->setExtraParam("rangerate", r_dot);
    pNode->setExtraParam("inc", inc);
    pNode->setExtraParam("raan", raan);
    pNode->setExtraParam("theta", theta);
}//====================================================

/**
 *  \ingroup 2bp
 *  \brief [brief description]
 *  @details [long description]
 * 
 *  \param pSys [description]
 *  \param a [description]
 *  \param e [description]
 *  \param argPeri [description]
 *  \param i [description]
 *  \param RAAN [description]
 *  \param TA [description]
 *  @return [description]
 */
std::vector<double> r2bp_stateFromKepler(const SysData_2bp *pSys, double a, double e, double argPeri, double i, double RAAN, double TA){
    if(a < 0)
        throw Exception("r2bp_stateFromKepler: a < 0 -> Code is not equiped to handle hyperbolic orbits");
    if(e >= 1)
        throw Exception("r2bp_stateFromKepler: e >=1 -> Code is not equipped to handle parabolic or hyperbolic orbits");

    double r = a*(1 - e*e)/(1 + e*cos(TA));             // radius, km
    double v = sqrt(pSys->getMu()*(2.0/r - 1.0/a));     // speed, km/s
    double h = sqrt(pSys->getMu()*a*(1 - e*e));         // specific angular momentum, km^2/s
    double fpa = acos(boundValue(h/(r*v), -1.0, 1.0));  // flight path angle, rad

    // double fpa = 0;
    // if(e != 0)  // When eccentricity is zero, can get nan from arccos
    //     fpa = acos(h/(r*v));                         

    // if true anomaly is between 180 and 360 degrees, make flight path angle negative
    if(std::fmod(TA, 2.0*PI) > PI)
        fpa *= -1;

    // Put velocity in terms of r-theta-h coordinates
    double v_r = v*sin(fpa);
    double v_theta = v*cos(fpa);

    double theta = TA+argPeri;

    // DCM that converts r-theta-h coordinates to x-y-z coordinates
    // The first row is (x dot r, x dot theta), second row is (y dot r, y dot theta),
    // third row is (z dot r, z dot theta)
    double dcm[][2] =   {{cos(RAAN)*cos(theta) - sin(RAAN)*cos(i)*sin(theta), -cos(RAAN)*sin(theta) - sin(RAAN)*cos(i)*cos(theta)},
                         {sin(RAAN)*cos(theta) + cos(RAAN)*cos(i)*sin(theta), -sin(RAAN)*sin(theta) + cos(RAAN)*cos(i)*cos(theta)},
                         {sin(i)*sin(theta), sin(i)*cos(theta)}
                        };
    std::vector<double> state(6,0);
    state[0] = r*dcm[0][0];
    state[1] = r*dcm[1][0];
    state[2] = r*dcm[2][0];
    state[3] = v_r*dcm[0][0] + v_theta*dcm[0][1];
    state[4] = v_r*dcm[1][0] + v_theta*dcm[1][1];
    state[5] = v_r*dcm[2][0] + v_theta*dcm[2][1];

    printf("r = %.4f km\n", r);
    printf("v = %.4f km/s\n", v);
    printf("h = %.4f km^2/s\n", h);
    printf("fpa = %.4f deg\n", fpa*180/PI);
    printf("theta = %.4f deg\n", theta*180/PI);

    return state;
}//====================================================

//-----------------------------------------------------
//      CR3BP Utility Functions
//-----------------------------------------------------

/**
 *  \ingroup cr3bp
 *  \brief Compute the magnitude of a velocity component given Jacobi Constant
 * 
 *  \param s state vector (non-dimensional); MUST contain at least the 6 position and velocity states.
 *  Note that the desired velocity component identified by velIxToFind is not required; put a placeholder
 *  zero or NAN (or anything, really) in its place; this value will not be used in computations.
 *  \param mu non-dimensional system mass ratio
 *  \param C Jacobi constant value
 *  \param velIxToFind index of the velocity component to compute (i.e. 3, 4, or 5)
 *  @return the magnitude of vx (3), vy (4), or vz (5).
 *  \throws Exception if <tt>velIxToFind</tt> is out of bounds
 */
double cr3bp_getVel_withC(const double s[], double mu, double C, int velIxToFind){
    double v_squared = 0;
    switch(velIxToFind){
        case 3: v_squared = s[4]*s[4] + s[5]*s[5]; break;
        case 4: v_squared = s[3]*s[3] + s[5]*s[5]; break;
        case 5: v_squared = s[3]*s[3] + s[4]*s[4]; break;
        default: throw Exception("Utilities::cr3bp_getVel_withC: velocity index is invalid\n");
    }
    
    double d = sqrt((s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2]);
    double r = sqrt((s[0] - 1 + mu)*(s[0] - 1 + mu) + s[1]*s[1] + s[2]*s[2]);
    double U = (1 - mu)/d + mu/r + 0.5*(s[0]*s[0] + s[1]*s[1]);

    // Solve for the desired velocity component
    return sqrt(2*U - v_squared - C);
}//=================================================

/**
 *  \ingroup cr3bp
 *  \brief Compute a periodic orbit in the CR3BP system
 *
 *  The initial and final states are constrained based on the mirror type, but
 *  all other states are allowed to vary during the corrections process. If you
 *  wish to constrain a specific states, see 
 *  cr3bp_getPeriodic(pSys, IC, period, numNodes, mirrorType, fixedStates)
 *  This function also uses only two nodes; to specify more, see the function above.
 *
 *  This function also assumes the order of the periodic orbit is 1
 *
 *  \param pSys the dynamical system
 *  \param IC non-dimensional initial state vector
 *  \param period non-dimensional period for the orbit
 *  \param mirrorType how this periodic orbit mirrors in the CR3BP
 *  \param tol tolerance to use in the corrections process
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 */
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp *pSys, std::vector<double> IC,
    double period, Mirror_tp mirrorType, double tol){
    
    std::vector<int> fixedStates;   // Initialize an empty vector
    return cr3bp_getPeriodic(pSys, IC, period, 2, 1, mirrorType, fixedStates, tol);
}//========================================

/**
 *  \ingroup cr3bp
 *  \brief Compute a periodic orbit in the CR3BP system
 *  @details This method ignores all crash events, so it is possible to compute a
 *  periodic orbit that passes through a primary
 *  
 *  \param pSys the dynamical system
 *  \param IC non-dimensional initial state vector
 *  \param period non-dimensional period for the orbit
 *  \param numNodes the number of nodes to use for HALF of the periodic orbit; more nodes 
 *  may result in a more robust correction
 *  \param order the number of revolutions about the system/primary this orbit completes before
 *  it repeats periodically. Think of a period-3 DRO (order = 3) or a butterfly (order = 2)
 *  \param mirrorType how this periodic orbit mirrors in the CR3BP
 *  \param fixedStates a vector containing the indices of which initial states
 *  we would like to fix. Not all states are possible for each mirror condition.
 *  See the enum definition for specific details.
 *  \param tol tolerance to use in the corrections process
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 */
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp *pSys, std::vector<double> IC,
    double period, int numNodes, int order, Mirror_tp mirrorType, std::vector<int> fixedStates,
    double tol){


    // return cr3bp_getPeriodic(pSys, IC, period, numNodes, order, mirrorType, fixedStates, tol, &itData);
    return cr3bp_getPeriodic(pSys, IC, period, numNodes, order, mirrorType, fixedStates, tol, nullptr);
}//====================================================================

/**
 *  \ingroup cr3bp
 *  \brief Compute a periodic orbit in the CR3BP system
 *  @details This method ignores all crash events, so it is possible to compute a
 *  periodic orbit that passes through a primary
 *  
 *  \param pSys the dynamical system
 *  \param IC non-dimensional initial state vector
 *  \param period non-dimensional period for the orbit
 *  \param numNodes the number of nodes to use for HALF of the periodic orbit; more nodes 
 *  may result in a more robust correction
 *  \param order the number of revolutions about the system/primary this orbit completes before
 *  it repeats periodically. Think of a period-3 DRO (order = 3) or a butterfly (order = 2)
 *  \param mirrorType how this periodic orbit mirrors in the CR3BP
 *  \param fixedStates a vector containing the indices of which initial states
 *  we would like to fix. Not all states are possible for each mirror condition.
 *  See the enum definition for specific details.
 *  \param tol tolerance to use in the corrections process
 *  \param pItData a pointer to an iteration data object that contains data from the
 *  multiple shooting run that corrects the periodic orbit; this is useful when
 *  attempting to determine how well (or poorly) the multiple shooting algorithm
 *  performed (e.g., for a variable step-size process)
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 *  \throws Exception if <tt>mirrorType</tt> is invalid
 *  \throws DivergeException if the multiple shooting algorithm cannot converge on a 
 *  mirrored solution.
 */
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp *pSys, std::vector<double> IC,
    double period, int numNodes, int order, Mirror_tp mirrorType, std::vector<int> fixedStates,
    double tol, MultShootData* pItData){

    SimEngine sim;    // Engine to perform simulation
    sim.setAbsTol(tol < 1e-12 ? 1e-15 : tol/1000.0);
    sim.setRelTol(sim.getAbsTol());
    sim.setMakeDefaultEvents(false);      // Ignore any crashes into the primaries
    std::vector<int> zeroStates;        // Which states must be zero to ensure a perpendicular crossing

    Event mirrorEvt;
    // Determine which states must be zero for mirroring
    switch(mirrorType){
        case Mirror_tp::MIRROR_XZ:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(3);    // x-dot
            zeroStates.push_back(5);    // z-dot
            mirrorEvt.createEvent(Event_tp::XZ_PLANE, 0, true);     // Tell the sim to quit once it reaches the XZ plane
            break;
        case Mirror_tp::MIRROR_YZ:
            zeroStates.push_back(0);    // x
            zeroStates.push_back(4);    // y-dot
            zeroStates.push_back(5);    // z-dot
            mirrorEvt.createEvent(Event_tp::YZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_XY:
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            zeroStates.push_back(4);    // y-dot
            mirrorEvt.createEvent(Event_tp::YZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_X_AX_H:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            mirrorEvt.createEvent(Event_tp::XZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_X_AX_V:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            mirrorEvt.createEvent(Event_tp::XY_PLANE, 0, true);
            break;
        default:
            throw Exception("Mirror type either not defined or not implemented");
    }
    // sim.setVerbosity(true);
    mirrorEvt.setStopCount(order);
    sim.addEvent(mirrorEvt);

    // Create a constraint to enforce mirror condition at the beginning and end of arc
    double mirrorCon0[] = {NAN,NAN,NAN,NAN,NAN,NAN};
    double mirrorCon1[] = {NAN,NAN,NAN,NAN,NAN,NAN};

    // Set the zero states to zero in both constraints
    for(unsigned int i = 0; i < zeroStates.size(); i++){
        mirrorCon0[zeroStates[i]] = 0;
        mirrorCon1[zeroStates[i]] = 0;
    }

    // Fix states at the initial point
    for(unsigned int i = 0; i < fixedStates.size(); i++){
        bool okToFix = true;
        for(unsigned int n = 0; n < zeroStates.size(); n++){
            if(fixedStates[i] == zeroStates[n]){
                astrohelion::printWarn("Cannot fix state %d; it must be zero for this mirror condition; ignoring\n", fixedStates[i]);
                okToFix = false;
                break;
            }
        }
        if(okToFix){
            mirrorCon0[fixedStates[i]] = IC[fixedStates[i]];
        }
    }

    Constraint initStateCon(Constraint_tp::STATE, 0, mirrorCon0, 6);
    Constraint finalStateCon(Constraint_tp::STATE, numNodes-1, mirrorCon1, 6);

    // Run the sim until the event is triggered
    Traj_cr3bp halfOrbArc(pSys);
    sim.runSim(IC, period, &halfOrbArc);
    
    // Check to make sure the simulation ended with the event (not running out of time)
    std::vector<Event> endEvts = sim.getEndEvents(&halfOrbArc);
    if(std::find(endEvts.begin(), endEvts.end(), mirrorEvt) == endEvts.end()){
        astrohelion::printErr("Calculations::cr3bp_getPeriodic: simulation of half-period orbit did not end in mirror event; may have diverged\n");
        halfOrbArc.saveToMat("halfOrbArc_failedMirror.mat");
    }
    // halfOrbArc.saveToMat("HalfOrbArc.mat");

    double halfOrbTOF = halfOrbArc.getTimeByIx(-1);
    double tofErr = 100*std::abs(halfOrbTOF-period/2.0)/(period/2.0);

    if(tofErr > 10)
        astrohelion::printWarn("Calculations::cr3bp_getPeriodic: Half-Period arc TOF varies from input half-period by more than 10%%\n");

    // Create a nodeset from arc
    Nodeset_cr3bp halfOrbNodes(halfOrbArc, numNodes, Nodeset::TIME);
    halfOrbNodes.addConstraint(initStateCon);
    halfOrbNodes.addConstraint(finalStateCon);

    // Use differential corrections to enforce the mirror conditions
    CorrectionEngine corrector;
    corrector.setTol(tol);
    corrector.setIgnoreCrash(true); // Corrector also ignores crash events
    corrector.setVarTime(true);
    corrector.setEqualArcTime(true);
    // corrector.setVerbosity(Verbosity_tp::ALL_MSG);

    try{
        Nodeset_cr3bp correctedHalfPer(pSys);

        // If the user passed in a null pointer, create a temporory data object
        // on the stack to avoid seg faults, then delete it before exiting to 
        // avoid memory leaks
        bool createdTempMSData = false;
        if(!pItData){
            createdTempMSData = true;
            pItData = new MultShootData(&halfOrbNodes);
        }

        *pItData = corrector.multShoot(&halfOrbNodes, &correctedHalfPer);
        // correctedHalfPer.saveToMat("temp_correctedHalfPer.mat");

        // Make the nodeset into a trajectory
        Traj_cr3bp halfPerTraj = Traj_cr3bp::fromNodeset(correctedHalfPer);
        double halfTOF = halfPerTraj.getTimeByIx(-1);
        double halfPerTraj_len = halfPerTraj.getNumNodes();
        MatrixXRd halfPerSTM = halfPerTraj.getSTMByIx(-1);
        
        // Use Mirror theorem to create the second half of the orbit
        MatrixXRd mirrorMat = getMirrorMat(mirrorType);
        int prevID = halfPerTraj.getNodeByIx(halfPerTraj_len-1).getID();
        // halfPerTraj.saveToMat("temp_halfPerTraj.mat");
        for(int i = halfPerTraj_len-2; i >= 0; i--){
            // Use mirroring to populate second half of the orbit
            std::vector<double> state = halfPerTraj.getStateByIx(i);
            Eigen::RowVectorXd stateVec = Eigen::Map<Eigen::RowVectorXd>(&(state[0]), 1, 6);
            Eigen::RowVectorXd newStateVec = stateVec*mirrorMat;

            Node node;
            node.setState(newStateVec.data(), newStateVec.cols());
            node.setEpoch(2*halfTOF - halfPerTraj.getTimeByIx(i));

            int id = halfPerTraj.addNode(node);
            // printf("Added node at epoch %.6f\n", node.getEpoch());
            // fprintf("Adding segment with tof = %.6f\n", halfPerTraj.getEpoch(id) - halfPerTraj.getEpoch(prevID))
            Segment seg(prevID, id, halfPerTraj.getEpoch(id) - halfPerTraj.getEpoch(prevID));
            prevID = id;
            halfPerTraj.addSeg(seg);
        }

        // waitForUser();
        // Compute the monodromy matrix from the half-period STM
        double M_data[] = { 0, 0, 0, -1, 0, 0,
                            0, 0, 0, 0, -1, 0,
                            0, 0, 0, 0, 0, -1,
                            1, 0, 0, 0, -2, 0,
                            0, 1, 0, 2, 0, 0,
                            0, 0, 1, 0, 0, 0};
        // Inverse of M
        double MI_data[] = {0, -2, 0, 1, 0, 0,
                            2, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, 1,
                            -1, 0, 0, 0, 0, 0,
                            0, -1, 0, 0, 0, 0,
                            0, 0, -1, 0, 0, 0};
        MatrixXRd M = Eigen::Map<MatrixXRd>(M_data, 6, 6);
        MatrixXRd MI = Eigen::Map<MatrixXRd>(MI_data, 6, 6);
        
        MatrixXRd monoMat(6,6);
        monoMat.noalias() = mirrorMat*M*halfPerSTM.transpose()*MI*mirrorMat*halfPerSTM;

        // Set final STM of mirrored trajectory to the one computed here
        halfPerTraj.setSTMByIx(-1, monoMat);
        
        if(createdTempMSData){
            delete pItData;
            pItData = nullptr;
        }
        return halfPerTraj;     // Now contains entire trajectory
    }catch(DivergeException &e){
        throw DivergeException("Calculations::cr3bp_getPeriodic: Could not converge half-period arc with mirroring condition");
    }
}//================================================

/**
 *  \ingroup cr3bp
 *  \brief Transition a trajectory from the EM system to the SE system
 *  
 *  The relative orientation between the two systems is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. These angles describe the
 *  system orientation at time t = 0, so if the start of the trajectory does not coincide
 *  with this initial time, be sure to shift the time coordinates to assure that the first
 *  value in the time vector reflects correct initial time.
 *
 *  \param EMTraj a CR3BP Earth-Moon trajectory
 *  \param pSESys a Sun-Earth CR3BP system data object
 *  \param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  \param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  \param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Traj_cr3bp cr3bp_EM2SE(Traj_cr3bp EMTraj, const SysData_cr3bp *pSESys, double thetaE0, double thetaM0, double gamma){
    printf("cr3bp_EM2SE\n");
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp seNodes = cr3bp_EM2SE(static_cast<Nodeset_cr3bp>(EMTraj), pSESys, thetaE0, thetaM0, gamma);
    return static_cast<Traj_cr3bp>(seNodes);
}//=========================================================

/**
 *  \ingroup cr3bp
 *  \brief Transition a nodeset from the EM system to the SE system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. Accordingly, the EM nodeset
 *  should have node epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire set may be accomplished via the updateEpochs() function.
 *
 *  \param EMNodes a CR3BP Earth-Moon nodeset
 *  \param pSESys a Sun-Earth CR3BP system data object
 *  \param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  \param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  \param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Nodeset_cr3bp cr3bp_EM2SE(Nodeset_cr3bp EMNodes, const SysData_cr3bp *pSESys, double thetaE0, double thetaM0,
    double gamma){

    // astrohelion::printColor(BLUE, "Converting EM to SE\nEM Sys:\n  %d Nodes\n  %d Segments\n", EMNodes.getNumNodes(),
    //     EMNodes.getNumSegs());

    Nodeset_cr3bp SENodes(pSESys);

    double charTE = EMNodes.getSysData()->getCharT();       // characteristic time in EM system
    double charLE = EMNodes.getSysData()->getCharL();       // characteristic length in EM system
    double charTS = pSESys->getCharT();                       // characteristic time in SE system
    double charLS = pSESys->getCharL();                       // characteristic length in SE system

    std::vector<int> map_oldID_to_newID(EMNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> state_SE;
    double epoch;

    for(int n = 0; n < EMNodes.getNumNodes(); n++){
        epoch = EMNodes.getEpochByIx(n);
        state_SE = cr3bp_EM2SE_state(EMNodes.getStateByIx(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, pSESys->getMu());
        
        map_oldID_to_newID[EMNodes.getNodeByIx(n).getID()] = SENodes.addNode(Node(state_SE, epoch*charTE/charTS));
    }

    // Convert time units for segment TOFs and update IDs of origin and terminus nodes
    for(int s = 0; s < EMNodes.getNumSegs(); s++){
        Segment seg = EMNodes.getSegByIx(s);
        seg.setTOF(seg.getTOF()*charTE/charTS);
        seg.setSTM(MatrixXRd::Identity(6, 6));

        // Remap the origin and terminus to the new IDs
        if(seg.getOrigin() != Linkable::INVALID_ID)
            seg.setOrigin(map_oldID_to_newID[seg.getOrigin()]);

        if(seg.getTerminus() != Linkable::INVALID_ID)
            seg.setTerminus(map_oldID_to_newID[seg.getTerminus()]);

        SENodes.addSeg(seg);
    }

    return SENodes;
}//=========================================================

/**
 *  \ingroup cr3bp
 *  \brief Transition a trajectory from the SE system to the EM system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. Accordingly, the SE trajectory
 *  should have epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire trajectory may be accomplished via the updateEpochs() function.
 *
 *  \param SETraj a CR3BP Sun-Earth trajectory
 *  \param pEMSys an Earth-Moon CR3BP system data object
 *  \param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  \param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  \param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Traj_cr3bp cr3bp_SE2EM(Traj_cr3bp SETraj, const SysData_cr3bp *pEMSys, double thetaE0, double thetaM0, double gamma){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp emNodes = cr3bp_SE2EM(static_cast<Nodeset_cr3bp>(SETraj), pEMSys, thetaE0, thetaM0, gamma);
    return static_cast<Traj_cr3bp>(emNodes);
}//=========================================================

/**
 *  \ingroup cr3bp
 *  \brief Transition a nodeset from the SE system to the EM system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. Accordingly, the SE nodeset
 *  should have node epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire set may be accomplished via the updateEpochs() function.
 *
 *  \param SENodes a CR3BP Sun-Earth nodeset
 *  \param pEMSys an Earth-Moon CR3BP system data object
 *  \param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  \param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  \param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Nodeset_cr3bp cr3bp_SE2EM(Nodeset_cr3bp SENodes, const SysData_cr3bp *pEMSys, double thetaE0, double thetaM0,
    double gamma){

    // astrohelion::printColor(BLUE, "Converting SE to EM\nSE Sys:\n  %d Nodes\n  %d Segments\n", SENodes.getNumNodes(),
    //     SENodes.getNumSegs());

    Nodeset_cr3bp EMNodes(pEMSys);

    double charTE = pEMSys->getCharT();                   // characteristic time in EM system
    double charLE = pEMSys->getCharL();                   // characteristic length in EM system
    double charTS = SENodes.getSysData()->getCharT();    // characteristic time in SE system
    double charLS = SENodes.getSysData()->getCharL();    // characteristic length in SE system

    const SysData_cr3bp *pSESys = static_cast<const SysData_cr3bp*>(SENodes.getSysData());

    std::vector<int> map_oldID_to_newID(SENodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> state_EM;
    double epoch;

    for(int n = 0; n < SENodes.getNumNodes(); n++){
        epoch = SENodes.getEpochByIx(n);
        // Transform a single node
        state_EM = cr3bp_SE2EM_state(SENodes.getStateByIx(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, pSESys->getMu());
        
        map_oldID_to_newID[SENodes.getNodeByIx(n).getID()] = EMNodes.addNode(Node(state_EM, epoch*charTS/charTE));
    }

    // Convert time units for segment TOFs and update IDs of origin and terminus nodes
    for(int s = 0; s < SENodes.getNumSegs(); s++){
        Segment seg = SENodes.getSegByIx(s);
        seg.setTOF(seg.getTOF()*charTS/charTE);
        seg.setSTM(MatrixXRd::Identity(6, 6));

        // Remap the origin and terminus to the new IDs
        if(seg.getOrigin() != Linkable::INVALID_ID)
            seg.setOrigin(map_oldID_to_newID[seg.getOrigin()]);

        if(seg.getTerminus() != Linkable::INVALID_ID)
            seg.setTerminus(map_oldID_to_newID[seg.getTerminus()]);

        EMNodes.addSeg(seg);
    }
    return EMNodes;
}//=========================================================

/**
 *  \ingroup cr3bp
 *  \brief Transform a single state from EM coordinates to SE coordinates
 *
 *  \param state_EM a 6- or 9-element state vector
 *  \param t Earth-Moon nondimensional time associated with the Earth-Moon state
 *  \param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  \param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  \param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 *  \param charLE EM characteristic length
 *  \param charTE EM characterstic time
 *  \param charLS SE characteristic length
 *  \param charTS SE characteristic time
 *  \param mu_SE SE mass ratio
 *
 *  @return a 6-element state vector in EM coordinates
 */
std::vector<double> cr3bp_EM2SE_state(std::vector<double> state_EM, double t, double thetaE0, double thetaM0,
    double gamma, double charLE, double charTE, double charLS, double charTS, double mu_SE){

    // Shift coordinates to SE barcyenter from EM barycenter
    Matrix3Rd I = Matrix3Rd::Identity(3,3);
    Eigen::RowVector3d posShift = I.row(0)*(1 - mu_SE);

    // Compute Earth's position in inertial frame at this time
    double thetaE_k = thetaE0 + t*charTE/charTS;
    double thetaM_k = thetaM0 + t;

    // Compute DCMs
    double inertToLunar[] = {cos(gamma), 0, sin(gamma), 0, 1, 0, -sin(gamma), 0, cos(gamma)};
    double inertToSE[] = {cos(thetaE_k), -sin(thetaE_k), 0, sin(thetaE_k), cos(thetaE_k), 0,
                                0, 0, 1};
    double lunarToEM[] = {cos(thetaM_k), -sin(thetaM_k), 0, sin(thetaM_k), cos(thetaM_k), 0,
                                0, 0, 1};

    Matrix3Rd DCM_I2L = Eigen::Map<Matrix3Rd>(inertToLunar);
    Matrix3Rd DCM_I2S = Eigen::Map<Matrix3Rd>(inertToSE);
    Matrix3Rd DCM_L2E = Eigen::Map<Matrix3Rd>(lunarToEM);
    
    Eigen::RowVector3d posEM = Eigen::Map<Eigen::Vector3d>(&(state_EM[0]));
    Eigen::RowVector3d velEM = Eigen::Map<Eigen::Vector3d>(&(state_EM[0])+3);

    // Rotate the position into SE frame and shift basepoint to SE barycenter (EM working frame)
    Eigen::RowVector3d posSE = posEM*DCM_L2E.transpose()*DCM_I2L.transpose()*DCM_I2S + posShift*charLS/charLE;
    

    // Angular velocity of SE frame in EM frame (working frame = EM)
    Eigen::RowVector3d omega = -1*charTE/charTS*I.row(2)*DCM_I2S.transpose()*DCM_I2L*DCM_L2E + 
        charTE/charTE*I.row(2);

    // Apply BKE to get velocity with SE observer, still in EM working frame
    Eigen::RowVector3d velSE = velEM + omega.cross(posEM);

    // Rotate the velocity into the SE working frame
    velSE *= DCM_L2E.transpose()*DCM_I2L.transpose()*DCM_I2S;

    // Units are still in non-dim EM, so change to SE non-dim
    posSE *= charLE/charLS;
    velSE *= (charLE/charTE)/(charLS/charTS);

    // Put new data into state vector
    std::vector<double> state_SE;
    state_SE.insert(state_SE.begin(), posSE.data(), posSE.data()+3);
    state_SE.insert(state_SE.begin()+3, velSE.data(), velSE.data()+3);

    return state_SE;
}//====================================================

/**
 *  \ingroup cr3bp
 *  \brief Transform a single state from SE coordinates to EM coordinates
 *
 *  \param state_SE a 6- or 9-element state vector
 *  \param t Sun-Earth nondimensional time associated with the state
 *  \param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  \param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  \param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 *  \param charLE EM characteristic length, km
 *  \param charTE EM characterstic time, sec
 *  \param charLS SE characteristic length, km
 *  \param charTS SE characteristic time, sec
 *  \param mu_SE SE mass ratio, nondimensional
 *
 *  @return a 6-element state vector in SE coordinates
 */
std::vector<double> cr3bp_SE2EM_state(std::vector<double> state_SE, double t, double thetaE0, double thetaM0,
    double gamma, double charLE, double charTE, double charLS, double charTS, double mu_SE){
    
    Matrix3Rd I = Matrix3Rd::Identity(3,3);
    Eigen::RowVector3d posShift = I.row(0)*(1 - mu_SE);

    // Compute Earth's position in inertial frame at this time
    double thetaE_k = thetaE0 + t;
    double thetaM_k = thetaM0 + t*charTS/charTE;

    // Compute DCMs
    double inertToLunar[] = {cos(gamma), 0, sin(gamma), 0, 1, 0, -sin(gamma), 0, cos(gamma)};
    double inertToSE[] = {cos(thetaE_k), -sin(thetaE_k), 0, sin(thetaE_k), cos(thetaE_k), 0,
                                0, 0, 1};
    double lunarToEM[] = {cos(thetaM_k), -sin(thetaM_k), 0, sin(thetaM_k), cos(thetaM_k), 0,
                                0, 0, 1};

    Matrix3Rd DCM_I2L = Eigen::Map<Matrix3Rd>(inertToLunar);
    Matrix3Rd DCM_I2S = Eigen::Map<Matrix3Rd>(inertToSE);
    Matrix3Rd DCM_L2E = Eigen::Map<Matrix3Rd>(lunarToEM);
    
    Eigen::RowVector3d posSE = Eigen::Map<Eigen::Vector3d>(&(state_SE[0]));
    Eigen::RowVector3d velSE = Eigen::Map<Eigen::Vector3d>(&(state_SE[0])+3);

    // Rotate the position into SE frame, coordinates are EM ND
    Eigen::RowVector3d posEM = (posSE - posShift)*DCM_I2S.transpose()*DCM_I2L*DCM_L2E;

    // Angular velocity of EM frame in SE frame (working frame = SE)
    Eigen::RowVector3d omega = charTS/charTS*I.row(2) - charTS/charTE*I.row(2)*DCM_I2L.transpose()*DCM_I2S;

    // Compute velocity in EM frame (working frame = SE)
    Eigen::RowVector3d velEM = velSE + omega.cross(posSE - posShift);

    // Rotate the velocity into the EM working frame
    velEM *= DCM_I2S.transpose()*DCM_I2L*DCM_L2E;

    // Units are still in non-dim SE, so change to EM non-dim
    posEM *= charLS/charLE;
    velEM *= (charLS/charTS)/(charLE/charTE);

    // Put new data into state vector
    std::vector<double> state_EM;
    state_EM.insert(state_EM.begin(), posEM.data(), posEM.data()+3);
    state_EM.insert(state_EM.begin()+3, velEM.data(), velEM.data()+3);

    return state_EM;
}//====================================================

/**
 *  \ingroup cr3bp
 *  \brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  \param traj CR3BP trajectory
 *  *  \param epoch0 Epoch associated with t = 0 in the CR3BP; epoch in seconds past J2000 (i.e., ephemeris time)
 *  \param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP trajectory is input, selecting 
 *  <tt>centerIx = 1</tt> will produce Earth-centered inertial coordinates and selecting
 *  <tt>centerIx = 2</tt> will produce Moon-centered inertial coordinates. A choice of
 *  <tt>centerIx = 0</tt> produces system barycenter-centered inertial coordinates.
 * 
 *  @return A trajectory centered around the specified index in inertial, dimensional coordinates
 */
Traj_cr3bp cr3bp_rot2inert(Traj_cr3bp traj, double epoch0, int centerIx){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp inertNodes = cr3bp_rot2inert(static_cast<Nodeset_cr3bp>(traj), epoch0, centerIx);
    return static_cast<Traj_cr3bp>(inertNodes);
}//==================================================================

/**
 *  \ingroup cr3bp
 *  \brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  \param nodes CR3BP nodeset
 *  \param epoch0 Epoch associated with t = 0 in the CR3BP; epoch in seconds past J2000 (i.e., ephemeris time)
 *  \param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP trajectory is input, selecting 
 *  <tt>centerIx = 1</tt> will produce Earth-centered inertial coordinates and selecting
 *  <tt>centerIx = 2</tt> will produce Moon-centered inertial coordinates. A choice of
 *  <tt>centerIx = 0</tt> produces system barycenter-centered inertial coordinates.
 * 
 *  @return A nodeset centered around the specified index in inertial, dimensional coordinates
 */
Nodeset_cr3bp cr3bp_rot2inert(Nodeset_cr3bp nodes, double epoch0, int centerIx){

    if(centerIx < 0 || centerIx > 2){
        throw Exception("Calculations::cr3bp_rot2inert: Invalid center index");
    }

    const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp *>(nodes.getSysData());

    Nodeset_cr3bp inertNodes(pSys);
    std::vector<int> map_oldID_to_newID(nodes.getNextNodeID(), Linkable::INVALID_ID);

    for(int i = 0; i < nodes.getNumNodes(); i++){
        std::vector<double> stateInert = cr3bp_rot2inert_state(nodes.getStateByIx(i), pSys, nodes.getEpochByIx(i), epoch0, centerIx);
        
        // Save the new ID and add the node, making sure to update Epoch as well!
        map_oldID_to_newID[nodes.getNodeByIx(i).getID()] = inertNodes.addNode(Node(stateInert,
            nodes.getEpochByIx(i)*pSys->getCharT() + epoch0));
    }

    for(int s = 0; s < nodes.getNumSegs(); s++){
        Segment seg = nodes.getSegByIx(s);
        seg.setTOF(seg.getTOF()*pSys->getCharT());

        // Remap the origin and terminus to the new IDs
        if(seg.getOrigin() != Linkable::INVALID_ID)
            seg.setOrigin(map_oldID_to_newID[seg.getOrigin()]);

        if(seg.getTerminus() != Linkable::INVALID_ID)
            seg.setTerminus(map_oldID_to_newID[seg.getTerminus()]);

        inertNodes.addSeg(seg);
    }

    return inertNodes;
}//==================================================================

/**
 *  \ingroup cr3bp
 *  \brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  \param stateRot a 6-element non-dimensional state in rotating coordinates
 *  \param pSys a pointer to the CR3BP system the state exists in
 *  \param t the non-dimensional time associated with the input state
 *  \param epoch0 Epoch associated with t = 0 in the CR3BP; epoch in seconds past J2000 (i.e., ephemeris time)
 *  \param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP trajectory is input, selecting 
 *  <tt>centerIx = 1</tt> will produce Earth-centered inertial coordinates and selecting
 *  <tt>centerIx = 2</tt> will produce Moon-centered inertial coordinates. A choice of
 *  <tt>centerIx = 0</tt> produces system barycenter-centered inertial coordinates.
 * 
 *  @return A state centered around the specified point in inertial, ecliptic J2000, dimensional coordinates
 *  @throw Exception if <tt>centerIx</tt> is out of bounds
 */
std::vector<double> cr3bp_rot2inert_state(std::vector<double> stateRot, const SysData_cr3bp *pSys, 
    double t, double epoch0, int centerIx){

    if(centerIx < 0 || centerIx > 2){
        throw Exception("Calculations::cr3bp_rot2inert: Invalid center index");
    }

    ConstSpiceChar *abcorr = "NONE";
    ConstSpiceChar *ref = "ECLIPJ2000";

    std::string str_obs = getNameFromSpiceID(pSys->getPrimID(0));
    std::string str_targ = getNameFromSpiceID(pSys->getPrimID(1));
    ConstSpiceChar *obs = str_obs.c_str();
    ConstSpiceChar *targ = str_targ.c_str();
    SpiceDouble lt, p2State[6];

    std::vector<double> stateInert;
    double r_p2[3], v_p2[3], angMom_p2[3], mag_angMom_p2, inst_charL, inst_charT;
    double unit_x[3], unit_y[3], unit_z[3];

    BodyData P1(pSys->getPrimID(0));
    BodyData P2(pSys->getPrimID(1));
    double sysGM = P1.getGravParam() + P2.getGravParam();

    // Locate P2 relative to P1 at the specified time
    spkezr_c(targ, epoch0 + t*pSys->getCharT(), ref, abcorr, obs, p2State, &lt);
    std::copy(p2State, p2State+3, r_p2);
    std::copy(p2State+3, p2State+6, v_p2);

    // printf("p2 w.r.t. p1 = [%f, %f, %f]\n", p2State[0], p2State[1], p2State[2]);

    // Define angular momentum vector of P2 about P1 (cross product: r x v )
    angMom_p2[0] = r_p2[1]*v_p2[2] - r_p2[2]*v_p2[1];
    angMom_p2[1] = r_p2[2]*v_p2[0] - r_p2[0]*v_p2[2];
    angMom_p2[2] = r_p2[0]*v_p2[1] - r_p2[1]*v_p2[0];
    mag_angMom_p2 = sqrt(angMom_p2[0]*angMom_p2[0] + angMom_p2[1]*angMom_p2[1] + angMom_p2[2]*angMom_p2[2]);

    // printf("angMom_p2 = [%f, %f, %f]\n", angMom_p2[0], angMom_p2[1], angMom_p2[2]);

    // Define instantaneous characteristic length and time
    inst_charL = sqrt(r_p2[0]*r_p2[0] + r_p2[1]*r_p2[1] + r_p2[2]*r_p2[2]);
    inst_charT = sqrt(pow(inst_charL,3)/sysGM);

    // Define instantaneous rotating frame in inertial coordinates
    for(unsigned int j = 0; j < 3; j++){
        unit_x[j] = r_p2[j]/inst_charL;
        unit_z[j] = angMom_p2[j]/mag_angMom_p2;
    }
    
    // y-unit vector is cross(z, x)
    unit_y[0] = unit_z[1]*unit_x[2] - unit_z[2]*unit_x[1];
    unit_y[1] = unit_z[2]*unit_x[0] - unit_z[0]*unit_x[2];
    unit_y[2] = unit_z[0]*unit_x[1] - unit_z[1]*unit_x[0];

    // printf("unit_y = [%f, %f, %f]\n", unit_y[0], unit_y[1], unit_y[2]);

    // Convert rotating state to dimensional coordinates
    for(unsigned int j = 0; j < 6; j++){
        // Shift coordinates to center at P1 or P2
        if(j == 0 && centerIx > 0)
            stateRot[j] -= centerIx == 1 ? -pSys->getMu() : 1 - pSys->getMu();

        // Dimensionalize in length
        stateRot[j] *= inst_charL;

        // Dimensionalize velocities in time
        if(j > 2)
            stateRot[j] /= inst_charT;
    }

    // Transform rotating state into inertial state: project position onto
    // inertial unit vectors; for velocity, make sure to include cross product
    // of angular velocity and position (BKE)
    stateInert.clear();
    stateInert.push_back( unit_x[0]*stateRot[0] + unit_y[0]*stateRot[1] + unit_z[0]*stateRot[2]);   // Inertial X
    stateInert.push_back( unit_x[1]*stateRot[0] + unit_y[1]*stateRot[1] + unit_z[1]*stateRot[2]);   // Inertial Y
    stateInert.push_back( unit_x[2]*stateRot[0] + unit_y[2]*stateRot[1] + unit_z[2]*stateRot[2]);   // Inertial Z
    stateInert.push_back( unit_x[0]*stateRot[3] + unit_y[0]*stateRot[4] + unit_z[0]*stateRot[5] + pow(inst_charL, -2)*(angMom_p2[1]*stateInert[2] - angMom_p2[2]*stateInert[1]));   // Inertial VX
    stateInert.push_back( unit_x[1]*stateRot[3] + unit_y[1]*stateRot[4] + unit_z[1]*stateRot[5] - pow(inst_charL, -2)*(angMom_p2[0]*stateInert[2] - angMom_p2[2]*stateInert[0]));   // Inertial VY
    stateInert.push_back( unit_x[2]*stateRot[3] + unit_y[2]*stateRot[4] + unit_z[2]*stateRot[5] + pow(inst_charL, -2)*(angMom_p2[0]*stateInert[1] - angMom_p2[1]*stateInert[0]));   // Inertial VZ

    // Nondimensionalize using instantaneous characteristic quantities
    // for(unsigned int j = 0; j < 6; j++){
    //     stateInert[j] /= inst_charL;
    //     if(j > 2)
    //         stateInert[j] *= inst_charT;
    // }

    // Matrix3Rd I = Matrix3Rd::Identity(3,3);

    // // Shift from Barycenter-centered coordinates to primary-centered coordinates
    // Eigen::RowVector3d posShift = centerIx == 0 ? I.row(0)*(0-pSys->getMu()) : I.row(0)*(1 - pSys->getMu());
    
    // // Angular velocity of rotating frame relative to inertial, non-dim units
    // Eigen::RowVector3d omegaECI = I.row(2);

    // double rotToInert[] = {cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1};
    // Matrix3Rd DCM_R2I = Eigen::Map<Matrix3Rd>(rotToInert);

    // Eigen::RowVector3d posRot = Eigen::Map<Eigen::Vector3d>(&(stateRot[0]));
    // Eigen::RowVector3d velRot = Eigen::Map<Eigen::Vector3d>(&(stateRot[0])+3);

    // // Shift, rotate, scale coordinates into inertial, dimensional
    // Eigen::RowVector3d posInert = (posRot - posShift)*DCM_R2I*pSys->getCharL();

    // // Transform velocity into inertial, dimensional
    // Eigen::RowVector3d velInert = (velRot + omegaECI.cross(posRot - posShift))*DCM_R2I*pSys->getCharL()/pSys->getCharT();

    // // Concatenate data into a state vector
    // std::vector<double> stateInert;
    // stateInert.insert(stateInert.begin(), posInert.data(), posInert.data()+3);
    // stateInert.insert(stateInert.end(), velInert.data(), velInert.data()+3);

    return stateInert;
}//==================================================================

//-----------------------------------------------------
//      BCR4BP Utility Functions
//-----------------------------------------------------

/**
 *  \brief Convert a CR3BP Sun-Earth trajectory into a BCR4BPR Sun-Earth-Moon trajectory
 *
 *  \param crTraj a CR3BP Sun-Earth trajectory
 *  \param pBCSys a BCR4BPR Sun-Earth-Moon system data object; contains information about system
 *  scaling and orientation at time t = 0
 *  \param nodeID ID of a node in the set for which the epoch is known
 *  \param t0 the epoch at the specified node in BCR4BPR non-dimensional time units
 *
 *  @return a BCR4BPR Trajectory object
 *  @throw Exception if <tt>crTraj</tt> is not a Sun-Earth trajectory or if <tt>pBCSys</tt> is
 *  not the Sun-Earth-Moon system.
 */
Traj_bc4bp bcr4bpr_SE2SEM(Traj_cr3bp crTraj, const SysData_bc4bp *pBCSys, int nodeID, double t0){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_bc4bp bcNodeset = bcr4bpr_SE2SEM(static_cast<Nodeset_cr3bp>(crTraj), pBCSys, nodeID, t0);
    return static_cast<Traj_bc4bp>(bcNodeset);
}//==================================================

/**
 *  \ingroup bc4bp
 *  \brief Convert a Sun-Earth CR3BP Nodeset to a Sun-Earth-Moon BCR4BPR Nodeset.
 *  
 *  \param crNodes a CR3BP Sun-Earth nodeset
 *  \param pBCSys a BCR4BPR Sun-Earth-Moon system data object; contains information about system
 *  scaling and orientation at time t = 0
 *  \param nodeID ID of a node in the set for which the epoch is known
 *  \param t0 the epoch at the specified node in BCR4BPR non-dimensional time units
 *
 *  @return a BCR4BPR nodeset
 *  @throw Exception if <tt>crNodes</tt> is not a Sun-Earth trajectory or if <tt>pBCSys</tt> is
 *  not the Sun-Earth-Moon system.
 */
Nodeset_bc4bp bcr4bpr_SE2SEM(Nodeset_cr3bp crNodes, const SysData_bc4bp *pBCSys, int nodeID, double t0){
    if(crNodes.getSysData()->getPrimID(0) != 10 || crNodes.getSysData()->getPrimID(1) != 399){
        throw Exception("CR3BP trajectory is not in the Sun-Earth System");
    }

    if(pBCSys->getPrimID(0) != 10 || pBCSys->getPrimID(1) != 399 || pBCSys->getPrimID(2) != 301){
        throw Exception("BCR4BPR system is not Sun-Earth-Moon");
    }

    // Create a BCR4BPR Nodeset
    Nodeset_bc4bp bcNodes(pBCSys);

    double charL2 = crNodes.getSysData()->getCharL();
    double charT2 = crNodes.getSysData()->getCharT();
    double charL3 = pBCSys->getCharL();
    double charT3 = pBCSys->getCharT();

    crNodes.updateEpochs(nodeID, t0*charT3/charT2);

    // A mapping vector: index is the old node ID, value is the new node ID
    // All new IDs are initialized to the default INVALID_ID value
    std::vector<int> map_oldID_to_newID(crNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> bcNodeState, crNodeState;
    double epoch;
    for(int n = 0; n < crNodes.getNumNodes(); n++){
        bcNodeState.clear();
        crNodeState = crNodes.getStateByIx(n);
        for(int r = 0; r < static_cast<int>(crNodeState.size()); r++){
            if(r == 0)  // Convert x-coordinate, shift base to P2/P3 Barycenter
                bcNodeState.push_back(crNodeState[r]*charL2/charL3 - (1.0/pBCSys->getK() - pBCSys->getMu()));
            else if(r < 3)   // Convert position
                bcNodeState.push_back(crNodeState[r]*charL2/charL3);
            else  // Convert velocity
                bcNodeState.push_back(crNodeState[r]*(charL2/charL3)*(charT3/charT2));
        }
        epoch = crNodes.getEpochByIx(n)*charT2/charT3;

        // Save the new ID
        map_oldID_to_newID[crNodes.getNodeByIx(n).getID()] = bcNodes.addNode(Node(bcNodeState, epoch));
    }

    for(int s = 0; s < crNodes.getNumSegs(); s++){
        Segment seg = crNodes.getSegByIx(s);
        seg.setTOF(seg.getTOF()*charT2/charT3);

        // Remap the origin and terminus to the new IDs
        if(seg.getOrigin() != Linkable::INVALID_ID)
            seg.setOrigin(map_oldID_to_newID[seg.getOrigin()]);

        if(seg.getTerminus() != Linkable::INVALID_ID)
            seg.setTerminus(map_oldID_to_newID[seg.getTerminus()]);

        bcNodes.addSeg(seg);
    }

    return bcNodes;
}//==================================================

/**
 *  \ingroup bc4bp
 *  \brief Transform a BCR4BPR Sun-Earth-Moon nodeset into a CR3BP Sun-Earth nodeset
 *
 *  \param bcNodes a BCR4BPR Sun-Earth-Moon nodeset
 *  \param pCRSys a CR3BP Sun-Earth system data object; contains information about system
 *  scaling
 *
 *  @return a CR3BP nodeset object
 *  @throw Exception if <tt>pCRSys</tt> is not a Sun-Earth system or if <tt>bcNodes</tt> is
 *  not in the Sun-Earth-Moon system.
 */
Nodeset_cr3bp bcr4bpr_SEM2SE(Nodeset_bc4bp bcNodes, const SysData_cr3bp *pCRSys){
    if(pCRSys->getPrimID(0) != 10 || pCRSys->getPrimID(1) != 399){
        throw Exception("CR3BP system is not the Sun-Earth System");
    }

    if(bcNodes.getSysData()->getPrimID(0) != 10 || bcNodes.getSysData()->getPrimID(1) != 399 ||
            bcNodes.getSysData()->getPrimID(2) != 301){
        throw Exception("BCR4BPR trajectory is not in Sun-Earth-Moon");
    }

    const SysData_bc4bp *pBCSys = static_cast<const SysData_bc4bp *>(bcNodes.getSysData());

    // Create a BCR4BPR Trajectory
    Nodeset_cr3bp crNodes(pCRSys);

    double charL2 = pCRSys->getCharL();
    double charT2 = pCRSys->getCharT();
    double charL3 = bcNodes.getSysData()->getCharL();
    double charT3 = bcNodes.getSysData()->getCharT();

    std::vector<int> map_oldID_to_newID(bcNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> bcState, crNodeState;
    for(int n = 0; n < bcNodes.getNumNodes(); n++){
        bcState = bcNodes.getStateByIx(n);
        crNodeState.clear();

        for(int r = 0; r < static_cast<int>(bcState.size()); r++){
            if(r == 0)  // Shift origin to P2-P3 barycenter
                crNodeState.push_back((bcState[r] + 1.0/pBCSys->getK() - pBCSys->getMu())*charL3/charL2);
            else if(r < 3)   // Convert position
                crNodeState.push_back(bcState[r]*charL3/charL2);
            else if(r < 6)  // Convert velocity
                crNodeState.push_back(bcState[r]*(charL3/charL2)*(charT2/charT3));
        }

        // Use time relative to Epoch0 defined in pBCSys
        Node node(crNodeState, bcNodes.getEpochByIx(n)*charT3/charT2);
        // Node node(crNodeState, (bcNodes.getEpochByIx(n)*charT3 + pBCSys->getEpoch0() - SysData_bc4bp::REF_EPOCH)/charT2);
        node.setExtraParam("J", DynamicsModel_cr3bp::getJacobi(&(crNodeState[0]), pCRSys->getMu()));
        
        // Save the new ID
        map_oldID_to_newID[bcNodes.getNodeByIx(n).getID()] = crNodes.addNode(node);
    }

    for(int s = 0; s < bcNodes.getNumSegs(); s++){
        Segment seg = bcNodes.getSegByIx(s);
        seg.setTOF(seg.getTOF()*charT3/charT2);

        // Remap the origin and terminus to the new IDs
        if(seg.getOrigin() != Linkable::INVALID_ID)
            seg.setOrigin(map_oldID_to_newID[seg.getOrigin()]);

        if(seg.getTerminus() != Linkable::INVALID_ID)
            seg.setTerminus(map_oldID_to_newID[seg.getTerminus()]);

        crNodes.addSeg(seg);
    }

    return crNodes;
}//==================================================

/**
 *  \ingroup bc4bp
 *  \brief Transform a BCR4BPR Sun-Earth-Moon trajectory into a CR3BP Sun-Earth nodeset
 *
 *  \param bcTraj a BCR4BPR Sun-Earth-Moon trajectory
 *  \param pCRSys a CR3BP Sun-Earth system data object; contains information about system
 *  scaling
 *
 *  @return a CR3BP Trajectory object
 *  @throw Exception if <tt>pCRSys</tt> is not a Sun-Earth system or if <tt>bcTraj</tt> is
 *  not in the Sun-Earth-Moon system.
 */
Traj_cr3bp bcr4bpr_SEM2SE(Traj_bc4bp bcTraj, const SysData_cr3bp *pCRSys){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp seNodeset = bcr4bpr_SEM2SE(static_cast<Nodeset_bc4bp>(bcTraj), pCRSys);
    return static_cast<Traj_cr3bp>(seNodeset);
}//==================================================

/**
 *  \ingroup bc4bp
 *  \brief Compute the location of the saddle point for a specific bicircular system
 *  and epoch
 *  @details This function uses a Newton-Raphson procedure to locate the zeros of the local
 *  acceleration field.
 *  \throws DivergeException if the Newton-Raphson procedure cannot converge
 *  \param pBCSys system data object describing the bicircular system
 *  \param t0 the epoch at which the saddle point's location is computed
 * 
 *  @return the location of the saddle pointin BCR4BP rotating coordinates
 *  @throw DivergeException if the multiple shooting algorithm cannot converge on the
 *  saddle point location
 */
Eigen::Vector3d bcr4bpr_getSPLoc(const SysData_bc4bp *pBCSys, double t0){

    // Compute approximate SP location using 3-body dynamics
    SysData_cr3bp subSys(pBCSys->getPrimary(0), pBCSys->getPrimary(1));
    double a = 2*subSys.getMu() - 1;
    double b = 4*subSys.getMu()*subSys.getMu() - 4*subSys.getMu() + 2;
    double c = 2*pow(subSys.getMu(), 3) - 3*subSys.getMu()*subSys.getMu() + 3*subSys.getMu() - 1;
    
    // x-coordinate in bcr4bpr
    Eigen::Vector3d spPos((-b + sqrt(b*b - 4*a*c))/(2*a*pBCSys->getK()) - (1/pBCSys->getK() - pBCSys->getMu()), 0, 0);
    
    // Get primary positions
    double primPos[9];
    DynamicsModel_bc4bp::getPrimaryPos(t0, pBCSys, primPos);
    
    double err = 1;
    double okErr = 1e-10;
    int maxIts = 20;
    int count = 0;
    while(count < maxIts && err > okErr){
        
        Eigen::Vector3d s(spPos(0) - primPos[0], spPos(1) - primPos[1], spPos(2) - primPos[2]);
        Eigen::Vector3d e(spPos(0) - primPos[3], spPos(1) - primPos[4], spPos(2) - primPos[5]);
        Eigen::Vector3d m(spPos(0) - primPos[6], spPos(1) - primPos[7], spPos(2) - primPos[8]);
        double ds = s.norm();
        double de = e.norm();
        double dm = m.norm();
        double k = pBCSys->getK();
        double mu = pBCSys->getMu();
        double nu = pBCSys->getNu();

        // Compute acceleration due to three primaries
        Eigen::Vector3d accel = -(1/k - mu)*s/pow(ds,3) - (mu - nu)*e/pow(de,3) - nu*m/pow(dm,3);

        err = accel.norm();
        // astrohelion::printColor(YELLOW, "Iteration %02d: ||A|| = %.6e\n", count, err);
        if(err > okErr){
            // Partial derivatives; Create Jacobian
            Matrix3Rd J;
            J(0,0) = -(1/k - mu)*(1/pow(ds,3) - 3*s(0)*s(0)/pow(ds,5)) -
                (mu-nu)*(1/pow(de,3) - 3*e(0)*e(0)/pow(de,5)) - nu*(1/pow(dm,3) - 3*m(0)*m(0)/pow(dm,5));
            J(0,1) = (1/k - mu)*3*s(0)*s(1)/pow(ds,5) + (mu - nu)*3*e(0)*e(1)/pow(de,5) +
                nu*3*m(0)*m(1)/pow(dm,5);
            J(0,2) = (1/k - mu)*3*s(0)*s(2)/pow(ds,5) + (mu - nu)*3*e(0)*e(2)/pow(de,5) +
                nu*3*m(0)*m(2)/pow(dm,5);
            J(1,1) = -(1/k - mu)*(1/pow(ds,3) - 3*s(1)*s(1)/pow(ds,5)) -
                (mu-nu)*(1/pow(de,3) - 3*e(1)*e(1)/pow(de,5)) - nu*(1/pow(dm,3) - 3*m(1)*m(1)/pow(dm,5));
            J(1,2) = (1/k - mu)*3*s(1)*s(2)/pow(ds,5) + (mu - nu)*3*e(1)*e(2)/pow(de,5) +
                nu*3*m(1)*m(2)/pow(dm,5);
            J(2,2) = -(1/k - mu)*(1/pow(ds,3) - 3*s(2)*s(2)/pow(ds,5)) -
                (mu-nu)*(1/pow(de,3) - 3*e(2)*e(2)/pow(de,5)) - nu*(1/pow(dm,3) - 3*m(2)*m(2)/pow(dm,5));
            J(1,0) = J(0,1);
            J(2,0) = J(0,2);
            J(2,1) = J(1,2);

            // Solve linear system of equations using LU Decomposition
            accel *= -1;
            Eigen::FullPivLU<Matrix3Rd> lu(J);
            Eigen::Vector3d posDiff = lu.solve(accel);
            spPos += posDiff;
        }
        count++;
    }
    
    if(err <= okErr)
        return spPos;
    else
        throw DivergeException("Calculations::bcr4bpr_getSPLoc: Could not converge on SP location");
}//=====================================================

/**
 *  \ingroup bc4bp
 *  \brief Compute coefficients for 2nd-order polynomials in Epoch time that
 *  describe the x, y, and z coordinates of the saddle point.
 *  @details This function employs least squares to compute the coefficients. The
 *  number of points used and time span searched are hard coded in the function.
 * 
 *  \param pBCSys data about the bicircular system
 *  \param T0 the "center" epoch; points are generated within +/- timeSpan of this
 *  epoch
 * 
 *  @return a 3x3 matrix of coefficients. The first column contains the coefficients
 *  for the x-position approximation, the second column contains the y-position 
 *  coefficients, etc. The first row holds the second-order coefficients, 
 *  the second row holds the first-order coefficients, and the final row holds the
 *  0th order coefficients. To compute the saddle point's location at an epoch T,
 *  you can solve for x, y, and z with the following equation:
 *  
 *      [x, y, z] = [T*T, T, 1]*[Coefficient Matrix]
 */
MatrixXRd bcr4bpr_spLoc_polyFit(const SysData_bc4bp *pBCSys, double T0){
    const int STATE_SIZE = 3;   // predicting three position states
    double tSpan = 24*3600/pBCSys->getCharT();   // Amount of time to span before and after T0
    int numPts = 100;
    double dT = 2*tSpan/(static_cast<double>(numPts-1));

    MatrixXRd indVarMat(numPts, 3); // i.e. Vandermonde matrix
    MatrixXRd depVarMat(numPts, STATE_SIZE);
    for(int row = 0; row < 100; row++){
        double T = -tSpan + row*dT;     // Let T0 = 0 to center data and hopefully avoid numerical issues
        Eigen::Vector3d spPos = bcr4bpr_getSPLoc(pBCSys, T + T0);
        indVarMat(row, 0) = T*T;
        indVarMat(row, 1) = T;
        indVarMat(row, 2) = 1;
        depVarMat.block(row, 0, 1, 3) = spPos.transpose();
        // row++;
    }

    // Create Gramm matrix
    MatrixXRd G(indVarMat.cols(), indVarMat.cols());
    G.noalias() = indVarMat.transpose()*indVarMat;

    // Use 2nd-order polynomial fit; solve GC = A.transpose()*B for C
    MatrixXRd C = G.fullPivLu().solve(indVarMat.transpose()*depVarMat);
    
    // C holds the coefficients for the centered data set, need to adjust 
    // these coefficients for non-centered data (see notebook notes from 2/15)
    Eigen::RowVector3d linTerms = C.row(1) - 2*T0*C.row(0);
    Eigen::RowVector3d shiftTerms = C.row(0)*T0*T0 - C.row(1)*T0 + C.row(2);
    C.block(1,0,1,3) = linTerms;
    C.block(2,0,1,3) = shiftTerms;
    
    return C;
}//=================================================

/** \} */   // END of util group

}// END of Astrohelion namespace



