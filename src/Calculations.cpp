/**
 *  @file Calculations.cpp
 *
 *  @brief  Contains non-member calculation functions
 *
 *   This library contains functions to perform miscellaneous calculations.
 *   Common calculations, such as EOM evaluation or primary location, should
 *   be included in the tpat_model class and its derivatives.
 *   
 *   @author Andrew Cox
 *  @version May 25, 2016
 *  @copyright GNU GPL v3.0
 */
/*
 *  Astrohelion 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
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
#include "Constants.hpp"
#include "CorrectionEngine.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "Nodeset_bc4bp.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_ltvp.hpp"
#include "Traj.hpp"
#include "Traj_bc4bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Utilities.hpp"

#include "cspice/SpiceUsr.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <typeinfo>


namespace astrohelion{
//-----------------------------------------------------
//      General Utility Functions
//-----------------------------------------------------

/**
 *  @brief convert a date to epoch time
 *  
 *  @param date a string representing the date. The string can be formatted in one
 *  of two ways. First, a Gregorian-style date: 'YYYY/MM/DD HH:II:SS' (UTC, 24-hour clock);
 *  The time 'HH:II:SS' can be ommited, and the function will assume the time is 0:0:00
 *  Second, a Julian date (UTC) can be input with the format 'jd #' where '#' represents
 *  the Julian date.
 *
 *  @return the J2000 epoch time, or number of seconds after Jan 1, 2000 at 0:0:00 UTC.
 *  @throws Exception if the SPICE kernels cannot be loaded: the kernel names and
 *  filepaths are located in the settings XML file
 */
double dateToEpochTime(const char *date){
    std::string spice_path = Core::initializer.settings.spice_data_filepath;
    std::string time_kernel = Core::initializer.settings.spice_time_kernel;

    char timeKernel[512];
    sprintf(timeKernel, "%s%s", spice_path.c_str(), time_kernel.c_str());

    furnsh_c(timeKernel);

    if(failed_c()){
        char errMsg[26];
        getmsg_c("short", 25, errMsg);
        astrohelion::printErr("Spice Error: %s\n", errMsg);
        reset_c();  // reset error status
        throw Exception("Could not load SPICE time kernel");
    }

    // Convert the date to epoch time
    double et = 0;
    str2et_c(date, &et);

    // Unload the kernel
    unload_c(timeKernel);
    if(failed_c()){
        char errMsg[26];
        getmsg_c("short", 25, errMsg);
        astrohelion::printErr("Spice Error: %s\n", errMsg);
        reset_c();  // reset error status
        throw Exception("Could not unload SPICE time kernel");
    }

    return et;
}//==========================================

/**
 *  @brief use least squares to predict new values of variables in a continuation process
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
 *  @param indVarIx the index of the state variable to use as the indpendent variable
 *  @param nextInd The next value of the independent variable
 *  @param depVars a vector specifying the indices of the states that will be
 *  dependent variables. The algorithm will predict fugure values for these
 *  variables based on how they have changed with the independent variable.
 *  @param varHistory a vector representing an n x 8 matrix which contains
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
 *  @throws Exception if the <tt>varHistory</tt> vector contains fewer than 
 *  8 elements or if the <tt>depVars</tt> vector has no data
 */
std::vector<double> familyCont_LS(int indVarIx, double nextInd, std::vector<int> depVars, std::vector<double> varHistory){
    const int STATE_SIZE = 8;
    const double EPS = 1e-14;

    if(varHistory.size() < STATE_SIZE)
        throw Exception("tpat_calculations::familyCont_LS: Not enough data to create A matrix\n");

    if(depVars.size() == 0)
        throw Exception("tpat_calculations::familyCont_LS: Not enough data to create B matrix\n");

    // Form A and B matrices
    std::vector<double> A_data;
    std::vector<double> B_data;

    for(size_t n = 0; n < varHistory.size()/STATE_SIZE; n++){
        double d = varHistory[n*STATE_SIZE + indVarIx];
        A_data.push_back(d*d);  // ind. var^2
        A_data.push_back(d);    // ind. var^1
        A_data.push_back(1);    // ind. var^0

        for(size_t p = 0; p < depVars.size(); p++){
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
        for(size_t n = 0; n < varHistory.size()/STATE_SIZE; n++){
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
    for(size_t i = 0; i < depVars.size(); i++)
        predicted[depVars[i]] = P(i);

    return predicted;
}//====================================================

/**
 *  @brief construct a matrix to mirror a 6-d state over the specified plane or axis
 *  @param mirrorType describes how to mirror a 6-d state
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
            astrohelion::printErr("tpat_calculations::getMirrorMat: Mirror type is not implemented; returning identiy\n");
            return Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Identity();
    }
}//===================================================

/**
 *  @brief Sort a list of eigenvalues
 *
 *  This algorithm is intended to sort eigenvalues from the monodromy matrices of a family of
 *  periodic orbits. Because of this, the code attempts to locate at least one pair of 
 *  eigenvalues that are equal to 1.0. Other eigenvalues come in either real, reciprocal pairs
 *  or in complex conjugate pairs. The first row is sorted such that the pairs of eigenvalues
 *  are located next to each other in the row. All subsequent rows are sorted so that the order
 *  of the eigenvalues remains the same
 *
 *
 *  @param eigVals vector of eigenvalues in row-major order. We assume that 6 eigenvalues are 
 *  included in each row.
 *  @param sortedIxs a pointer to an integer vector that will store the indices of the original
 *  eigenvalues in their new order (after sorting)
 *  @return the sorted eigenvalues, again in row-major order in a vector
 */
std::vector<cdouble> sortEig(std::vector<cdouble> eigVals, std::vector<int> *sortedIxs){
    if(eigVals.size() == 0){
        astrohelion::printErr("tpat_calculations::sortEig: Cannot sort eigenvalues: there are no family members\n");
        return eigVals;
    }

    if(eigVals.size() % 6 != 0){
        astrohelion::printErr("tpat_calculations::sortEig: Must have 6n eigenvalues\n");
        return eigVals;
    }

    const double MAX_ONES_ERR = 5e-5;
    std::vector<cdouble> sortedEigs;

    // Figure out which eigenvalues are closest to one, save their indices
    std::vector<double> onesErr;
    std::vector<cdouble> e1(eigVals.begin(), eigVals.begin()+6);
    for(int i = 0; i < 6; i++){
        onesErr.push_back(std::abs(std::real(e1[i])-1) + std::abs(std::imag(e1[i])));
    }
    
    std::vector<double>::iterator smallestErr = std::min_element(onesErr.begin(), onesErr.end());

    int smallIx = smallestErr - onesErr.begin();
    double saveSmallestErrVal = *smallestErr;
    // artifically inflate the error on the smallest so we can find the second smallest
    onesErr[smallIx] += 1e20;
    std::vector<double>::iterator secondSmallestErr = std::min_element(onesErr.begin(), onesErr.end());

    int otherSmallIx = secondSmallestErr - onesErr.begin();

    int onesIx[] = {smallIx, otherSmallIx};   // Indices of the eigenvalues that (nearly) exactly 1.0
    if(saveSmallestErrVal > MAX_ONES_ERR || *secondSmallestErr > MAX_ONES_ERR){
        astrohelion::printWarn("tpat_calculations::sortEigs: Eigenvalues closest to one have errors of %.4e and %.4e > %.4e\n",
            saveSmallestErrVal, *secondSmallestErr, MAX_ONES_ERR);
    }

    // Generate all permutations of the indices 0 through 5
    std::vector<int> vals {0,1,2,3,4,5};
    std::vector<int> ixPerms = astrohelion::generatePerms<int>(vals);
    cdouble predict[6];
    std::copy(e1.begin(), e1.end(), predict);

    for(size_t m = 0; m < eigVals.size()/6; m++){
        if(m == 1){
            if(saveSmallestErrVal < MAX_ONES_ERR){
                // Update the indices of the ones in case they got moved in the first sorting
                onesIx[0] = sortedIxs->at(onesIx[0]);
                onesIx[1] = sortedIxs->at(onesIx[1]);
                printf("Updated unit eigenvalue positions to %d and %d\n", onesIx[0], onesIx[1]);
            }
            std::copy(sortedEigs.begin(), sortedEigs.begin()+6, predict);
        }else if(m > 1){
            // Use linear extrapolation of m-1 and m-2 to predict m
            for(int i = 0; i < 6; i++){
                cdouble deltaEig = sortedEigs[(m-1)*6 + i] - sortedEigs[(m-2)*6 + i];
                predict[i] = sortedEigs[(m-1)*6 + i] + deltaEig;
            }
        }

        // astrohelion::printColor(RED, "Unsorted Set %03d: [%s %s %s %s %s %s]\n", m, complexToStr(eigVals[m*6+0]).c_str(),
        //     complexToStr(eigVals[m*6+1]).c_str(), complexToStr(eigVals[m*6+2]).c_str(), complexToStr(eigVals[m*6+3]).c_str(),
        //     complexToStr(eigVals[m*6+4]).c_str(), complexToStr(eigVals[m*6+5]).c_str());
        
        std::vector<int> killers(ixPerms.size()/6, 0);  // Keep track of which cost killed each permutation possibility
        std::vector<double> cost;
        cost.assign(ixPerms.size()/6, 0);
        cdouble one(1,0);
        for(size_t p = 0; p < ixPerms.size()/6; p++){

            // rearrange the splitOrig vector using the current permutation
            cdouble swappedOrig[6];
            for(int i = 0; i < 6; i++){
                swappedOrig[i] = eigVals[m*6 + ixPerms[6*p + i]];
            }

            // Compute difference between this permutation and prediction
            for(int i = 0; i < 6; i++){
                cost[p] += std::abs(swappedOrig[i] - predict[i]);
            }
            
            // Make sure ones end up next to each other the first time around
            bool disqual = false;

            if(!disqual){
                if(m == 0){
                    // Determine the new indices of the eigenvalues identified as ones
                    std::vector<int>::iterator a = std::find(ixPerms.begin()+6*p, ixPerms.begin()+6*(p+1), onesIx[0]);
                    std::vector<int>::iterator b = std::find(ixPerms.begin()+6*p, ixPerms.begin()+6*(p+1), onesIx[1]);
                    int ix0 = a - (ixPerms.begin()+6*p);
                    int ix1 = b - (ixPerms.begin()+6*p);

                    // If the two unit eigenvalues are not next to each other, toss this option
                    if(std::abs(ix0 - ix1) > 1){
                        cost[p] += 1e20;
                        killers[p] = 1;
                        disqual = true; // This permuation has been disqualified
                    }else{
                        // astrohelion::printColor(GREEN, "Permutation %03d has ones in indices %d and %d\n", p, ix0, ix1);
                    }
                }
            }

            if(!disqual){
                // Add infinite cost if the ones eigenvalues change spots
                if(saveSmallestErrVal < MAX_ONES_ERR){
                    double distToOne[6];
                    for(int i = 0; i < 6; i++){
                        distToOne[i] = std::abs(swappedOrig[i]-one);
                    }
                    double *closestToOne = std::min_element(distToOne, distToOne+6);
                    int closestIx = closestToOne - distToOne;

                    if( (closestIx != onesIx[0] && closestIx != onesIx[1]) ){
                        cost[p] += 1e20;
                        killers[p] = 2;
                        disqual = true;
                    }                    
                    // distToOne[closestIx] += 1e20;
                    // double *secondClosest = std::min_element(distToOne, distToOne+6);
                    // int secClosestIx = secondClosest - distToOne;

                    // if( (closestIx != onesIx[0] && closestIx != onesIx[1]) ||
                    //     (secClosestIx != onesIx[0] && secClosestIx != onesIx[1]) ){
                    //     cost[p] += 1e20;
                    //     killers[p] = 2;
                    //     disqual = true;
                    // }
                }
            }

            if(!disqual){
                // Add infinite cost if each consecutive pair are not inverses
                // or complex conjugates
                for(int i = 0; i < 6; i+=2){
                    /* don't consider the eigenvalues that are closest to one
                        because they may have small innaccuracies that trigger
                        these costs, screwing up the algorithm */
                    if(i != onesIx[0]){
                        // Complex parts don't sum to zero (conjugates will, two
                        // real eigenvalues will)
                        if(std::abs(std::imag(swappedOrig[i]) + std::imag(swappedOrig[i+1]))){
                            cost[p] += 1e20;
                            killers[p] = 3;
                            disqual = true;
                            break;
                        }

                        // We have established that the ones are in the same place; Compute
                        // the error in a 1.0 eigenvalue; this error is acceptable in the
                        // reciprocal calculation (error tends to get bad at the larger orbits)
                        double okErr = std::abs(swappedOrig[onesIx[0]] - one);

                        // Both are real but are not reciprocals
                        if(std::imag(swappedOrig[i]) == 0 && std::imag(swappedOrig[i+1]) == 0){
                            double recipErr = 1.0 - std::real(swappedOrig[i])*std::real(swappedOrig[i+1]);
                            if(recipErr > okErr){
                                cost[p] += 1e20;
                                killers[p] = 4;
                                disqual = true;
                                break;
                            }
                        }
                    }
                }// end of checking consecutive pairs
            }
        }// end of loop through permutations

        // Find the minimum cost
        std::vector<double>::iterator minCost = std::min_element(cost.begin(), cost.end());
        int ix = minCost - cost.begin();
        for(int i = 0; i < 6; i++){
            sortedEigs.push_back(eigVals[m*6 + ixPerms[ix*6 + i]]);
        }

        if(*minCost > 1e10){
            astrohelion::printErr("Minimum cost on set #%03d is %e - perm %03d, killed by cost %d - probably a bug!\n", m, ix, killers[ix], *minCost);
        }

        // printf("Chose Perm %04d: [%d %d %d %d %d %d]\n", ix, ixPerms[6*ix+0], ixPerms[6*ix+1], ixPerms[6*ix+2],
        //     ixPerms[6*ix+3], ixPerms[6*ix+4], ixPerms[6*ix+5]);

        sortedIxs->insert(sortedIxs->end(), ixPerms.begin() + ix*6, ixPerms.begin() + (ix+1)*6);

        // astrohelion::printColor(GREEN, "  Sorted Set %03d: [%s %s %s %s %s %s]\n", m, complexToStr(sortedEigs[m*6+0]).c_str(),
        //     complexToStr(sortedEigs[m*6+1]).c_str(), complexToStr(sortedEigs[m*6+2]).c_str(), complexToStr(sortedEigs[m*6+3]).c_str(),
        //     complexToStr(sortedEigs[m*6+4]).c_str(), complexToStr(sortedEigs[m*6+5]).c_str());
        // waitForUser();
    }// end of loop through all members/eigenvalue sets

    return sortedEigs;
}//=====================================================

/**
 *  @brief Compute manifolds of a CR3BP trajectory
 *  @details [long description]
 * 
 *  @param type The type of manifolds to generate
 *  @param perOrbit A periodic, CR3BP orbit. No checks are made to ensure periodicity,
 *  so this function also performs well for nearly periodic segments of quasi-periodic
 *  arcs. If an arc that is not approximately periodic is input, the behavior may be... strange.
 *  @param numMans The number of manifolds to generate
 *  @param tof Time-of-flight along each manifold arc after stepping off the specified orbit
 *  @return a vector of trajectory objects, one for each manifold arc.
 *  @throws Exception if the eigenvalues cannot be computed, or if only one
 *  stable or unstable eigenvalue is computed (likely because of an impropper monodromy matrix)
 */
std::vector<Traj_cr3bp> getManifolds(Manifold_tp type, const Traj_cr3bp *perOrbit, int numMans, double tof){
    // Get eigenvalues of monodromy matrix
    MatrixXRd mono = perOrbit->getSTMByIx(-1);

    Eigen::EigenSolver<MatrixXRd> eigensolver(mono);
    if(eigensolver.info() != Eigen::Success)
        throw Exception("tpat_calculations::getManifolds: Could not compute eigenvalues of monodromy matrix");

    Eigen::VectorXcd vals = eigensolver.eigenvalues();
    MatrixXRcd eigVecs = eigensolver.eigenvectors();
    std::vector<cdouble> eigData(vals.data(), vals.data()+6);

    // Sort eigenvalues to put them in a "propper" order, get indices
    // to sort the eigenvectors to match
    std::vector<int> sortedIx;
    std::vector<cdouble> sortedEig = sortEig(eigData, &sortedIx);

    // Figure out which eigenvalues are the stable and unstable ones,
    // delete the rest
    std::vector<cdouble> nonCenterVals;
    std::vector<cdouble> nonCenterVecs;
    for(size_t c = 0; c < sortedEig.size(); c++){
        double realErr = std::real(sortedEig[c]) - 1.0;
        double imagErr = std::imag(sortedEig[c]);

        if(std::abs(realErr) > 1e-5 && std::abs(imagErr) < 1e-5){
            // Keep this eigenvalue/eigenvector pair
            nonCenterVals.push_back(sortedEig[c]);
            int vecIx = sortedIx[c];
            nonCenterVecs.insert(nonCenterVecs.end(),
                eigVecs.data()+vecIx*6, eigVecs.data()+(vecIx+1)*6); 
        }
    }

    std::vector<Traj_cr3bp> allManifolds;

    if(nonCenterVals.size() == 0){
        astrohelion::printWarn("tpat_calculations::getManifolds: No stable/unstable eigenvalues were found\n");
        return allManifolds;
    }

    if(nonCenterVals.size() == 1){
        throw Exception("tpat_calculations::getManifolds: Only found one stable/unstable eigenvalue");
    }

    if(nonCenterVals.size() > 2){
        astrohelion::printWarn("tpat_calculations::getManifolds: Stable/Unstable subspace is larger than 2D. Only the first pair will be considered\n");
    }

    /** TODO: If and when I can flexibily integrate generic trajectories and nodesets,
     *  I would like to replace this with a nodeset to space the points evenly
     *  in time and/or arclength
     */
    // Get a bunch of points to use as starting guesses for the manifolds
    if(numMans > perOrbit->getNumNodes()){
        astrohelion::printWarn("tpat_calculations::getManifolds: Requested too many manifolds... will return fewer\n");
        numMans = perOrbit->getNumNodes();
    }

    double stepSize = ((double)perOrbit->getNumNodes())/((double)numMans);
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
    double charL = perOrbit->getSysData()->getCharL();
    for(int n = 0; n < numMans; n++){
        // Transform the eigenvectors to this updated time
        MatrixXRd newVecs = perOrbit->getSTMByIx(pointIx[n])*vecs;

        // Pick the direction from one of the transformed eigenvectors
        Eigen::VectorXd direction(6);
        for(size_t v = 0; v < 2; v++){
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
        std::vector<double> state = perOrbit->getStateByIx(pointIx[n]);
        Eigen::VectorXd q0 = Eigen::Map<Eigen::VectorXd>(&(state[0]), 6, 1);
        q0 += stepDist/charL * direction;

        // Simulate for some time to generate a manifold arc
        Traj_cr3bp traj(static_cast<const SysData_cr3bp *>(perOrbit->getSysData()));
        sim.runSim(q0.data(), tof, &traj);
        allManifolds.push_back(traj);
    }

    return allManifolds;
}//====================================================

/**
 *  @brief Compute the stability index of a periodic orbit from a set of eigenvalues
 *  @details This algorithm assumes the orbit is periodic and that the eigenvalues 
 *  have been sorted (so they come in pairs).
 * 
 *  @param eigs A 6-element vector of eigenvalues associated with a periodic orbit
 *  @return the stability index, or NAN if no real, reciprocal eigenvalue pair is found
 *  @throws Exception if <tt>eigs</tt> does not have six elements
 */
double getStabilityIndex(std::vector<cdouble> eigs){
    if(eigs.size() != 6)
        throw Exception("tpat_calculations::getStabilityIndex: Must input 6 eigenvalues!");

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
 *  @brief Compute the total delta-V along a corrected nodeset
 * 
 *  @param it pointer to an MultShootData object associated with a corrections process
 *  @return the total delta-V, units consistent with the nodeset's stored velocity states
 */
double getTotalDV(const MultShootData *it){
    double total = 0;
    for(size_t n = 0; n < it->deltaVs.size()/3; n++){
        total += sqrt(it->deltaVs[3*n + 0]*it->deltaVs[3*n + 0] +
            it->deltaVs[3*n + 1]*it->deltaVs[3*n + 1] + 
            it->deltaVs[3*n + 2]*it->deltaVs[3*n + 2]);
    }
    return total;
}//=====================================================

/**
 *  @brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  @details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  @param nodeset A nodeset with some constraints
 */
void finiteDiff_checkMultShoot(const Nodeset *nodeset){
    CorrectionEngine engine;  // Create engine with default settings
    finiteDiff_checkMultShoot(nodeset, engine);
}//====================================================

/**
 *  @brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  @details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  @param nodeset A nodeset with some constraints
 *  @param engine correction engine object configured with the appropriate settings (i.e.,
 *  equal arc time, scaling variables, etc.). Note that the maxIts, verbosity, and ignoreDiverge
 *  attributes of the engine will be overridden by this function.
 */
void finiteDiff_checkMultShoot(const Nodeset *nodeset, CorrectionEngine engine){
    printf("Finite Diff: Checking DF matrix... ");
    // Create multiple shooter that will only do 1 iteration
    CorrectionEngine corrector(engine);
    corrector.setMaxIts(1);
    corrector.setVerbose(Verbosity_tp::NO_MSG);
    corrector.setIgnoreDiverge(true);

    // Run multiple shooter to get X, FX, and DF
    MultShootData it = corrector.multShoot(nodeset, NULL);
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

    astrohelion::toCSV(DF, "FiniteDiff_DF.csv");
    astrohelion::toCSV(DFest, "FiniteDiff_DFest.csv");
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
    // astrohelion::toCSV(diff, "FiniteDiff_Diff.csv");

    Eigen::VectorXd rowMax = diff.rowwise().maxCoeff();
    Eigen::RowVectorXd colMax = diff.colwise().maxCoeff();

    double rowMaxMax = rowMax.maxCoeff();
    double colMaxMax = colMax.maxCoeff();
    int errScalar = 10000;

    if(rowMaxMax < errScalar*pertSize && colMaxMax < errScalar*colMaxMax){
        astrohelion::printColor(BOLDGREEN, "No significant errors!\n");
    }else{
        astrohelion::printColor(BOLDRED, "Significant errors!\n");
        printf("Maximum relative difference between computed DF and estimated DF\n");
        int conCount = 0;
        for(long r = 0; r < rowMax.size(); r++){
            if(r == 0 && it.totalCons > 0){
                printf("Node %d %s Constraint:\n", it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
            }else if(conCount < (int)(it.allCons.size()) && r >= it.conRows[conCount+1]){
                conCount++;
                printf("Node %d %s Constraint:\n", it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
            }
            astrohelion::printColor(rowMax[r] > errScalar*pertSize || std::isnan(rowMax[r]) ? RED : GREEN, "  row %03zu: %.6e\n", r, rowMax[r]);
        }
        for(long c = 0; c < colMax.size(); c++){
            astrohelion::printColor(colMax[c] > errScalar*pertSize || std::isnan(colMax[c]) ? RED : GREEN, "Free Var %03zu: %.6e\n", c, colMax[c]);
        }
    }
}//====================================================

//-----------------------------------------------------
//      CR3BP Utility Functions
//-----------------------------------------------------

/**
 *  @brief Compute the magnitude of a velocity component given Jacobi Constant
 * 
 *  @param s state vector (non-dimensional); MUST contain at least the 6 position and velocity states.
 *  Note that the desired velocity component identified by velIxToFind is not required; put a placeholder
 *  zero or NAN (or anything, really) in its place; this value will not be used in computations.
 *  @param mu non-dimensional system mass ratio
 *  @param C Jacobi constant value
 *  @param velIxToFind index of the velocity component to compute (i.e. 3, 4, or 5)
 *  @return the magnitude of vx (3), vy (4), or vz (5).
 *  @throws Exception if <tt>velIxToFind</tt> is out of bounds
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
 *  @brief Compute a periodic orbit in the CR3BP system
 *
 *  The initial and final states are constrained based on the mirror type, but
 *  all other states are allowed to vary during the corrections process. If you
 *  wish to constrain a specific states, see 
 *  cr3bp_getPeriodic(sys, IC, period, numNodes, mirrorType, fixedStates)
 *  This function also uses only two nodes; to specify more, see the function above.
 *
 *  This function also assumes the order of the periodic orbit is 1
 *
 *  @param sys the dynamical system
 *  @param IC non-dimensional initial state vector
 *  @param period non-dimensional period for the orbit
 *  @param mirrorType how this periodic orbit mirrors in the CR3BP
 *  @param tol tolerance to use in the corrections process
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 */
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp *sys, std::vector<double> IC,
    double period, Mirror_tp mirrorType, double tol){
    
    std::vector<int> fixedStates;   // Initialize an empty vector
    return cr3bp_getPeriodic(sys, IC, period, 2, 1, mirrorType, fixedStates, tol);
}//========================================

/**
 *  @brief Compute a periodic orbit in the CR3BP system
 *  @details This method ignores all crash events, so it is possible to compute a
 *  periodic orbit that passes through a primary
 *  
 *  @param sys the dynamical system
 *  @param IC non-dimensional initial state vector
 *  @param period non-dimensional period for the orbit
 *  @param numNodes the number of nodes to use for HALF of the periodic orbit; more nodes 
 *  may result in a more robust correction
 *  @param order the number of revolutions about the system/primary this orbit completes before
 *  it repeats periodically. Think of a period-3 DRO (order = 3) or a butterfly (order = 2)
 *  @param mirrorType how this periodic orbit mirrors in the CR3BP
 *  @param fixedStates a vector containing the indices of which initial states
 *  we would like to fix. Not all states are possible for each mirror condition.
 *  See the enum definition for specific details.
 *  @param tol tolerance to use in the corrections process
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 */
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp *sys, std::vector<double> IC,
    double period, int numNodes, int order, Mirror_tp mirrorType, std::vector<int> fixedStates,
    double tol){

    // MultShootData itData;
    // return cr3bp_getPeriodic(sys, IC, period, numNodes, order, mirrorType, fixedStates, tol, &itData);
    return cr3bp_getPeriodic(sys, IC, period, numNodes, order, mirrorType, fixedStates, tol, NULL);
}//====================================================================

/**
 *  @brief Compute a periodic orbit in the CR3BP system
 *  @details This method ignores all crash events, so it is possible to compute a
 *  periodic orbit that passes through a primary
 *  
 *  @param sys the dynamical system
 *  @param IC non-dimensional initial state vector
 *  @param period non-dimensional period for the orbit
 *  @param numNodes the number of nodes to use for HALF of the periodic orbit; more nodes 
 *  may result in a more robust correction
 *  @param order the number of revolutions about the system/primary this orbit completes before
 *  it repeats periodically. Think of a period-3 DRO (order = 3) or a butterfly (order = 2)
 *  @param mirrorType how this periodic orbit mirrors in the CR3BP
 *  @param fixedStates a vector containing the indices of which initial states
 *  we would like to fix. Not all states are possible for each mirror condition.
 *  See the enum definition for specific details.
 *  @param tol tolerance to use in the corrections process
 *  @param itData a pointer to an iteration data object that contains data from the
 *  multiple shooting run that corrects the periodic orbit; this is useful when
 *  attempting to determine how well (or poorly) the multiple shooting algorithm
 *  performed (e.g., for a variable step-size process)
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 *  @throws Exception if <tt>mirrorType</tt> is invalid
 *  @throws DivergeException if the multiple shooting algorithm cannot converge on a 
 *  mirrored solution.
 */
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp *sys, std::vector<double> IC,
    double period, int numNodes, int order, Mirror_tp mirrorType, std::vector<int> fixedStates,
    double tol, MultShootData* itData){

    SimEngine sim;    // Engine to perform simulation
    sim.setAbsTol(tol < 1e-12 ? 1e-15 : tol/1000.0);
    sim.setRelTol(sim.getAbsTol());
    sim.clearEvents();                  // Ignore any crashes into the primaries
    std::vector<int> zeroStates;        // Which states must be zero to ensure a perpendicular crossing

    Event mirrorEvt(sys);
    // Determine which states must be zero for mirroring
    switch(mirrorType){
        case Mirror_tp::MIRROR_XZ:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(3);    // x-dot
            zeroStates.push_back(5);    // z-dot
            mirrorEvt.createEvent(Event_Tp::XZ_PLANE, 0, true);     // Tell the sim to quit once it reaches the XZ plane
            break;
        case Mirror_tp::MIRROR_YZ:
            zeroStates.push_back(0);    // x
            zeroStates.push_back(4);    // y-dot
            zeroStates.push_back(5);    // z-dot
            mirrorEvt.createEvent(Event_Tp::YZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_XY:
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            zeroStates.push_back(4);    // y-dot
            mirrorEvt.createEvent(Event_Tp::YZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_X_AX_H:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            mirrorEvt.createEvent(Event_Tp::XZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_X_AX_V:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            mirrorEvt.createEvent(Event_Tp::XY_PLANE, 0, true);
            break;
        default:
            throw Exception("Mirror type either not defined or not implemented");
    }
    // sim.setVerbose(true);
    mirrorEvt.setStopCount(order);
    sim.addEvent(mirrorEvt);

    // Create a constraint to enforce mirror condition at the beginning and end of arc
    double mirrorCon0[] = {NAN,NAN,NAN,NAN,NAN,NAN};
    double mirrorCon1[] = {NAN,NAN,NAN,NAN,NAN,NAN};

    // Set the zero states to zero in both constraints
    for(size_t i = 0; i < zeroStates.size(); i++){
        mirrorCon0[zeroStates[i]] = 0;
        mirrorCon1[zeroStates[i]] = 0;
    }

    // Fix states at the initial point
    for(size_t i = 0; i < fixedStates.size(); i++){
        bool okToFix = true;
        for(size_t n = 0; n < zeroStates.size(); n++){
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
    Traj_cr3bp halfOrbArc(sys);
    sim.runSim(IC, period, &halfOrbArc);
    
    // Check to make sure the simulation ended with the event (not running out of time)
    std::vector<Event> endEvts = sim.getEndEvents(&halfOrbArc);
    if(std::find(endEvts.begin(), endEvts.end(), mirrorEvt) == endEvts.end()){
        astrohelion::printErr("tpat_calculations::cr3bp_getPeriodic: simulation of half-period orbit did not end in mirror event; may have diverged\n");
    }
    // halfOrbArc.saveToMat("HalfOrbArc.mat");

    double halfOrbTOF = halfOrbArc.getTimeByIx(-1);
    double tofErr = 100*std::abs(halfOrbTOF-period/2.0)/(period/2.0);

    if(tofErr > 10)
        astrohelion::printWarn("tpat_calculations::cr3bp_getPeriodic: Half-Period arc TOF varies from input half-period by more than 10%%\n");

    // Create a nodeset from arc
    Nodeset_cr3bp halfOrbNodes(halfOrbArc, numNodes, Nodeset::DISTRO_TIME);
    halfOrbNodes.addConstraint(initStateCon);
    halfOrbNodes.addConstraint(finalStateCon);

    // Use differential corrections to enforce the mirror conditions
    CorrectionEngine corrector;
    corrector.setTol(tol);
    corrector.setIgnoreCrash(true); // Corrector also ignores crash events
    corrector.setVarTime(true);
    corrector.setEqualArcTime(true);

    try{
        Nodeset_cr3bp correctedHalfPer(sys);
        *itData = corrector.multShoot(&halfOrbNodes, &correctedHalfPer);

        // Make the nodeset into a trajectory
        Traj_cr3bp halfPerTraj = Traj_cr3bp::fromNodeset(correctedHalfPer);
        double halfTOF = halfPerTraj.getTimeByIx(-1);
        double halfPerTraj_len = halfPerTraj.getNumNodes();
        MatrixXRd halfPerSTM = halfPerTraj.getSTMByIx(-1);
        
        // Use Mirror theorem to create the second half of the orbit
        MatrixXRd mirrorMat = getMirrorMat(mirrorType);
        for(int i = halfPerTraj_len-2; i >= 0; i--){
            // Use mirroring to populate second half of the orbit
            std::vector<double> state = halfPerTraj.getStateByIx(i);
            Eigen::RowVectorXd stateVec = Eigen::Map<Eigen::RowVectorXd>(&(state[0]), 1, 6);
            Eigen::RowVectorXd newStateVec = stateVec*mirrorMat;

            Node node;
            node.setState(newStateVec.data());
            node.setEpoch(2*halfTOF - halfPerTraj.getTimeByIx(i));

            int id = halfPerTraj.addNode(node);
            Segment seg(id-1, id, halfPerTraj.getTimeByIx(id) - halfPerTraj.getTimeByIx(id-1));
            halfPerTraj.addSeg(seg);
        }

        // Compute the monodromy matrix from the half-period STM
        double M_data[] = {   0, 0, 0, -1, 0, 0,
                            0, 0, 0, 0, -1, 0,
                            0, 0, 0, 0, 0, -1,
                            1, 0, 0, 0, -2, 0,
                            0, 1, 0, 2, 0, 0,
                            0, 0, 1, 0, 0, 0};
        // Inverse of M
        double MI_data[] = {  0, -2, 0, 1, 0, 0,
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
        
        return halfPerTraj;     // Now contains entire trajectory
    }catch(DivergeException &e){
        throw DivergeException("tpat_calculations::cr3bp_getPeriodic: Could not converge half-period arc with mirroring condition");
    }
}//================================================

/**
 *  @brief Transition a trajectory from the EM system to the SE system
 *  
 *  The relative orientation between the two systems is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. These angles describe the
 *  system orientation at time t = 0, so if the start of the trajectory does not coincide
 *  with this initial time, be sure to shift the time coordinates to assure that the first
 *  value in the time vector reflects correct initial time.
 *
 *  @param EMTraj a CR3BP Earth-Moon trajectory
 *  @param SESys a Sun-Earth CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Traj_cr3bp cr3bp_EM2SE(Traj_cr3bp EMTraj, const SysData_cr3bp *SESys, double thetaE0, double thetaM0, double gamma){
    printf("cr3bp_EM2SE\n");
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp seNodes = cr3bp_EM2SE(static_cast<Nodeset_cr3bp>(EMTraj), SESys, thetaE0, thetaM0, gamma);
    return static_cast<Traj_cr3bp>(seNodes);
}//=========================================================

/**
 *  @brief Transition a nodeset from the EM system to the SE system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. Accordingly, the EM nodeset
 *  should have node epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire set may be accomplished via the updateEpochs() function.
 *
 *  @param EMNodes a CR3BP Earth-Moon nodeset
 *  @param SESys a Sun-Earth CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Nodeset_cr3bp cr3bp_EM2SE(Nodeset_cr3bp EMNodes, const SysData_cr3bp *SESys, double thetaE0, double thetaM0,
    double gamma){

    // astrohelion::printColor(BLUE, "Converting EM to SE\nEM Sys:\n  %d Nodes\n  %d Segments\n", EMNodes.getNumNodes(),
    //     EMNodes.getNumSegs());

    Nodeset_cr3bp SENodes(SESys);

    double charTE = EMNodes.getSysData()->getCharT();       // characteristic time in EM system
    double charLE = EMNodes.getSysData()->getCharL();       // characteristic length in EM system
    double charTS = SESys->getCharT();                       // characteristic time in SE system
    double charLS = SESys->getCharL();                       // characteristic length in SE system

    std::vector<int> map_oldID_to_newID(EMNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> state_SE;
    double epoch;

    for(int n = 0; n < EMNodes.getNumNodes(); n++){
        epoch = EMNodes.getEpochByIx(n);
        state_SE = cr3bp_EM2SE_state(EMNodes.getStateByIx(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, SESys->getMu());
        
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
 *  @brief Transition a trajectory from the SE system to the EM system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. Accordingly, the SE trajectory
 *  should have epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire trajectory may be accomplished via the updateEpochs() function.
 *
 *  @param SETraj a CR3BP Sun-Earth trajectory
 *  @param EMSys an Earth-Moon CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Traj_cr3bp cr3bp_SE2EM(Traj_cr3bp SETraj, const SysData_cr3bp *EMSys, double thetaE0, double thetaM0, double gamma){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp emNodes = cr3bp_SE2EM(static_cast<Nodeset_cr3bp>(SETraj), EMSys, thetaE0, thetaM0, gamma);
    return static_cast<Traj_cr3bp>(emNodes);
}//=========================================================

/**
 *  @brief Transition a nodeset from the SE system to the EM system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. Accordingly, the SE nodeset
 *  should have node epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire set may be accomplished via the updateEpochs() function.
 *
 *  @param SENodes a CR3BP Sun-Earth nodeset
 *  @param EMSys an Earth-Moon CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Nodeset_cr3bp cr3bp_SE2EM(Nodeset_cr3bp SENodes, const SysData_cr3bp *EMSys, double thetaE0, double thetaM0,
    double gamma){

    // astrohelion::printColor(BLUE, "Converting SE to EM\nSE Sys:\n  %d Nodes\n  %d Segments\n", SENodes.getNumNodes(),
    //     SENodes.getNumSegs());

    Nodeset_cr3bp EMNodes(EMSys);

    double charTE = EMSys->getCharT();                   // characteristic time in EM system
    double charLE = EMSys->getCharL();                   // characteristic length in EM system
    double charTS = SENodes.getSysData()->getCharT();    // characteristic time in SE system
    double charLS = SENodes.getSysData()->getCharL();    // characteristic length in SE system

    const SysData_cr3bp *SESys = static_cast<const SysData_cr3bp*>(SENodes.getSysData());

    std::vector<int> map_oldID_to_newID(SENodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> state_EM;
    double epoch;

    for(int n = 0; n < SENodes.getNumNodes(); n++){
        epoch = SENodes.getEpochByIx(n);
        // Transform a single node
        state_EM = cr3bp_SE2EM_state(SENodes.getStateByIx(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, SESys->getMu());
        
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
 *  @brief Transform a single state from EM coordinates to SE coordinates
 *
 *  @param state_EM a 6- or 9-element state vector
 *  @param t Earth-Moon nondimensional time associated with the Earth-Moon state
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 *  @param charLE EM characteristic length
 *  @param charTE EM characterstic time
 *  @param charLS SE characteristic length
 *  @param charTS SE characteristic time
 *  @param mu_SE SE mass ratio
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
 *  @brief Transform a single state from SE coordinates to EM coordinates
 *
 *  @param state_SE a 6- or 9-element state vector
 *  @param t Sun-Earth nondimensional time associated with the state
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 *  @param charLE EM characteristic length, km
 *  @param charTE EM characterstic time, sec
 *  @param charLS SE characteristic length, km
 *  @param charTS SE characteristic time, sec
 *  @param mu_SE SE mass ratio, nondimensional
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
 *  @brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  @param traj CR3BP trajectory
 *  @param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP trajectory is input, selecting 
 *  <tt>centerIx = 0</tt> will produce Earth-centered inertial coordinates and selecting
 *  <tt>centerIx = 1</tt> will produce Moon-centered inertial coordinates.
 * 
 *  @return A trajectory centered around the specified index in inertial, dimensional coordinates
 */
Traj_cr3bp cr3bp_rot2inert(Traj_cr3bp traj, int centerIx){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp inertNodes = cr3bp_rot2inert(static_cast<Nodeset_cr3bp>(traj), centerIx);
    return static_cast<Traj_cr3bp>(inertNodes);
}//==================================================================

/**
 *  @brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  @param nodes CR3BP nodeset
 *  @param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP nodeset is input, selecting 
 *  <tt>centerIx = 0</tt> will produce Earth-centered inertial coordinates and selecting
 *  <tt>centerIx = 1</tt> will produce Moon-centered inertial coordinates.
 * 
 *  @return A nodeset centered around the specified index in inertial, dimensional coordinates
 */
Nodeset_cr3bp cr3bp_rot2inert(Nodeset_cr3bp nodes, int centerIx){

    if(centerIx < 0 || centerIx > 1){
        throw Exception("tpat_calculations::cr3bp_rot2inert: Invalid center index");
    }

    const SysData_cr3bp *sys = static_cast<const SysData_cr3bp *>(nodes.getSysData());

    Nodeset_cr3bp inertNodes(sys);
    std::vector<int> map_oldID_to_newID(nodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> stateInert;
    
    for(int i = 0; i < nodes.getNumNodes(); i++){
        stateInert = cr3bp_rot2inert_state(nodes.getStateByIx(i), sys, nodes.getEpochByIx(i), centerIx);
        
        // Save the new ID    
        map_oldID_to_newID[nodes.getNodeByIx(i).getID()] = inertNodes.addNode(Node(stateInert, nodes.getEpochByIx(i)*sys->getCharT()));
    }

    for(int s = 0; s < nodes.getNumSegs(); s++){
        Segment seg = nodes.getSegByIx(s);
        seg.setTOF(seg.getTOF()*sys->getCharT());

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
 *  @brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  @param state_rot a 6-element non-dimensional state in rotating coordinates
 *  @param sys a pointer to the CR3BP system the state exists in
 *  @param t the non-dimensional time associated with the input state
 *  @param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP trajectory is input, selecting 
 *  <tt>centerIx = 0</tt> will produce Earth-centered inertial coordinates and selecting
 *  <tt>centerIx = 1</tt> will produce Moon-centered inertial coordinates.
 * 
 *  @return A state centered around the specified index in inertial, dimensional coordinates
 *  @throw Exception if <tt>centerIx</tt> is out of bounds
 */
std::vector<double> cr3bp_rot2inert_state(std::vector<double> state_rot, const SysData_cr3bp *sys, 
    double t, int centerIx){

    if(centerIx < 0 || centerIx > 1){
        throw Exception("tpat_calculations::cr3bp_rot2inert: Invalid center index");
    }

    Matrix3Rd I = Matrix3Rd::Identity(3,3);

    // Shift from Barycenter-centered coordinates to primary-centered coordinates
    Eigen::RowVector3d posShift = centerIx == 0 ? I.row(0)*(0-sys->getMu()) : I.row(0)*(1 - sys->getMu());
    
    // Angular velocity of rotating frame relative to inertial, non-dim units
    Eigen::RowVector3d omegaECI = I.row(2);

    double rotToInert[] = {cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1};
    Matrix3Rd DCM_R2I = Eigen::Map<Matrix3Rd>(rotToInert);

    Eigen::RowVector3d posRot = Eigen::Map<Eigen::Vector3d>(&(state_rot[0]));
    Eigen::RowVector3d velRot = Eigen::Map<Eigen::Vector3d>(&(state_rot[0])+3);

    // Shift, rotate, scale coordinates into inertial, dimensional
    Eigen::RowVector3d posInert = (posRot - posShift)*DCM_R2I*sys->getCharL();

    // Transform velocity into inertial, dimensional
    Eigen::RowVector3d velInert = (velRot + omegaECI.cross(posRot - posShift))*DCM_R2I*sys->getCharL()/sys->getCharT();

    // Concatenate data into a state vector
    std::vector<double> stateInert;
    stateInert.insert(stateInert.begin(), posInert.data(), posInert.data()+3);
    stateInert.insert(stateInert.end(), velInert.data(), velInert.data()+3);

    return stateInert;
}//==================================================================

//-----------------------------------------------------
//      BCR4BP Utility Functions
//-----------------------------------------------------

/**
 *  @brief Convert a CR3BP Sun-Earth trajectory into a BCR4BPR Sun-Earth-Moon trajectory
 *
 *  @param crTraj a CR3BP Sun-Earth trajectory
 *  @param bcSys a BCR4BPR Sun-Earth-Moon system data object; contains information about system
 *  scaling and orientation at time t = 0
 *  @param nodeID ID of a node in the set for which the epoch is known
 *  @param t0 the epoch at the specified node in BCR4BPR non-dimensional time units
 *
 *  @return a BCR4BPR Trajectory object
 *  @throw Exception if <tt>crTraj</tt> is not a Sun-Earth trajectory or if <tt>bcSys</tt> is
 *  not the Sun-Earth-Moon system.
 */
Traj_bc4bp bcr4bpr_SE2SEM(Traj_cr3bp crTraj, const SysData_bc4bp *bcSys, int nodeID, double t0){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_bc4bp bcNodeset = bcr4bpr_SE2SEM(static_cast<Nodeset_cr3bp>(crTraj), bcSys, nodeID, t0);
    return static_cast<Traj_bc4bp>(bcNodeset);
}//==================================================

/**
 *  @brief Convert a Sun-Earth CR3BP Nodeset to a Sun-Earth-Moon BCR4BPR Nodeset.
 *  
 *  @param crNodes a CR3BP Sun-Earth nodeset
 *  @param bcSys a BCR4BPR Sun-Earth-Moon system data object; contains information about system
 *  scaling and orientation at time t = 0
 *  @param nodeID ID of a node in the set for which the epoch is known
 *  @param t0 the epoch at the specified node in BCR4BPR non-dimensional time units
 *
 *  @return a BCR4BPR nodeset
 *  @throw Exception if <tt>crNodes</tt> is not a Sun-Earth trajectory or if <tt>bcSys</tt> is
 *  not the Sun-Earth-Moon system.
 */
Nodeset_bc4bp bcr4bpr_SE2SEM(Nodeset_cr3bp crNodes, const SysData_bc4bp *bcSys, int nodeID, double t0){
    if(crNodes.getSysData()->getPrimID(0) != 10 || crNodes.getSysData()->getPrimID(1) != 399){
        throw Exception("CR3BP trajectory is not in the Sun-Earth System");
    }

    if(bcSys->getPrimID(0) != 10 || bcSys->getPrimID(1) != 399 || bcSys->getPrimID(2) != 301){
        throw Exception("BCR4BPR system is not Sun-Earth-Moon");
    }

    // Create a BCR4BPR Nodeset
    Nodeset_bc4bp bcNodes(bcSys);

    double charL2 = crNodes.getSysData()->getCharL();
    double charT2 = crNodes.getSysData()->getCharT();
    double charL3 = bcSys->getCharL();
    double charT3 = bcSys->getCharT();

    crNodes.updateEpochs(nodeID, t0*charT3/charT2);

    // A mapping vector: index is the old node ID, value is the new node ID
    // All new IDs are initialized to the default INVALID_ID value
    std::vector<int> map_oldID_to_newID(crNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> bcNodeState, crNodeState;
    double epoch;
    for(int n = 0; n < crNodes.getNumNodes(); n++){
        bcNodeState.clear();
        crNodeState = crNodes.getStateByIx(n);
        for(int r = 0; r < ((int)crNodeState.size()); r++){
            if(r == 0)  // Convert x-coordinate, shift base to P2/P3 Barycenter
                bcNodeState.push_back(crNodeState[r]*charL2/charL3 - (1.0/bcSys->getK() - bcSys->getMu()));
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
 *  @brief Transform a BCR4BPR Sun-Earth-Moon nodeset into a CR3BP Sun-Earth nodeset
 *
 *  @param bcNodes a BCR4BPR Sun-Earth-Moon nodeset
 *  @param crSys a CR3BP Sun-Earth system data object; contains information about system
 *  scaling
 *
 *  @return a CR3BP nodeset object
 *  @throw Exception if <tt>crSys</tt> is not a Sun-Earth system or if <tt>bcNodes</tt> is
 *  not in the Sun-Earth-Moon system.
 */
Nodeset_cr3bp bcr4bpr_SEM2SE(Nodeset_bc4bp bcNodes, const SysData_cr3bp *crSys){
    if(crSys->getPrimID(0) != 10 || crSys->getPrimID(1) != 399){
        throw Exception("CR3BP system is not the Sun-Earth System");
    }

    if(bcNodes.getSysData()->getPrimID(0) != 10 || bcNodes.getSysData()->getPrimID(1) != 399 ||
            bcNodes.getSysData()->getPrimID(2) != 301){
        throw Exception("BCR4BPR trajectory is not in Sun-Earth-Moon");
    }

    const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(bcNodes.getSysData());

    // Create a BCR4BPR Trajectory
    Nodeset_cr3bp crNodes(crSys);

    double charL2 = crSys->getCharL();
    double charT2 = crSys->getCharT();
    double charL3 = bcNodes.getSysData()->getCharL();
    double charT3 = bcNodes.getSysData()->getCharT();

    std::vector<int> map_oldID_to_newID(bcNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> bcState, crNodeState;
    for(int n = 0; n < bcNodes.getNumNodes(); n++){
        bcState = bcNodes.getStateByIx(n);
        crNodeState.clear();

        for(int r = 0; r < ((int)bcState.size()); r++){
            if(r == 0)  // Shift origin to P2-P3 barycenter
                crNodeState.push_back((bcState[r] + 1.0/bcSys->getK() - bcSys->getMu())*charL3/charL2);
            else if(r < 3)   // Convert position
                crNodeState.push_back(bcState[r]*charL3/charL2);
            else if(r < 6)  // Convert velocity
                crNodeState.push_back(bcState[r]*(charL3/charL2)*(charT2/charT3));
        }

        Node node(crNodeState, bcNodes.getEpochByIx(n)*charT3/charT2);
        node.setExtraParam(0, DynamicsModel_cr3bp::getJacobi(&(crNodeState[0]), crSys->getMu()));
        
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
 *  @brief Transform a BCR4BPR Sun-Earth-Moon trajectory into a CR3BP Sun-Earth nodeset
 *
 *  @param bcTraj a BCR4BPR Sun-Earth-Moon trajectory
 *  @param crSys a CR3BP Sun-Earth system data object; contains information about system
 *  scaling
 *
 *  @return a CR3BP Trajectory object
 *  @throw Exception if <tt>crSys</tt> is not a Sun-Earth system or if <tt>bcTraj</tt> is
 *  not in the Sun-Earth-Moon system.
 */
Traj_cr3bp bcr4bpr_SEM2SE(Traj_bc4bp bcTraj, const SysData_cr3bp *crSys){
    // Process is identical for nodesets and trajectories, so cast to nodeset, perform transformation, and cast back
    Nodeset_cr3bp seNodeset = bcr4bpr_SEM2SE(static_cast<Nodeset_bc4bp>(bcTraj), crSys);
    return static_cast<Traj_cr3bp>(seNodeset);
}//==================================================

/**
 *  @brief Compute the location of the saddle point for a specific bicircular system
 *  and epoch
 *  @details This function uses a Newton-Raphson procedure to locate the zeros of the local
 *  acceleration field.
 *  @throws DivergeException if the Newton-Raphson procedure cannot converge
 *  @param bcSys system data object describing the bicircular system
 *  @param t0 the epoch at which the saddle point's location is computed
 * 
 *  @return the location of the saddle pointin BCR4BP rotating coordinates
 *  @throw DivergeException if the multiple shooting algorithm cannot converge on the
 *  saddle point location
 */
Eigen::Vector3d bcr4bpr_getSPLoc(const SysData_bc4bp *bcSys, double t0){

    // Compute approximate SP location using 3-body dynamics
    SysData_cr3bp subSys(bcSys->getPrimary(0), bcSys->getPrimary(1));
    double a = 2*subSys.getMu() - 1;
    double b = 4*subSys.getMu()*subSys.getMu() - 4*subSys.getMu() + 2;
    double c = 2*pow(subSys.getMu(), 3) - 3*subSys.getMu()*subSys.getMu() + 3*subSys.getMu() - 1;
    
    // x-coordinate in bcr4bpr
    Eigen::Vector3d spPos((-b + sqrt(b*b - 4*a*c))/(2*a*bcSys->getK()) - (1/bcSys->getK() - bcSys->getMu()), 0, 0);
    
    // Get primary positions
    double primPos[9];
    DynamicsModel_bc4bp::getPrimaryPos(t0, bcSys, primPos);
    
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
        double k = bcSys->getK();
        double mu = bcSys->getMu();
        double nu = bcSys->getNu();

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
        throw DivergeException("tpat_calculations::bcr4bpr_getSPLoc: Could not converge on SP location");
}//=====================================================

/**
 *  @brief Compute coefficients for 2nd-order polynomials in Epoch time that
 *  describe the x, y, and z coordinates of the saddle point.
 *  @details This function employs least squares to compute the coefficients. The
 *  number of points used and time span searched are hard coded in the function.
 * 
 *  @param bcSys data about the bicircular system
 *  @param T0 the "center" epoch; points are generated within +/- timeSpan of this
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
MatrixXRd bcr4bpr_spLoc_polyFit(const SysData_bc4bp *bcSys, double T0){
    const int STATE_SIZE = 3;   // predicting three position states
    double tSpan = 24*3600/bcSys->getCharT();   // Amount of time to span before and after T0
    int numPts = 100;
    double dT = 2*tSpan/((double)(numPts-1));

    MatrixXRd indVarMat(numPts, 3); // i.e. Vandermonde matrix
    MatrixXRd depVarMat(numPts, STATE_SIZE);
    for(int row = 0; row < 100; row++){
        double T = -tSpan + row*dT;     // Let T0 = 0 to center data and hopefully avoid numerical issues
        Eigen::Vector3d spPos = bcr4bpr_getSPLoc(bcSys, T + T0);
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

}// END of Astrohelion namespace


