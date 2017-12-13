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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstring>
#include <numeric>
#include <string>
#include <typeinfo>

#include "AsciiOutput.hpp"
#include "Arcset.hpp"
#include "Arcset_2bp.hpp"
#include "Arcset_bc4bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_periodic.hpp"
#include "BaseArcset.hpp"
#include "BodyData.hpp"
#include "Common.hpp"
#include "MultShootEngine.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "SysData_2bp.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Utilities.hpp"

#include <cspice/SpiceUsr.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

namespace astrohelion{

//-----------------------------------------------------
//      General Utility Functions
//-----------------------------------------------------
/**
 * @addtogroup util
 * \{
 */


/**
 *  @brief convert a date to ephemeris time
 *  
 *  @param pDate a string representing the date. The string can be formatted in one
 *  of two ways. First, a Gregorian-style date: 'YYYY/MM/DD HH:II:SS' (UTC, 24-hour clock);
 *  The time 'HH:II:SS' can be ommited, and the function will assume the time is 0:0:00
 *  Second, a Julian date (UTC) can be input with the format 'jd #' where '#' represents
 *  the Julian date.
 *
 *  @return the J2000 ephemeris time, or number of seconds after Jan 1, 2000 at 0:0:00 UTC.
 *  @throws Exception if the SPICE kernels cannot be loaded: the kernel names and
 *  filepaths are located in the settings XML file
 */
double dateToEphemerisTime(const char *pDate){
    // Convert the date to ephemeris time
    double et = 0;
    str2et_c(pDate, &et);
    checkAndReThrowSpiceErr("dateToEphemerisTime error");

    return et;
}//==========================================

/**
 *  @brief Convert a gregorian date to Julian Date
 * 
 *  @param yr Year
 *  @param mo Month (Jan = 1, ..., Dec = 12)
 *  @param d Day of mongth
 *  @param h Hour (24-hr clock)
 *  @param m minute
 *  @param s second
 *  @return Julian Date (days)
 */
double gregorianToJulian(double yr, double mo, double d, double h, double m, double s){
    return 367.0*yr - floor(7.0*(yr + floor((mo + 9.0)/12.0))/4.0) +
        floor(275.0*mo/9.0) + d + 1721013.5 + ((5.0/60.0 + m)/60.0 + h + s/3600.0)/24.0;
}//==========================================

/**
 *  @brief Determine the Greenwich Sidereal Time (i.e., angle) at the specified date
 *  @details Input date must be in UT1 time (always within +/- 0.9 seconds of UTC thanks
 *  to leapseconds)
 *  
 *  @param yr Year
 *  @param mo Month (Jan = 1, ..., Dec = 12)
 *  @param d Day of mongth
 *  @param h Hour (24-hr clock)
 *  @param m minute
 *  @param s second
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
            astrohelion::printErr("Calculations::getMirrorMat: Mirror type is not implemented; returning identiy\n");
            return Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Identity();
    }
}//====================================================

/**
 *  @brief Compute the eigendata (values and, optionally, vectors) of the matrix A
 *  @details This method leverages matrix balancing to reduce round-off errors in
 *  the eigendata computation
 * 
 *  @param A a real matrix (passed by value because it will be manipulated in-place)
 *  @param pVecs pointer to a matrix of complex doubles; eigenvectors are stored here as columns,
 *  each column is associated with the eigenvalue of the same index
 * 
 *  @return A vector of complex doubles that contains the eigenvalues
 */
std::vector<cdouble> getBalancedEigData(MatrixXRd A, MatrixXRcd *pVecs){
    bool compVecs = !(pVecs == nullptr);

    // Balance the matrix
    unsigned int low = 0, hi = 0;
    std::vector<double> perms {};
    balanceMat(A, low, hi, perms);

    // Compute the eigenvalues of the balanced matrix (and maybe eigenvectors)
    Eigen::EigenSolver<MatrixXRd> eigensolver(A, compVecs);
    if(eigensolver.info() != Eigen::Success)
        throw Exception("Calculations:getEigData: Could not compute eigenvalues of A matrix");

    // Extract the eigen data to return to the user
    std::vector<cdouble> eigVals {};
    Eigen::VectorXcd vals = eigensolver.eigenvalues();
    eigVals.insert(eigVals.begin(), vals.data(), vals.data() + vals.size());

    if(pVecs != nullptr){
        // Compute the eigenvectors of the balanced matrix
        *pVecs = eigensolver.eigenvectors();

        // backtransform the eigenvectors
        eigVec_backTrans(low, hi, perms, *pVecs);

    }

    return eigVals;
}//====================================================

/**
 *  @brief Balance a matrix, `A`, to decrease rounding errors when computing eigenvalues
 *  
 *  @details The eigenvalues of the balanced matrix are equivalent to the
 *  eigenvalues of the un-balanced matrix if no rounding errors are present. However,
 *  the eigenvectors of the balanced matrix may be very different. Use the `eigVec_backTrans`
 *  function to retrieve the eigenvectors of the original matrix from the eigenvectors of the
 *  balanced matrix.
 *  
 *  The balancing algorithm may time-out in some cases, thus, a maximum of 500 iterations are
 *  allowed to balance the STM. If this limit is reached, the function returns as usual but the
 *  balanced A matrix may not be balanced as well as it could be. (Reason for infinite loop is
 *  unknown; could be a coding error, though unlikely)
 * 
 *  @param A The matrix to balance. This matrix is manipulated in place, i.e., the original
 *  values are not preserved.
 *  @param low An integer that describes the form of the balanced (and permuted) matrix.
 *  The value passed in will be overwritten and filled with the value for the balanced matrix.
 *  @param hi An integer that describes the form of the balanced (and permuted) matrix.
 *  The value passed in will be overwritten and filled with the value for the balanced matrix.
 *  The values `hi` and `low` are defined such that A(i,j) = 0 if i > j AND
 *  j = 1, ... , (low-1) or i = (hi+1), ... , n
 *  
 *  @param perms A vector that stores the scaling and permutation information for the balanced
 *  matrix. The vector passed in will be overwritten with the data for the balanced matrix.
 */
void balanceMat(MatrixXRd& A, unsigned int& low, unsigned int& hi, std::vector<double>& perms){
    unsigned int n = A.rows();
    low = 0;
    hi = n-1;

    double radix = static_cast<double>(RADIX);
    double radix2 = radix*radix;        // Radix squared
    perms.clear();
    perms.assign(n, 0);

    // -------------------------------------------
    //    Step 1: Preliminary Permutations
    // -------------------------------------------
    /*
     *   After the operations below, the matrix is in this form:
     *       ________________________
     *   1   | x x x :       :       |
     *       | 0 x x :   X   :   Y   |
     *       | 0 0 x :       :       |
     *       |.......................|
     *   low |       :       :       |    
     *       |   0   :   Z   :   W   |
     *   hi  |       :       :       |
     *       |.......................|
     *       |       :       : x x x |
     *       |   0   :   0   : 0 x x |
     *   n   |       :       : 0 0 x |
     *       -------------------------
     *
     *   "x" indicates elements which are not necessarily zero
     *   low =/= hi
     */

    // Search for rows isolating an eigenvalue and push them down
    bool doLoop = true;
    while(doLoop){
        doLoop = false;
        for(int r = hi; r >= 0; r--){
            // Sum non-diagonal elements
            double rowSum = 0;
            for(int c = 0; c <= static_cast<int>(hi); c++){   // Only look at columns between 1 and hi
                if(c != r)
                    rowSum += std::abs(A(r,c));
            }

            // If all the non-diagonal elements sum to zero, the element on the
            // diagonal is an eigenvalue
            if(std::abs(rowSum) < 1e-12){
                exchange(r, hi, perms, hi, low, A);
                if(hi > 0){
                    doLoop = true;
                    hi--;
                }
                break;
            }
        }
    }

    // Search for columns isolating an eigenvalue and push them left
    doLoop = true;
    while(doLoop){
        doLoop = false;
        for(unsigned int c = low; c <= hi; c++){
            // Sum non-diagonal elements
            double colSum = 0;
            for(unsigned int r = low; r <= hi; r++){
                if(r != c)
                    colSum += std::abs(A(r,c));
            }

            // If all the non-diagonal elements sum to zero, the element on the
            // diagonal is an eigenvalue
            if(std::abs(colSum) < 1e-12){
                exchange(c, low, perms, hi, low, A);
                low++;
                // doLoop = true;
                doLoop = low < A.rows();
                break;
            }
        }
    }

    // These elements have not been permuted
    for(unsigned int i = low; i <= hi; i++){
        perms[i] = 1;
    }

    bool noConv = false;
    unsigned int count = 0, maxCount = 500;
    while(!noConv && count < maxCount){
        for(unsigned int i = low; i <= hi; i++){
            double c = 0, r = 0;
            for(unsigned int j = low; j <= hi; j++){
                if(i != j){
                    c += std::abs(A(j,i));
                    r += std::abs(A(i,j));
                }
            }

            double g = r/radix;
            double f = 1.0, s = c + r;

            while(c < g){
                f *= radix;
                c *= radix2;
            }

            g = r*radix;

            while( c >= g){
                f /= radix;
                c /= radix2;
            }

            if( (c+r)/f < 0.95*s ){
                g = 1.0/f;
                perms[i] *= f;
                noConv = true;

                for(unsigned int j = low; j < n; j++){
                    A(i,j) *= g;
                }
                for(unsigned int j = 0; j <= hi; j++){
                    A(j,i) *= f;
                }
            }
        }
        count++;
    }
    // printf("  Count = %u\n", count);
}//====================================================

/**
 *  @brief Exchange rows and columns to reach the appropriate matrix structure
 *  @details This sub-routine is leveraged in the `balanceMat()` algorithm
 *  to arrange the matrix, A, which is being balanced, in a way that isolates trivial
 *  rows and columns to simplify the computations.
 *  
 *  Source: Wilkinson, Reinsch: Handbook for Automatic Computation, pp. 315-326
 * 
 *  @param ix index of a row/column of A
 *  @param ix_move index of a row/column of A 
 *  @param perms vector containing permutation and scaling information
 *  @param hi the current value of `hi` from the `balanceMat()` algorithm
 *  @param low the current value of `low` from the `balanceMat()` algorithm
 *  @param A The matrix that is being balanced
 */
void exchange(unsigned int ix, unsigned int ix_move, std::vector<double>& perms, unsigned int hi,
    unsigned int low, MatrixXRd& A){

    unsigned int n = A.rows();

    char msg[128];
    if(ix >= A.cols() || ix >= n){
        sprintf(msg, "Calculations:exchange: ix = %u is out of range for A (%ux%ld)", ix, n, A.cols());
        throw Exception(msg);
    }

    if(ix_move >= A.cols() || ix_move >= n){
        sprintf(msg, "Calculations:exchange: ix_move = %u is out of range for A (%ux%ld)", ix_move, n, A.cols());
        throw Exception (msg);
    }

    // Remember this exchange
    perms[ix_move] = ix;
    double temp = 0;

    if(ix != ix_move){
        // Exchange entries between columns ix and ix_move over rows 1:hi
        for(unsigned int i = 0; i <= hi; i++){
            temp = A(i, ix);
            A(i, ix) = A(i, ix_move);
            A(i, ix_move) = temp;
        }

        // Exchange entries between rows ix and ix_move over columns low:n
        for(unsigned int i = low; i < n; i++){
            temp = A(ix, i);
            A(ix, i) = A(ix_move, i);
            A(ix_move, i) = temp;
        }
    }
}//====================================================

/**
 *  @brief Backtransform right-hand eigenvectors of a balanced matrix to
 *  retrieve the original eigenvectors
 *  @details When a matrix is balanced via `balanceMat()`, the eigenvalues are
 *  preserved but the eigenvectors are not due to permutations and scaling. This method
 *  transforms the right-hand eigenvectors of the balanced matrix, `vecs`, to 
 *  the eigenvectors associated with the original, un-balanced matrix
 * 
 *  @param low value output from `balanceMat()`
 *  @param hi value output from `balanceMat()`
 *  @param perms vector containing permutation and scaling information; this data is used to 
 *  transform the eigenvectors back to the original set.
 *  @param vecs matrix of eigenvectors (stored as columns) associated with a balanced matrix.
 *  These vectors are manipulated in place during the backtransformation
 */
void eigVec_backTrans(unsigned int low, unsigned int hi, const std::vector<double>& perms, MatrixXRcd& vecs){

    unsigned int n = vecs.rows(), m = vecs.cols();

    // For debugging:
    // printf("eigVec_backTrans:\n  low = %u\n  hi = %u\n  perms = [", low, hi);
    // for(unsigned int i = 0; i < perms.size(); i++){
    //     printf("%.4f%s", perms[i], i < perms.size() - 1 ? ", " : "]\n");
    // }
    // printf("  vecs: %ux%u matrix\n", n, m);

    double s = 0;
    for(unsigned int r = low; r <= hi; r++){
        s = perms[r];   // for left-hand eigenvectors, scaling factor is s = 1/perms[r]

        // Scale these elements by the value in perms for the row
        for(unsigned int c = 0; c < m; c++){
            vecs(r, c) *= s;
        }
    }

    // Swap elements in the opposite order as they were originally permuted
    int k = 0;
    for(int i = low-1; i >= 0; i--){
        k = static_cast<int>(perms[i]);

        if(k != i){
            for(unsigned int j = 0; j < m; j++){
                cdouble temp = vecs(i,j);  // store value temporarily
                vecs(i,j) = vecs(k,j);
                vecs(k,j) = temp;
            }
        }
    }
    for(int i = hi+1; i < static_cast<int>(n); i++){
        k = static_cast<int>(perms[i]);

        if(k != i){
            for(unsigned int j = 0; j < m; j++){
                cdouble temp = vecs(i,j);  // store value temporarily
                vecs(i,j) = vecs(k,j);
                vecs(k,j) = temp;
            }
        }
    }
}//====================================================

/**
 *  @brief Sort eigenvalues
 *  @details This algorithm leverages three simple axioms to determine if the eigenvalues
 *  are in a consistent order.
 *  
 *      1. Pairs of eigenvalues should occur in sequence ( [0,1] or [2,3], etc.) and have
 *      a product equal to one (reciprocal pairs)
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
 *  @param eigVals Vector of all eigenvalues along a family of trajectories. It is assumed that there
 *  are six eigenvalues per trajectory.
 *  
 *  @param eigVecs Vector of 6x6 eigenvector matrices; the eigenvectors are columns of the matrix and their
 *  order is consistent with the order of the six eigenvalues associated with the same trajectory
 * 
 *  @return Rearranged indices that describe the correct order of the eigenvalues/vectors
 */
std::vector<unsigned int> sortEig(std::vector<cdouble> eigVals, std::vector<MatrixXRcd> eigVecs){
    std::vector<unsigned int> sortedIxs(eigVals.size(), 0);

    if(eigVals.size() == 0 || eigVecs.size() == 0){
        astrohelion::printErr("Calculations::sortEig: No eigenvalues and/or eigenvectors - easy to sort! :P\n");
        return sortedIxs;
    }

    const unsigned int nE = eigVecs[0].rows();

    if(eigVals.size() % nE != 0){
        printErr("Calculations::sortEig: Must have %un eigenvalues!\n", nE);
        return sortedIxs;
    }

    // Generate all permutations of the indices 0 through 5
    std::vector<unsigned int> vals(nE);
    std::iota(std::begin(vals), std::end(vals), 0); // Fill with 0, 1, ..., nE
    std::vector<unsigned int> ixPerms = generatePerms<unsigned int>(vals);
    std::vector<float> cost(ixPerms.size()/nE, 0);

    // Sort the first set of eigenvalues so that reciprocal pairs occur near one another
    unsigned int s = 0;
    for(unsigned int p = 0; p < ixPerms.size()/nE; p++){
        for(unsigned int i = 0; i < nE/2; i++){
            // Penalize configurations with pairs that are not reciprocal
            cost[p] += std::abs(1.0 - eigVals[s*nE + ixPerms[p*nE + 2*i]] * eigVals[s*nE + ixPerms[p*nE + 2*i + 1]]);
        }
    }

    // Find the minimum cost
    std::vector<float>::iterator minCostIt = std::min_element(cost.begin(), cost.end());
    unsigned int minCostIx = minCostIt - cost.begin();

    // Sort the eigenvectors and eigenvalues according to the minimum cost permutation
    // printf("Minimum cost for member %u is %f on permutation %u\n", s, *minCostIt, minCostIx);
    std::vector<cdouble> prevSortVal(nE,0);
    MatrixXRcd prevSortVec = MatrixXRcd::Zero(nE,nE);
    for(unsigned int i = 0; i < nE; i++){
        sortedIxs[s*nE + i] = ixPerms[minCostIx*nE + i];
        prevSortVal[i] = eigVals[s*nE + ixPerms[minCostIx*nE + i]];
        prevSortVec.col(i) = eigVecs[s].col(ixPerms[minCostIx*nE+i]).normalized();
    }

    // printf("[%f%+fi, %f%+fi, %f%+fi, %f%+fi, %f%+fi, %f%+fi]\n",
    //     std::real(eigVals[sortedIxs[s*nE+0]]), std::imag(eigVals[sortedIxs[s*nE+0]]),
    //     std::real(eigVals[sortedIxs[s*nE+1]]), std::imag(eigVals[sortedIxs[s*nE+1]]),
    //     std::real(eigVals[sortedIxs[s*nE+2]]), std::imag(eigVals[sortedIxs[s*nE+2]]),
    //     std::real(eigVals[sortedIxs[s*nE+3]]), std::imag(eigVals[sortedIxs[s*nE+3]]),
    //     std::real(eigVals[sortedIxs[s*nE+4]]), std::imag(eigVals[sortedIxs[s*nE+4]]),
    //     std::real(eigVals[sortedIxs[s*nE+5]]), std::imag(eigVals[sortedIxs[s*nE+5]]));

    // Analyze all other eigenvalue sets using axioms while maintaining order of the
    // first set of eigenvalues/vectors
    for(s = 1; s < eigVals.size()/nE; s++){
        MatrixXRd dp_err = MatrixXRd::Zero(nE,nE);
        for(unsigned int i = 0; i < nE; i++){
            for(unsigned int j = i; j < nE; j++){
                dp_err(i,j) = std::abs(1.0 - std::abs(prevSortVec.col(i).dot(eigVecs[s].col(j).normalized())));
                dp_err(j,i) = dp_err(i,j);  // symmetric
            }
        }

        cost.clear();
        cost.assign(ixPerms.size()/nE, 0);
        for(unsigned int p = 0; p < ixPerms.size()/nE; p++){
            // bool printOut = s < 5 && (ixPerms[nE*p+0] == 0 && ixPerms[nE*p+1] == 1 && ixPerms[nE*p+2] == 2 && ixPerms[nE*p+3] == 3 &&
            //     ixPerms[nE*p+4] == 4 && ixPerms[nE*p+5] == 5);

            // if(printOut)
            //     printf("Eigenvalue set %u:\n", s);

            for(unsigned int i = 0; i < nE; i++){
                // Assign cost based on the dot product between this new arrangement
                // of eigenvectors and the previous arrangement
                // Ranges from 0 to 1
                cost[p] += dp_err(i, ixPerms[nE*p + i]);
                // if(printOut)
                //     printf("  dp cost %u: %.4e\n", i, dp_err(i, ixPerms[nE*p + i]));

                // Assign cost based on the distance from the previous sorted eigenvalues
                // Scale to give it more weight than some of the others
                cost[p] += 50*std::abs( (eigVals[nE*s + ixPerms[nE*p + i]] - prevSortVal[i])/prevSortVal[i] );
                // if(printOut)
                //     printf("  Dist cost %u: %.4e\n", i, 50*std::abs((eigVals[nE*s + ixPerms[nE*p + i]] - prevSortVal[i])/prevSortVal[i]));

                // Assign cost based on the reciprocal nature of the eigenvalues
                // This test could be removed to sort eigenvalues that don't occur in pairs
                if(i%2 == 1){
                    cdouble v1 = eigVals[s*nE + ixPerms[p*nE + i-1]];
                    cdouble v2 = eigVals[s*nE + ixPerms[p*nE + i]];
                    double minVal = std::abs(v1) < std::abs(v2) ? std::abs(v1) : std::abs(v2);
                    cost[p] += std::abs(1.0 - std::abs(v1 * v2) ) / minVal;
                    // if(printOut)
                    //     printf("  Recip cost: %u: %.4e\n", i, sqrt(std::abs(1.0 - std::abs(eigVals[s*nE + ixPerms[p*nE + i-1]] * eigVals[s*nE + ixPerms[p*nE + i]]))));
                }
            }
        }

        // Find the minimum cost
        minCostIt = std::min_element(cost.begin(), cost.end());
        minCostIx = minCostIt - cost.begin();

        // printf("Minimum cost for member %u is %f on permutation %u\n", s, *minCostIt, minCostIx);
        for(unsigned int i = 0; i < nE; i++){
            sortedIxs[s*nE + i] = ixPerms[minCostIx*nE + i];
            prevSortVal[i] = eigVals[s*nE + sortedIxs[s*nE + i]];
            prevSortVec.col(i) = eigVecs[s].col(sortedIxs[s*nE + i]).normalized();
        }
    }

    return sortedIxs;
}//====================================================

/**
 *  @brief Compute the stability index of a periodic orbit from a set of eigenvalues
 *  @details This algorithm assumes the orbit is periodic and that the eigenvalues 
 *  have been sorted (so they come in pairs).
 * 
 *  @param eigs A 6-element vector of eigenvalues associated with a periodic orbit
 *  @return the stability index, or NAN if no real, reciprocal eigenvalue pair is found
 *  @throws Exception if `eigs` does not have six elements
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
 *  @brief Use numerical integration to find a node on a trajectory at
 *  the specified time
 *  @details 
 * 
 *  @param traj An arcset 
 *  @param t The time at which the node is located
 * 
 *  @return A node on the specified trajectory at time `t`
 */
Node interpPointAtTime(const Arcset *traj, double t){
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

    Arcset *temp = new Arcset(traj->getSysData());
    SimEngine sim;
    sim.setVerbosity(Verbosity_tp::NO_MSG);
    sim.setNumSteps(2);
    sim.setVarStepSize(false);
    sim.setRevTime(t - node.getEpoch() < 0);
    sim.runSim(node.getState(), node.getEpoch(), t - node.getEpoch(), temp);

    Node node2 = temp->getNodeByIx(-1);
    delete(temp);
    return node2;
}//====================================================

void reconstructArc(const Arcset *pArcIn, Arcset *pArcOut){
    if(pArcIn == nullptr || pArcOut == nullptr)
        throw Exception("Calculations:reconstructArc: pArcIn or pArcOut are nullptrs; cannot proceed");

    const DynamicsModel *pModel = pArcIn->getSysData()->getDynamicsModel();

    // Construct a MultShootData from the input arcset
    MultShootData it(pArcIn);
    it.pArcOut = pArcOut;
    it.allCons = pArcIn->getAllConstraints();   // Copy for output
    pModel->multShoot_initDesignVec(it);
    
    // Propagate all segments from the design vector
    SimEngine sim;
    sim.setVerbosity(Verbosity_tp::NO_MSG);
    sim.setMakeDefaultEvents(false);    // No events - just propagate
    MultShootEngine::propSegsFromFreeVars(it, sim);

    // Create a full arcset from the propagated segments
    pModel->multShoot_createOutput(it);
}//====================================================

//-----------------------------------------------------
//      Orbit Determination Utility Functions
//-----------------------------------------------------


/**
 *  @brief Determine a set of spherical coordinates the describe 
 *  the position vector specified by x, y, and z
 *  @details The longitude angle measures in the in-plane angle, measured 
 *  from the +x-axis in a right-handed rotationg about the +z-axis. Latitude
 *  represents the out-of-plane angle measured from the xy-plane, with positive
 *  angles representing points with positive z-coordinates.
 * 
 *  @param x x-coordinate
 *  @param y y-coordinate
 *  @param z z-coordinate
 *  
 *  @return a vector containing {lat, long, R} in radians
 *  and the distance units of x, y, and z. Latitude takes a 
 *  value between -pi/2 and pi/2 while longitude takes a value
 *  between -pi and pi
 */
std::vector<double> getSpherical(double x, double y, double z){
    double R = sqrt(x*x + y*y + z*z);   // Distance from center of body to point
    if(R == 0)
        return std::vector<double>(3,0);

    double unit[] = {x/R, y/R, z/R};      // Unit vector pointing from center to point

    double lon = atan2(unit[1], unit[0]);           // longitude
    double lat = atan2(unit[2], sqrt(unit[0]*unit[0] + unit[1]*unit[1]));    // latitude

    std::vector<double> values {lat, lon, R};
    return values;
}//====================================================

/**
 *  @brief Convert inertial coordinates to local tangent coordinates
 * 
 *  @param inertPos Object position in inertial, cartesian (x,y,z) coordinates
 *  @param lat Latitude angle, radians, measured from the inertial xy plane;
 *  above the plane (z > 0) is a positive latitude.
 *  @param lon Longitude angle, radians, measured from reference meridian in a positive,
 *  right-handed rotation about inertial z
 *  @param theta_mer Meridian longitude, radians, measured from inertial x in a positive,
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
 *  @brief Convert local tangent coordinates to inertial coordinates
 * 
 *  @param localPos Object position in local tangent, cartesian (East, North, Up) coordinates
 *  @param lat Latitude angle, radians, measured from the inertial xy plane;
 *  above the plane (z > 0) is a positive latitude.
 *  @param lon Longitude angle, radians, measured from reference meridian in a positive,
 *  right-handed rotation about inertial z
 *  @param theta_mer Meridian longitude, radians, measured from inertial x in a positive,
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
 *  @brief Get the local tangent coordinates of an object given its 
 *  azimuth, elevation, and range
 *  @details [long description]
 * 
 *  @param s range distance
 *  @param az azimuth, measured from North toward East, radians
 *  @param el elevation, measured from local horizontal, radians
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
 *  @brief Compute the Keplarian elements at all steps/nodes
 *  of a 2-body trajectory or arcset
 * 
 *  @param pSet Pointer to a arcset or trajectory
 */
void r2bp_computeAllKepler(BaseArcset *pSet){
    if(pSet->getSysData()->getType() == SysData_tp::R2BP_SYS){
        const SysData_2bp *pSys = static_cast<const SysData_2bp*>(pSet->getSysData());
        
        for(unsigned int n = 0; n < pSet->getNumNodes(); n++){
            Node& node = pSet->getNodeRefByIx(n);
            r2bp_computeKepler(pSys, &node);
        }
    }else{
        printErr("r2bp_computeAllKepler: Arcset object is not associated with the 2BP");
    }
}//====================================================

/**
 *  \ingroup 2bp
 *  @brief Compute the Keplarian orbital elements at a 
 *  specific node
 * 
 *  @param pSys Pointer to the dynamical system the node exists in
 *  @param pNode Pointer to the node
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
    pNode->setExtraParam(PARAMKEY_SMA, a);
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
 *  @brief [brief description]
 *  @details [long description]
 * 
 *  @param pSys [description]
 *  @param a [description]
 *  @param e [description]
 *  @param argPeri [description]
 *  @param i [description]
 *  @param RAAN [description]
 *  @param TA [description]
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
 *  @brief Compute the magnitude of a velocity component given Jacobi Constant
 * 
 *  @param s state vector (non-dimensional); MUST contain at least the 6 position and velocity states.
 *  Note that the desired velocity component identified by velIxToFind is not required; put a placeholder
 *  zero or NAN (or anything, really) in its place; this value will not be used in computations.
 *  @param mu non-dimensional system mass ratio
 *  @param C Jacobi constant value
 *  @param velIxToFind index of the velocity component to compute (i.e. 3, 4, or 5)
 *  @return the magnitude of vx (3), vy (4), or vz (5).
 *  @throws Exception if `velIxToFind` is out of bounds
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
 *  @brief Generate an arcset that is a suitable guess for a symmetric periodic orbit
 *  in the CR3BP
 *  
 *  @details The output half-period arcset can be passed to cr3bp_correctHalfPerSymPO for corrections.
 *  
 *  Note: This method ignores all crash events, so it is possible to compute a
 *  periodic orbit that passes through a primary. 
 *  
 *  @param pSys the dynamical system
 *  @param ic non-dimensional initial state vector
 *  @param halfPeriod non-dimensional orbit period
 *  @param halfPerNumNodes the number of nodes to use for HALF of the periodic orbit; more nodes 
 *  may result in a more robust correction
 *  @param order the number of revolutions about the system/primary this orbit completes before
 *  it repeats periodically. Think of a period-3 DRO (order = 3) or a butterfly (order = 2)
 *  @param mirrorType how this periodic orbit mirrors in the CR3BP
 *  @param fixedStates a vector containing the indices of which initial states
 *  we would like to fix. Not all states are possible for each mirror condition.
 *  See the enum definition for specific details.
 *  @param tol tolerance to use in the corrections process
 *  
 *  @return An arcset that represents half of the desired periodic orbit and is constrained
 *  according to the specified mirror conditions
 *  
 *  @throws Exception if `mirrorType` is invalid
 *  @throws DivergeException if the multiple shooting algorithm cannot converge on a 
 *  mirrored solution.
 */
Arcset_cr3bp cr3bp_propHalfPerSymPO(const SysData_cr3bp *pSys, std::vector<double> ic, double halfPeriod, unsigned int halfPerNumNodes,
    unsigned int order, Mirror_tp mirrorType, std::vector<unsigned int> fixedStates, double tol){

    SimEngine sim;                              // Engine to perform simulation
    sim.setAbsTol(tol < 1e-12 ? 1e-15 : tol/1000.0);
    sim.setRelTol(sim.getAbsTol());
    sim.setMakeDefaultEvents(false);            // Ignore any crashes into the primaries
    
    // Determine which states must be zero for mirroring
    Event mirrorEvt;
    switch(mirrorType){
        case Mirror_tp::MIRROR_XZ:
            mirrorEvt.createEvent(Event_tp::XZ_PLANE, 0, true);     // Tell the sim to quit once it reaches the XZ plane
            break;
        case Mirror_tp::MIRROR_YZ:
            mirrorEvt.createEvent(Event_tp::YZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_XY:
            mirrorEvt.createEvent(Event_tp::YZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_X_AX_H:
            mirrorEvt.createEvent(Event_tp::XZ_PLANE, 0, true);
            break;
        case Mirror_tp::MIRROR_X_AX_V:
            mirrorEvt.createEvent(Event_tp::XY_PLANE, 0, true);
            break;
        default:
            throw Exception("Mirror type either not defined or not implemented");
    }

    mirrorEvt.setStopCount(order);
    sim.addEvent(mirrorEvt);

    // Run the sim until the event is triggered
    Arcset_cr3bp halfOrbArc(pSys);
    sim.runSim_manyNodes(ic, halfPeriod, halfPerNumNodes, &halfOrbArc);

    // Check to make sure the simulation ended with the event (not running out of time)
    if(halfOrbArc.getNodeByIx(-1).getTriggerEvent() != mirrorEvt.getType()){
        printErr("Calculations::cr3bp_getPeriodic: simulation of half-period orbit did not end in mirror event; may have diverged\n");
        // halfOrbArc.saveToMat("halfOrbArc_failedMirror.mat");
    }

    // Create constraints to ensure the mirroring is preserved and the fixedStates are fixed
    cr3bp_addMirrorCons(&halfOrbArc, mirrorType, fixedStates);

    // Arcset has been created and constrained
    return halfOrbArc;
}//====================================================

/**
 *  @brief Correct an arcset to be periodic using symmetry properties
 * 
 *  @param halfPerGuess an approximation of the first half of the periodic orbit. This arcset
 *  must include the constraints necessary to enforce the mirror conditions at either end of the half-period
 *  arc.
 *  
 *  @param halfPerCorrected the corrected half-period solution is stored here. Set to nullptr to
 *  discard this solution once the corrections process is complete.
 *  @param mirrorType Describes how the orbit is symmetric
 *  @param tol tolerance to enforce during the multiple-shooting process
 *  @param pItData iteration data from the multiple-shooting process is stored here. Set to nullptr
 *  to discard this data once the corrections process is complete.
 *  
 *  @return A full revolution of the periodic orbit. Symmetry properties are leveraged to construct
 *  the second half of the orbit from the first half, which saves computation time and improves
 *  numerical accuracy.
 */
Arcset_periodic cr3bp_correctHalfPerSymPO(const Arcset_cr3bp *halfPerGuess, Arcset_cr3bp *halfPerCorrected, Mirror_tp mirrorType,
    double tol, MultShootData *pItData){
    
    if(halfPerGuess == nullptr)
        throw Exception("Calculations::cr3bp_correctHalfPerSymPO: Cannot proceed with a nullptr input pointer for halfPerGuess");

    const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp *>(halfPerGuess->getSysData());

    bool createCorrectedArcset = halfPerCorrected == nullptr;
    if(createCorrectedArcset){
        halfPerCorrected = new Arcset_cr3bp(pSys);
    }

    if(halfPerCorrected->getSysData() != halfPerGuess->getSysData()){
        if(createCorrectedArcset){
            delete halfPerCorrected;
            halfPerCorrected = nullptr;
        }
        throw Exception("Calculations::cr3bp_correctHalfPerSymPO: input arcsets have different system data objects");
    }

    // Use differential corrections to enforce the mirror conditions
    MultShootEngine corrector;
    corrector.setTol(tol);
    corrector.setIgnoreCrash(true); // Corrector also ignores crash events
    corrector.setTOFType(MSTOF_tp::VAR_EQUALARC);
    // corrector.setVerbosity(Verbosity_tp::ALL_MSG);

    // If the user passed in a nullptr pointer, create a temporory data object
    // on the stack to avoid seg faults, then delete it before exiting to 
    // avoid memory leaks
    bool createdTempMSData = false;
    if(!pItData){
        createdTempMSData = true;
        pItData = new MultShootData(halfPerGuess);
    }


    // printf("\n************************************\nHalf Period Guess:\n************************************\n");
    // halfPerGuess->print();
    // halfPerGuess->saveToMat("data/temp_halfPerGuess.mat");

    *halfPerCorrected = Arcset_cr3bp(pSys);   // reset
    
    try{
        // Try once with the standard multiple-shooting parameters
        *pItData = corrector.multShoot(halfPerGuess, halfPerCorrected);
        // halfPerCorrected.print();
        // halfPerCorrected.saveToMat("data/temp_halfPerCorrected.mat");
    }catch(DivergeException &e){
        // Try again using a line search
        corrector.setDoLineSearch(true);
        corrector.setMaxIts(250);
        *halfPerCorrected = Arcset_cr3bp(pSys);   // reset
        
        // printf("\n************************************\nHalf Period Guess:\n************************************\n");
        // halfPerGuess->print();
        // halfPerGuess->saveToMat("data/temp_halfPerGuess.mat");

        try{
            corrector.multShoot(halfPerGuess, halfPerCorrected);
            
            // halfPerCorrected->print();
        }catch(DivergeException &ee){
            if(createCorrectedArcset){
                delete halfPerCorrected;
                halfPerCorrected = nullptr;
            }
            throw DivergeException("Calculations::cr3bp_getPeriodic: Could not converge half-period arc with mirroring condition");
        }
    }

    // Construct the full periodic orbit from the half-period arc
    Arcset_periodic po(pSys);
    cr3bp_halfPO2fullPO(halfPerCorrected, &po, mirrorType);
    
    // Free the iteration data if it was created on the stack
    if(createdTempMSData){
        delete pItData;
        pItData = nullptr;
    }

    // Free the arccset data if it was created on the stack
    if(createCorrectedArcset){
        delete halfPerCorrected;
        halfPerCorrected = nullptr;
    }

    return po;     // Now contains entire trajectory
}//====================================================

/**
 *  @brief Construct a full periodic orbit (PO) from a half-period representation of the orbit
 *  @details The half-period arc must start and end at mirror conditions. This function leverages
 *  the symmetry properties of the CR3BP to construct the second half of the full PO without 
 *  numerical integration.
 * 
 *  @param halfPerArc An arc that represents exactly half of the full periodic orbit
 *  @param pFullPO A pointer to an object that stores the full periodic orbit
 *  @param mirrorTp Describes how the orbit is symmetric
 *  
 *  @throws Exception if either arcset pointer is null
 */
void cr3bp_halfPO2fullPO(const Arcset_cr3bp *halfPerArc, Arcset_periodic *pFullPO, Mirror_tp mirrorTp){
    if(halfPerArc == nullptr)
        throw Exception("cr3bp_halfPO2fullPO: Cannot construct an orbit from a null half-period arc");

    if(pFullPO == nullptr)
        throw Exception("cr3bp_halfPO2fullPO: Cannot store full periodic orbit in null Arcset");

    // Make a copy for the final output periodic orbit
    *pFullPO = Arcset_periodic(*halfPerArc);

    // printf("\n************************************\nHalf Period Corrected:\n************************************\n");
    // halfPerArc->print();
    // halfPerArc->saveToMat("data/temp_halfPerCorrected.mat");
    // waitForUser();

    double halfTOF = pFullPO->getTimeByIx(-1);
    int halfPerTraj_len = static_cast<int>(pFullPO->getNumNodes());
    
    // Use Mirror theorem to create the second half of the orbit
    double I[] = {  1, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0,
                    0, 0, 1, 0, 0, 0,
                    0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 1};
    // Symplecticity matrix
    double S_data[] = {0, -2, 0, 1, 0, 0,
                        2, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 1,
                        -1, 0, 0, 0, 0, 0,
                        0, -1, 0, 0, 0, 0,
                        0, 0, -1, 0, 0, 0};
    MatrixXRd S = Eigen::Map<MatrixXRd>(S_data, 6, 6);

    // Inverse of Symplecticity matrix
    double SI_data[] = { 0, 0, 0, -1, 0, 0,
                        0, 0, 0, 0, -1, 0,
                        0, 0, 0, 0, 0, -1,
                        1, 0, 0, 0, -2, 0,
                        0, 1, 0, 2, 0, 0,
                        0, 0, 1, 0, 0, 0};
    MatrixXRd SI = Eigen::Map<MatrixXRd>(SI_data, 6, 6);
    MatrixXRd G = getMirrorMat(mirrorTp);

    MatrixXRd Phi_total_mirror = MatrixXRd::Identity(6,6);    // cumulative STM from the halfway point forward on the mirrored side
    MatrixXRd Phi_total_prop = MatrixXRd::Identity(6,6);      // cumulative STM to the halfway point on the propagated side

    int prevID = pFullPO->getNodeByIx(halfPerTraj_len-1).getID();
    for(int i = halfPerTraj_len-2; i >= 0; i--){
        // Use mirroring to populate second half of the orbit
        std::vector<double> state = pFullPO->getStateByIx(i);
        Eigen::RowVectorXd stateVec = Eigen::Map<Eigen::RowVectorXd>(&(state[0]), 1, 6);
        Eigen::RowVectorXd newStateVec = stateVec*G;

        // Create a node with the mirrored state and appropriate epoch
        Node node;
        node.setState(newStateVec.data(), newStateVec.cols());
        node.setEpoch(2*halfTOF - pFullPO->getTimeByIx(i));
        int id = pFullPO->addNode(node);

        Segment seg(prevID, id, pFullPO->getEpoch(id) - pFullPO->getEpoch(prevID));
        seg.setStateWidth(pFullPO->getSegRefByIx(0).getStateWidth());

        Segment &mirrorSeg = pFullPO->getSegRefByIx(i);   // Assumes segments all move linearly in time (simple arcset)
        
        // Get first and last state vectors on the segment (minimum required to save family members)
        std::vector<double> q0 = mirrorSeg.getStateByRow(0);
        std::vector<double> qf = mirrorSeg.getStateByRow(-1);

        // Map to Eigen vector objects and rotate via mirror matrix
        Eigen::RowVectorXd v0 = Eigen::Map<Eigen::RowVectorXd>(&(q0[0]), 1, 6);
        Eigen::RowVectorXd vf = Eigen::Map<Eigen::RowVectorXd>(&(qf[0]), 1, 6);
        Eigen::RowVectorXd v0_m = v0*G;
        Eigen::RowVectorXd vf_m = vf*G;

        // Update the cumulative STM on the propagated half
        Phi_total_prop = Phi_total_prop * pFullPO->getSTMByIx(i);
        // Compute the STM on the mirrored side
        MatrixXRd stm = SI * G * Phi_total_prop.transpose() * G * Phi_total_mirror.transpose() * S;
        // Update the cumulate STM on the mirrored half
        Phi_total_mirror = stm*Phi_total_mirror;

        unsigned int ctrl_dim = mirrorSeg.getCtrlLaw() ? mirrorSeg.getCtrlLaw()->getNumStates() : 0;
        unsigned int extra_dim = pFullPO->getSysData()->getDynamicsModel()->getExtraStateSize();

        // Create state "filler" to take up extra space (any states not included in the core 6)
        std::vector<double> filler;
        if(ctrl_dim > 0){
            std::vector<double> ctrl(q0.begin()+6, q0.begin()+6+ctrl_dim);
            filler.insert(filler.end(), ctrl.begin(), ctrl.end());
        }
        filler.insert(filler.end(), I, I+36);    // Add initial STM data
        if(extra_dim > 0){
            for(unsigned int k = 0; k < extra_dim; k++){
                filler.push_back(0);    // add a zero for the extra state
            }
        }

        // Save states to the segment
        seg.appendState(v0_m.data(), v0_m.cols());
        seg.appendState(filler);

        // Repeat for final state
        filler.clear();
        if(ctrl_dim > 0){
            std::vector<double> ctrl(qf.begin()+6, qf.begin()+6+ctrl_dim);
            filler.insert(filler.end(), ctrl.begin(), ctrl.end());
        }
        filler.insert(filler.end(), stm.data(), stm.data() + stm.rows()*stm.cols());    // Add STM data
        if(extra_dim > 0){
            for(unsigned int k = 0; k < extra_dim; k++){
                filler.push_back(0);    // add a zero for the extra state
            }
        }

        seg.appendState(vf_m.data(), vf_m.cols());
        seg.appendState(filler);
        seg.setSTM(stm);

        // Save times
        seg.appendTime(pFullPO->getNodeRef(prevID).getEpoch());
        seg.appendTime(node.getEpoch());
        seg.updateTOF();

        prevID = id;    // Update variable
        pFullPO->addSeg(seg); // Add the segment
    }
}//====================================================

/**
 *  \ingroup cr3bp
 *  @brief Add mirror and state constraints to the first and last nodes on an arcset
 *  @details A mirror constraint is constructed based on the mirrorType specified in the
 *  input arguments. Additionally, a state constraint is constructed and applied to only 
 *  the first node based on the fixedStates specified in the input arguments. Any pre-existing
 *  constraints on the initial and final nodes are deleted
 * 
 *  @param pArc pointer to the arcset
 *  @param mirrorType Describes how the trajectory is mirrored at both the initial and 
 *  final nodes.
 *  @param fixedStates The indices of any states at the initial node that should be fixed
 */
void cr3bp_addMirrorCons(Arcset_cr3bp *pArc, Mirror_tp mirrorType, std::vector<unsigned int> fixedStates){
    // Determine which states must be zero for mirroring
    std::vector<double> conData0(6, NAN);
    switch(mirrorType){
        case Mirror_tp::MIRROR_XZ:
            conData0[1] = 0; // y
            conData0[3] = 0; // x-dot
            conData0[5] = 0; // z-dot
            break;
        case Mirror_tp::MIRROR_YZ:
            conData0[0] = 0;    // x
            conData0[4] = 0;    // y-dot
            conData0[5] = 0;    // z-dot
            break;
        case Mirror_tp::MIRROR_XY:
            conData0[2] = 0;    // z
            conData0[3] = 0;    // x-dot
            conData0[4] = 0;    // y-dot
            break;
        case Mirror_tp::MIRROR_X_AX_H:
            conData0[1] = 0;    // y
            conData0[2] = 0;    // z
            conData0[3] = 0;    // x-dot
            break;
        case Mirror_tp::MIRROR_X_AX_V:
            conData0[1] = 0;    // y
            conData0[2] = 0;    // z
            conData0[3] = 0;    // x-dot
            break;
        default:
            throw Exception("Calculations::cr3bp_addMirrorCons: Mirror type not defined");
    }

    std::vector<double> conDataf = conData0;    // Same states must be zero at the halfway point

    // Constrain initial node as desired
    std::vector<double> ic = pArc->getStateByIx(0);
    for(unsigned int i = 0; i < fixedStates.size(); i++){
        bool okToFix = true;
        for(unsigned int n = 0; n < conData0.size(); n++){
            if(!std::isnan(conData0[fixedStates[i]])){
                printWarn("Cannot fix state %d; it must be zero for this mirror condition; ignoring\n", fixedStates[i]);
                okToFix = false;
                break;
            }
        }
        if(okToFix){
            conData0[fixedStates[i]] = ic[fixedStates[i]];
        }
    }

    Constraint con_0(Constraint_tp::STATE, pArc->getNodeRefByIx(0).getID(), conData0);
    Constraint con_f(Constraint_tp::STATE, pArc->getNodeRefByIx(-1).getID(), conDataf);
    pArc->getNodeRefByIx(0).clearConstraints();
    pArc->getNodeRefByIx(-1).clearConstraints();
    pArc->addConstraint(con_0);
    pArc->addConstraint(con_f);
}//====================================================

/**
 *  \ingroup cr3bp
 *  @brief Transition a arcset from the EM system to the SE system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  `thetaE0`, `thetaM0`, and `gamma`. Accordingly, the EM arcset
 *  should have node epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire set may be accomplished via the updateEpochs() function.
 *
 *  @param EMNodes a CR3BP Earth-Moon arcset
 *  @param pSESys a Sun-Earth CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Arcset_cr3bp cr3bp_EM2SE(Arcset_cr3bp EMNodes, const SysData_cr3bp *pSESys, double thetaE0, double thetaM0,
    double gamma){

    // astrohelion::printColor(BLUE, "Converting EM to SE\nEM Sys:\n  %d Nodes\n  %d Segments\n", EMNodes.getNumNodes(),
    //     EMNodes.getNumSegs());

    Arcset_cr3bp SENodes(pSESys);

    double charTE = EMNodes.getSysData()->getCharT();       // characteristic time in EM system
    double charLE = EMNodes.getSysData()->getCharL();       // characteristic length in EM system
    double charTS = pSESys->getCharT();                       // characteristic time in SE system
    double charLS = pSESys->getCharL();                       // characteristic length in SE system

    std::vector<int> map_oldID_to_newID(EMNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> state_SE;
    double epoch = 0;
    unsigned int core_dim = EMNodes.getSysData()->getDynamicsModel()->getCoreStateSize();
    unsigned int full_dim = (core_dim + 1)*core_dim + EMNodes.getSysData()->getDynamicsModel()->getExtraStateSize();

    for(unsigned int n = 0; n < EMNodes.getNumNodes(); n++){
        epoch = EMNodes.getEpochByIx(n);
        state_SE = cr3bp_EM2SE_state(EMNodes.getStateByIx(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, pSESys->getMu());
        
        map_oldID_to_newID[EMNodes.getNodeByIx(n).getID()] = SENodes.addNode(Node(state_SE, epoch*charTE/charTS));
    }

    // Convert time units for segment TOFs and update IDs of origin and terminus nodes
    for(unsigned int s = 0; s < EMNodes.getNumSegs(); s++){
        Segment seg = EMNodes.getSegByIx(s);
        seg.setSTM(MatrixXRd::Identity(core_dim, core_dim));

        // Get a copy of the state and time vectors and modify them to the proper coordinates / dimensions
        std::vector<double> segStates = seg.getStateVector();
        std::vector<double> segTimes = seg.getTimeVector();

        if(segTimes.size() != segStates.size()/full_dim)
            throw Exception("cr3bp_EM2SE::Segment time and state vectors are not the expected sizes.");

        for(unsigned int i = 0; i < segTimes.size(); i++){
            std::vector<double> newState = cr3bp_EM2SE_state(seg.getStateByRow(i), epoch, thetaE0, thetaM0,
                gamma, charLE, charTE, charLS, charTS, pSESys->getMu())
            ;
            for(unsigned int j = 0; j < newState.size(); j++){
                segStates[full_dim*i+j] = newState[j];
            }

            segTimes[i] *= charTE/charTS;
        }

        seg.setStateVector(segStates);
        seg.setTimeVector(segTimes);
        seg.updateTOF();

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
 *  @brief Transition a arcset from the SE system to the EM system
 *  
 *  The relative orientation between the two systems at time t = 0 is described by the three angles
 *  `thetaE0`, `thetaM0`, and `gamma`. Accordingly, the SE arcset
 *  should have node epochs such that t = 0 corresponds to the desired geometry; adjusting
 *  the epoch for the entire set may be accomplished via the updateEpochs() function.
 *
 *  @param SENodes a CR3BP Sun-Earth arcset
 *  @param pEMSys an Earth-Moon CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
Arcset_cr3bp cr3bp_SE2EM(Arcset_cr3bp SENodes, const SysData_cr3bp *pEMSys, double thetaE0, double thetaM0,
    double gamma){

    // astrohelion::printColor(BLUE, "Converting SE to EM\nSE Sys:\n  %d Nodes\n  %d Segments\n", SENodes.getNumNodes(),
    //     SENodes.getNumSegs());

    Arcset_cr3bp EMNodes(pEMSys);

    double charTE = pEMSys->getCharT();                   // characteristic time in EM system
    double charLE = pEMSys->getCharL();                   // characteristic length in EM system
    double charTS = SENodes.getSysData()->getCharT();    // characteristic time in SE system
    double charLS = SENodes.getSysData()->getCharL();    // characteristic length in SE system

    const SysData_cr3bp *pSESys = static_cast<const SysData_cr3bp*>(SENodes.getSysData());
    unsigned int core_dim = pSESys->getDynamicsModel()->getCoreStateSize();
    unsigned int full_dim = (core_dim + 1)*core_dim + pSESys->getDynamicsModel()->getExtraStateSize();

    std::vector<int> map_oldID_to_newID(SENodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> state_EM;
    double epoch = 0;

    for(unsigned int n = 0; n < SENodes.getNumNodes(); n++){
        epoch = SENodes.getEpochByIx(n);
        // Transform a single node
        state_EM = cr3bp_SE2EM_state(SENodes.getStateByIx(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, pSESys->getMu());
        
        map_oldID_to_newID[SENodes.getNodeByIx(n).getID()] = EMNodes.addNode(Node(state_EM, epoch*charTS/charTE));
    }

    // Convert time units for segment TOFs and update IDs of origin and terminus nodes
    for(unsigned int s = 0; s < SENodes.getNumSegs(); s++){
        Segment seg = SENodes.getSegByIx(s);
        seg.setSTM(MatrixXRd::Identity(core_dim, core_dim));

        // Get a copy of the state and time vectors and modify them to the proper coordinates / dimensions
        std::vector<double> segStates = seg.getStateVector();
        std::vector<double> segTimes = seg.getTimeVector();

        if(segTimes.size() != segStates.size()/full_dim)
            throw Exception("cr3bp_SE2EM::Segment time and state vectors are not the expected sizes.");

        for(unsigned int i = 0; i < segTimes.size(); i++){
            std::vector<double> newState = cr3bp_SE2EM_state(seg.getStateByRow(i), epoch, thetaE0, thetaM0,
                gamma, charLE, charTE, charLS, charTS, pSESys->getMu())
            ;
            for(unsigned int j = 0; j < newState.size(); j++){
                segStates[full_dim*i+j] = newState[j];
            }

            segTimes[i] *= charTS/charTE;
        }

        seg.setStateVector(segStates);
        seg.setTimeVector(segTimes);
        seg.updateTOF();

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
 *  \ingroup cr3bp
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
 *  \ingroup cr3bp
 *  @brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  @param arcset CR3BP arcset
 *  @param epoch0 Epoch associated with t = 0 in the CR3BP; epoch in seconds past J2000 (i.e., ephemeris time)
 *  @param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP trajectory is input, selecting 
 *  `centerIx = 1` will produce Earth-centered inertial coordinates and selecting
 *  `centerIx = 2` will produce Moon-centered inertial coordinates. A choice of
 *  `centerIx = 0` produces system barycenter-centered inertial coordinates.
 * 
 *  @return A arcset centered around the specified index in inertial, dimensional coordinates
 */
Arcset_cr3bp cr3bp_rot2inert(Arcset_cr3bp arcset, double epoch0, int centerIx){
    if(centerIx < 0 || centerIx > 2){
        throw Exception("Calculations::cr3bp_rot2inert: Invalid center index");
    }

    const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp *>(arcset.getSysData());
    unsigned int core_dim = pSys->getDynamicsModel()->getCoreStateSize();
    unsigned int full_dim = core_dim + core_dim*core_dim + pSys->getDynamicsModel()->getExtraStateSize();

    Arcset_cr3bp inertNodes(pSys);
    std::vector<int> map_oldID_to_newID(arcset.getNextNodeID(), Linkable::INVALID_ID);

    for(unsigned int i = 0; i < arcset.getNumNodes(); i++){
        std::vector<double> stateInert = cr3bp_rot2inert_state(arcset.getStateByIx(i), pSys, arcset.getEpochByIx(i), epoch0, centerIx);
        
        // Save the new ID and add the node, making sure to update Epoch as well!
        Node newNode(stateInert, arcset.getEpochByIx(i)*pSys->getCharT() + epoch0);
        map_oldID_to_newID[arcset.getNodeByIx(i).getID()] = inertNodes.addNode(newNode);
    }

    for(unsigned int s = 0; s < arcset.getNumSegs(); s++){
        Segment seg = arcset.getSegByIx(s);

        // Get a copy of the state and time vectors and modify them to the proper coordinates / dimensions
        std::vector<double> segStates = seg.getStateVector();
        std::vector<double> segTimes = seg.getTimeVector();

        if(segTimes.size() != segStates.size()/full_dim)
            throw Exception("cr3bp_rot2inert::Segment time and state vectors are not the expected sizes.");

        for(unsigned int i = 0; i < segTimes.size(); i++){
            std::vector<double> inertState = cr3bp_rot2inert_state(seg.getStateByRow(i), pSys, segTimes[i], epoch0, centerIx);
            for(unsigned int j = 0; j < inertState.size(); j++){
                segStates[full_dim*i+j] = inertState[j];
            }

            segTimes[i] *= pSys->getCharT();
        }

        seg.setStateVector(segStates);
        seg.setTimeVector(segTimes);
        seg.updateTOF();

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
 *  @brief Transform CR3BP rotating coordinates to inertial, dimensional coordinates
 * 
 *  @param stateRot a 6-element non-dimensional state in rotating coordinates
 *  @param pSys a pointer to the CR3BP system the state exists in
 *  @param t the non-dimensional time associated with the input state
 *  @param epoch0 Epoch associated with t = 0 in the CR3BP; epoch in seconds past J2000 (i.e., ephemeris time)
 *  @param centerIx The index of the primary that will be the center of the inertial
 *  frame. For example, if an Earth-Moon CR3BP trajectory is input, selecting 
 *  `centerIx = 1` will produce Earth-centered inertial coordinates and selecting
 *  `centerIx = 2` will produce Moon-centered inertial coordinates. A choice of
 *  `centerIx = 0` produces system barycenter-centered inertial coordinates.
 * 
 *  @return A state centered around the specified point in inertial, ecliptic J2000, dimensional coordinates
 *  @throw Exception if `centerIx` is out of bounds
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

    std::vector<double> stateInert(6,0);
    double r_p2[3], v_p2[3], angMom_p2[3], mag_angMom_p2, inst_charL, inst_charT;
    double unit_x[3], unit_y[3], unit_z[3];

    BodyData P1(pSys->getPrimID(0));
    BodyData P2(pSys->getPrimID(1));
    double sysGM = P1.getGravParam() + P2.getGravParam();

    // Locate P2 relative to P1 at the specified time
    spkezr_c(targ, epoch0 + t*pSys->getCharT(), ref, abcorr, obs, p2State, &lt);

    checkAndReThrowSpiceErr("cr3bp_rot2inert_state error");
    
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
    stateInert[0] = ( unit_x[0]*stateRot[0] + unit_y[0]*stateRot[1] + unit_z[0]*stateRot[2]);   // Inertial X
    stateInert[1] = ( unit_x[1]*stateRot[0] + unit_y[1]*stateRot[1] + unit_z[1]*stateRot[2]);   // Inertial Y
    stateInert[2] = ( unit_x[2]*stateRot[0] + unit_y[2]*stateRot[1] + unit_z[2]*stateRot[2]);   // Inertial Z
    stateInert[3] = ( unit_x[0]*stateRot[3] + unit_y[0]*stateRot[4] + unit_z[0]*stateRot[5] + pow(inst_charL, -2)*(angMom_p2[1]*stateInert[2] - angMom_p2[2]*stateInert[1]));   // Inertial VX
    stateInert[4] = ( unit_x[1]*stateRot[3] + unit_y[1]*stateRot[4] + unit_z[1]*stateRot[5] - pow(inst_charL, -2)*(angMom_p2[0]*stateInert[2] - angMom_p2[2]*stateInert[0]));   // Inertial VY
    stateInert[5] = ( unit_x[2]*stateRot[3] + unit_y[2]*stateRot[4] + unit_z[2]*stateRot[5] + pow(inst_charL, -2)*(angMom_p2[0]*stateInert[1] - angMom_p2[1]*stateInert[0]));   // Inertial VZ

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
 *  \ingroup bc4bp
 *  @brief Convert a Sun-Earth CR3BP Arcset to a Sun-Earth-Moon BCR4BPR Arcset.
 *  
 *  @param crNodes a CR3BP Sun-Earth arcset
 *  @param pBCSys a BCR4BPR Sun-Earth-Moon system data object; contains information about system
 *  scaling and orientation at time t = 0
 *  @param nodeID ID of a node in the set for which the epoch is known
 *  @param t0 the epoch at the specified node in BCR4BPR non-dimensional time units
 *
 *  @return a BCR4BPR arcset
 *  @throw Exception if `crNodes` is not a Sun-Earth trajectory or if `pBCSys` is
 *  not the Sun-Earth-Moon system.
 */
Arcset_bc4bp bcr4bpr_SE2SEM(Arcset_cr3bp crNodes, const SysData_bc4bp *pBCSys, int nodeID, double t0){
    if(crNodes.getSysData()->getPrimID(0) != 10 || crNodes.getSysData()->getPrimID(1) != 399){
        throw Exception("CR3BP trajectory is not in the Sun-Earth System");
    }

    if(pBCSys->getPrimID(0) != 10 || pBCSys->getPrimID(1) != 399 || pBCSys->getPrimID(2) != 301){
        throw Exception("BCR4BPR system is not Sun-Earth-Moon");
    }

    // Create a BCR4BPR Arcset
    Arcset_bc4bp bcNodes(pBCSys);

    double charL2 = crNodes.getSysData()->getCharL();
    double charT2 = crNodes.getSysData()->getCharT();
    double charL3 = pBCSys->getCharL();
    double charT3 = pBCSys->getCharT();
    unsigned int core_dim = crNodes.getSysData()->getDynamicsModel()->getCoreStateSize();
    unsigned int full_dim = (core_dim + 1)*core_dim + crNodes.getSysData()->getDynamicsModel()->getExtraStateSize();

    crNodes.updateEpochs(nodeID, t0*charT3/charT2);

    // A mapping vector: index is the old node ID, value is the new node ID
    // All new IDs are initialized to the default INVALID_ID value
    std::vector<int> map_oldID_to_newID(crNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> bcNodeState, crNodeState;
    double epoch;
    for(unsigned int n = 0; n < crNodes.getNumNodes(); n++){
        bcNodeState.clear();
        crNodeState = crNodes.getStateByIx(n);
        for(unsigned int r = 0; r < core_dim; r++){
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

    for(unsigned int s = 0; s < crNodes.getNumSegs(); s++){
        Segment seg = crNodes.getSegByIx(s);

        // Get a copy of the state and time vectors and modify them to the proper coordinates / dimensions
        std::vector<double> segStates = seg.getStateVector();
        std::vector<double> segTimes = seg.getTimeVector();

        if(segTimes.size() != segStates.size()/full_dim)
            throw Exception("cr3bp_rot2inert::Segment time and state vectors are not the expected sizes.");

        for(unsigned int i = 0; i < segTimes.size(); i++){
            std::vector<double> crSegState = seg.getStateByRow(i);

            for(unsigned int r = 0; r < core_dim; r++){
                if(r == 0)  // Convert x-coordinate, shift base to P2/P3 Barycenter
                    segStates[full_dim*i+r] = (crSegState[r]*charL2/charL3 - (1.0/pBCSys->getK() - pBCSys->getMu()));
                else if(r < 3)   // Convert position
                    segStates[full_dim*i+r] = (crSegState[r]*charL2/charL3);
                else  // Convert velocity
                    segStates[full_dim*i+r] = (crSegState[r]*(charL2/charL3)*(charT3/charT2));
            }

            segTimes[i] *= charT2/charT3;
        }

        seg.setStateVector(segStates);
        seg.setTimeVector(segTimes);
        seg.updateTOF();

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
 *  @brief Transform a BCR4BPR Sun-Earth-Moon arcset into a CR3BP Sun-Earth arcset
 *
 *  @param bcNodes a BCR4BPR Sun-Earth-Moon arcset
 *  @param pCRSys a CR3BP Sun-Earth system data object; contains information about system
 *  scaling
 *
 *  @return a CR3BP arcset object
 *  @throw Exception if `pCRSys` is not a Sun-Earth system or if `bcNodes` is
 *  not in the Sun-Earth-Moon system.
 */
Arcset_cr3bp bcr4bpr_SEM2SE(Arcset_bc4bp bcNodes, const SysData_cr3bp *pCRSys){
    if(pCRSys->getPrimID(0) != 10 || pCRSys->getPrimID(1) != 399){
        throw Exception("CR3BP system is not the Sun-Earth System");
    }

    if(bcNodes.getSysData()->getPrimID(0) != 10 || bcNodes.getSysData()->getPrimID(1) != 399 ||
            bcNodes.getSysData()->getPrimID(2) != 301){
        throw Exception("BCR4BPR trajectory is not in Sun-Earth-Moon");
    }

    const SysData_bc4bp *pBCSys = static_cast<const SysData_bc4bp *>(bcNodes.getSysData());

    unsigned int core_dim = pBCSys->getDynamicsModel()->getCoreStateSize();
    unsigned int full_dim = (core_dim + 1)*core_dim + pBCSys->getDynamicsModel()->getExtraStateSize();

    // Create a BCR4BPR Trajectory
    Arcset_cr3bp crNodes(pCRSys);

    double charL2 = pCRSys->getCharL();
    double charT2 = pCRSys->getCharT();
    double charL3 = bcNodes.getSysData()->getCharL();
    double charT3 = bcNodes.getSysData()->getCharT();

    std::vector<int> map_oldID_to_newID(bcNodes.getNextNodeID(), Linkable::INVALID_ID);
    std::vector<double> bcState, crNodeState;
    for(unsigned int n = 0; n < bcNodes.getNumNodes(); n++){
        bcState = bcNodes.getStateByIx(n);
        crNodeState.clear();

        for(unsigned int r = 0; r < core_dim; r++){
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
        node.setExtraParam(PARAMKEY_JACOBI, DynamicsModel_cr3bp::getJacobi(&(crNodeState[0]), pCRSys->getMu()));
        
        // Save the new ID
        map_oldID_to_newID[bcNodes.getNodeByIx(n).getID()] = crNodes.addNode(node);
    }

    for(unsigned int s = 0; s < bcNodes.getNumSegs(); s++){
        Segment seg = bcNodes.getSegByIx(s);


        // Get a copy of the state and time vectors and modify them to the proper coordinates / dimensions
        std::vector<double> segStates = seg.getStateVector();
        std::vector<double> segTimes = seg.getTimeVector();

        if(segTimes.size() != segStates.size()/full_dim)
            throw Exception("cr3bp_rot2inert::Segment time and state vectors are not the expected sizes.");

        for(unsigned int i = 0; i < segTimes.size(); i++){
            std::vector<double> bcSegState = seg.getStateByRow(i);

            for(int r = 0; r < 6; r++){
                if(r == 0)  // Shift origin to P2-P3 barycenter
                    segStates[full_dim*i+r] = ((bcSegState[r] + 1.0/pBCSys->getK() - pBCSys->getMu())*charL3/charL2);
                else if(r < 3)   // Convert position
                    segStates[full_dim*i+r] = (bcSegState[r]*charL3/charL2);
                else if(r < 6)  // Convert velocity
                    segStates[full_dim*i+r] = (bcSegState[r]*(charL3/charL2)*(charT2/charT3));
            }

            segTimes[i] *= charT3/charT2;
        }

        seg.setStateVector(segStates);
        seg.setTimeVector(segTimes);
        seg.updateTOF();

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
 *  @brief Compute the location of the saddle point for a specific bicircular system
 *  and epoch
 *  @details This function uses a Newton-Raphson procedure to locate the zeros of the local
 *  acceleration field.
 *  @throws DivergeException if the Newton-Raphson procedure cannot converge
 *  @param pBCSys system data object describing the bicircular system
 *  @param t0 the epoch at which the saddle point's location is computed
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
    pBCSys->getDynamicsModel()->getPrimPos(t0, pBCSys, -1, primPos);
    
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
 *  @brief Compute coefficients for 2nd-order polynomials in Epoch time that
 *  describe the x, y, and z coordinates of the saddle point.
 *  @details This function employs least squares to compute the coefficients. The
 *  number of points used and time span searched are hard coded in the function.
 * 
 *  @param pBCSys data about the bicircular system
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



