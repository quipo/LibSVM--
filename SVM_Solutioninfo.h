/* 
 * File:   SVM_SolutionInfo.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 5:30 PM
 */

#ifndef SVM_SOLUTIONINFO_H
#define	SVM_SOLUTIONINFO_H

struct SVM_SolutionInfo {
    double obj;
    double rho;
    double upper_bound_p;
    double upper_bound_n;
    double r;	// for Solver_NU
};

#endif	/* SVM_SOLUTIONINFO_H */
