/** 
 * File:   SVDD.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "SVDD.h"

void SVM_SVMType_SVDD::solve(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    std::vector<double> & alpha,
    SVM_SolutionInfo & si
) const {
    unsigned int l = prob.dimension;
    std::vector<schar> ones(l, 1);
    unsigned int i;

    unsigned int n = (unsigned int) (param.nu * prob.dimension); // # of alpha's at upper bound

    for (i = 0; i < n; i++) {
        alpha.push_back(1);
    }
    if (n < prob.dimension) {
        alpha.push_back(param.nu * prob.dimension - n); // alpha[n]
    }
    for (i = n + 1; i < l; i++) {
        alpha.push_back(0);
    }

    SVM_SVMType * svmtype = new SVM_SVMType_ONE_CLASS(prob, param);
    // diag depends on the kernel: default to [-0.5] if kernel != (Gaussian|Exponential|Laplace)
    // otherwise diag[i] *= Kernel::k_function(prob.x[i], prob.x[i], *param);
    std::vector<double> diag = svmtype->get_diag();

    SVM_Solver s;
    s.solve(l, svmtype, diag, ones,
        alpha, 1.0, 1.0, param.eps, si, param.shrinking);
    delete svmtype;
    diag.clear();
    ones.clear();
}


SVM_SVMType_SVDD::~SVM_SVMType_SVDD() {
    delete cache;
    QD.clear();
}