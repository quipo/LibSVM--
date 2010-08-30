/* 
 * File:   NU_SVC.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "NU_SVC.h"

void SVM_SVMType_NU_SVC::solve(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    std::vector<double> & alpha,
    SVM_SolutionInfo & si
) const {
    unsigned int i;
    unsigned int l = prob.dimension;
    double nu = param.nu;

    std::vector<schar> y(l);
    for (i = 0; i < l; i++) {
        if (prob.labels[i] > 0) {
            y.push_back(+1);
        } else {
            y.push_back(-1);
        }
    }

    double sum_pos = nu * l / 2;
    double sum_neg = sum_pos;

    double a;
    for (i = 0; i < l; i++) {
        if (y[i] == +1) {
            a = min(1.0, sum_pos);
            sum_pos -= a;
        } else {
            a = min(1.0, sum_neg);
            sum_neg -= a;
        }
        alpha.push_back(a);
    }

    std::vector<double> zeros(l, 0);

    SVM_Solver_NU s;
    SVM_SVMType * svmtype = new SVM_SVMType_NU_SVC(prob, param, y);
    s.solve(l, svmtype, zeros, y,
        alpha, 1.0, 1.0, param.eps, si, param.shrinking);
    delete svmtype;

    double r = si.r;

    std::cout << "C = " << (1 / r) << "\n";

    for (i = 0; i < l; i++) {
        alpha[i] *= y[i] / r;
    }

    si.rho /= r;
    si.obj /= (r * r);
    si.upper_bound_p = 1 / r;
    si.upper_bound_n = 1 / r;

    y.clear();
    zeros.clear();
}


SVM_SVMType_NU_SVC::~SVM_SVMType_NU_SVC() {
    y.clear();
    delete cache;
    QD.clear();
}
