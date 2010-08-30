/* 
 * File:   Kernel.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "NU_SVR.h"

void SVM_SVMType_NU_SVR::solve(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    std::vector<double> & alpha,
    SVM_SolutionInfo & si
) const {
    int l = prob.dimension;
    double C = param.C;
    std::vector<double> alpha2(2 * l, 0);
    std::vector<double> linear_term(2 * l, 1);
    std::vector<schar> y(2 * l, 1);
    int i;

    double sum = C * param.nu * l / 2;
    for (i = 0; i < l; i++) {
        alpha2[i] = alpha2[i + l] = min(sum, C);
        sum -= alpha2[i];

        linear_term[i]     = -prob.labels[i];
        linear_term[i + l] = prob.labels[i];
        y[i + l] = -1;
    }

    SVM_Solver_NU s;
    SVM_SVMType * svmtype = new SVM_SVMType_SVR(prob, param);
    s.solve(2 * l, svmtype, linear_term, y,
        alpha2, C, C, param.eps, si, param.shrinking);
    delete svmtype;

    std::cout << "epsilon = " << (-si.r) << "\n";

    for (i = 0; i < l; i++) {
        alpha.push_back(alpha2[i] - alpha2[i + l]);
    }

    alpha2.clear();
    linear_term.clear();
    y.clear();
}

SVM_SVMType_NU_SVR::~SVM_SVMType_NU_SVR() {
    delete cache;
    sign.clear();
    index.clear();
    buffer.clear();
    QD.clear();
}
