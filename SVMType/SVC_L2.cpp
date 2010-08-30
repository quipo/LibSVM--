/* 
 * File:   SVC_L2.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "SVC_L2.h"

std::vector<Qfloat> SVM_SVMType_SVC_L2::get_Q(unsigned int column, unsigned int len) const {
    std::vector<Qfloat> data(len, 0);
    //unsigned int start = (cache->find(column, data)) ? data.size() : 0;
    unsigned int start = 0;
    if (start < len) {
        data.resize(len, 0);
        for (unsigned int j=start; j<len; ++j) {
            data[j] = y[column] * y[j] * kernel->kernel(column, j);
        }
        if (column >= start && column < len) {
            data[column] += 1/C;
        }
        //cache->insert(column, data);
    }
    return data;
}

void SVM_SVMType_SVC_L2::solve(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    std::vector<double> & alpha,
    SVM_SolutionInfo & si,
    double Cp,
    double Cn
) const {
    unsigned int l = prob.dimension;
    std::vector<double> minus_ones(l, -1);
    std::vector<schar> labels(l, 1);
    unsigned int i;

    for (i = 0; i < l; i++) {
        alpha.push_back(0);
        if (prob.labels[i] > 0) {
            labels[i] = +1;
        } else {
            labels[i] = -1;
        }
    }

    SVM_Solver s;
    SVM_SVMType * svmtype = new SVM_SVMType_SVC_L2(prob, param, labels);
    s.solve(l, svmtype, minus_ones, labels,
        alpha, Cp, Cn, param.eps, si, param.shrinking);
    delete svmtype;

    double sum_alpha = 0;
    for (i = 0; i < l; i++) {
        sum_alpha += alpha[i];
    }

    if (Cp == Cn) {
        std::cout << "nu = " << (sum_alpha / (Cp * prob.dimension)) << "\n";
    }

    for (i = 0; i < l; i++) {
        alpha[i] *= labels[i];
    }

    minus_ones.clear();
    labels.clear();
}


SVM_SVMType_SVC_L2::~SVM_SVMType_SVC_L2() {
    y.clear();
    QD.clear();
    delete cache;
}
