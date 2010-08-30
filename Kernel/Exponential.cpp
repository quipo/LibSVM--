/* 
 * File:   Kernel.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Exponential.h"

SVM_Kernel_Exponential::SVM_Kernel_Exponential(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
    x_square.assign(dimension, 0);
    for (unsigned int i = 0; i < dimension; ++i) {
        x_square.push_back(dot(x[i], x[i]));
    }
}

SVM_Kernel_Exponential::~SVM_Kernel_Exponential() {
    x_square.clear();
}

Qfloat SVM_Kernel_Exponential::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    return exp(-param.gamma * sqrt(dist_2_sqr(x, y)));
}
