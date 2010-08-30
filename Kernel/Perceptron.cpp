/* 
 * File:   Perceptron.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Perceptron.h"

SVM_Kernel_Perceptron::SVM_Kernel_Perceptron(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
    x_square.assign(dimension, 0);
    for (unsigned int i = 0; i < dimension; ++i) {
        x_square.push_back(dot(x[i], x[i]));
    }
}

SVM_Kernel_Perceptron::~SVM_Kernel_Perceptron() {
    x_square.clear();
}

Qfloat SVM_Kernel_Perceptron::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    return -sqrt(dist_2_sqr(x, y)) + param.coef0;
}

std::vector<double> SVM_Kernel_Perceptron::get_diag(double c) const {
    std::vector<double> diag(dimension, -0.5);
    for (unsigned int i = 0; i < dimension; ++i) {
        diag[i] *= this->kernel(x[i], x[i], param) + c;
    }
    return diag;
}
