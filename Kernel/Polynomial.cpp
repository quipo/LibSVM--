/* 
 * File:   Polynomial.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Polynomial.h"

SVM_Kernel_Polynomial::SVM_Kernel_Polynomial(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
}

SVM_Kernel_Polynomial::~SVM_Kernel_Polynomial() {

}

Qfloat SVM_Kernel_Polynomial::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    return powi(param.gamma * dot(x, y) + param.coef0, param.degree);
}

std::vector<double> SVM_Kernel_Polynomial::get_diag(double c) const {
    std::vector<double> diag(dimension, -0.5);
    for (unsigned int i = 0; i < dimension; ++i) {
        diag[i] *= this->kernel(x[i], x[i], param) + c;
    }
    return diag;
}
