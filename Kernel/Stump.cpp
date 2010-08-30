/* 
 * File:   Stump.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Stump.h"

SVM_Kernel_Stump::SVM_Kernel_Stump(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
}

SVM_Kernel_Stump::~SVM_Kernel_Stump() {
    
}

Qfloat SVM_Kernel_Stump::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    return -dist_1(x, y) + param.coef0;
}

std::vector<double> SVM_Kernel_Stump::get_diag(double c) const {
    std::vector<double> diag(dimension, -0.5);
    for (unsigned int i = 0; i < dimension; ++i) {
        diag[i] *= this->kernel(x[i], x[i], param) + c;
    }
    return diag;
}
