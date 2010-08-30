/* 
 * File:   Kernel.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Laplace.h"

SVM_Kernel_Laplace::SVM_Kernel_Laplace(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
}

SVM_Kernel_Laplace::~SVM_Kernel_Laplace() {
    x_square.clear();
}

Qfloat SVM_Kernel_Laplace::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    return exp(-param.gamma * dist_1(x, y));
}
