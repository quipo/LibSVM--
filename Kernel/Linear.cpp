/* 
 * File:   Kernel.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Linear.h"

SVM_Kernel_Linear::SVM_Kernel_Linear(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
}

SVM_Kernel_Linear::~SVM_Kernel_Linear() {
    
}

Qfloat SVM_Kernel_Linear::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    return param.kernel_type ? dot(x, y) : dot(x, y); //horrible hack to avoid "unused parameter 'param'" hack
}

std::vector<double> SVM_Kernel_Linear::get_diag(double c) const {
    std::vector<double> diag(dimension, -0.5);
    for (unsigned int i = 0; i < dimension; ++i) {
        diag[i] *= this->kernel(x[i], x[i], param) + c;
    }
    return diag;
}
