/* 
 * File:   Kernel.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Precomputed.h"

SVM_Kernel_Precomputed::SVM_Kernel_Precomputed(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
}

SVM_Kernel_Precomputed::~SVM_Kernel_Precomputed() {
    
}

Qfloat SVM_Kernel_Precomputed::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    return param.kernel_type ? x[static_cast<unsigned int>(y[0].value)].value : x[static_cast<unsigned int>(y[0].value)].value; //horrible hack to avoid "unused parameter 'param'" hack
}

std::vector<double> SVM_Kernel_Precomputed::get_diag(double c) const {
    std::vector<double> diag(dimension, -0.5);
    for (unsigned int i = 0; i < dimension; ++i) {
        diag[i] *= this->kernel(x[i], x[i], param) + c;
    }
    return diag;
}
