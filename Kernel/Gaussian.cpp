/* 
 * File:   Gaussian.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Gaussian.h"

SVM_Kernel_Gaussian::SVM_Kernel_Gaussian(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    x_square.reserve(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        x_square.push_back(dot(x[i], x[i]));
    }
}

SVM_Kernel_Gaussian::~SVM_Kernel_Gaussian() {
    x_square.clear();
}

//Qfloat SVM_Kernel_Gaussian::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
//    //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(1) << dist_2_sqr(x, y) << " ";
//    return exp(-param.gamma * dist_2_sqr(x, y));
//}
