/**
 * Two layered neural net: tanh(a x*y + b)
 *  
 * File:   Neural.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "Neural.h"

SVM_Kernel_Neural::SVM_Kernel_Neural(const node_matrix & x, const SVM_Parameter & param)
  : SVM_Kernel(x, param)
{
    //clone(x_, x);
    x_square.assign(dimension, 0);
    for (unsigned int i = 0; i < dimension; ++i) {
        x_square.push_back(dot(x[i], x[i]));
    }
}

SVM_Kernel_Neural::~SVM_Kernel_Neural() {
    x_square.clear();
}
