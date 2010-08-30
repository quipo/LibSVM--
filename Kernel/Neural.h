/* 
 * File:   Neural.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEL_NEURAL_H
#define	SVM_KERNEL_NEURAL_H

#include "../SVM_Kernel.h"

class SVM_Kernel_Neural : public SVM_Kernel {
public:
    SVM_Kernel_Neural(const SVM_Parameter & param) : SVM_Kernel(param) {};
    SVM_Kernel_Neural(const node_matrix & x, const SVM_Parameter & param);
    ~SVM_Kernel_Neural();

    Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
        return tanh(param.gamma * dot(x, y) + param.coef0);
    }
};

#endif	/* SVM_KERNEL_NEURAL_H */

