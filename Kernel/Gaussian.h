/* 
 * File:   Gaussian.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEL_GAUSSIAN_H
#define	SVM_KERNEL_GAUSSIAN_H

#include "../SVM_Kernel.h"
#include <math.h>

class SVM_Kernel_Gaussian : public SVM_Kernel {
public:
	SVM_Kernel_Gaussian(const SVM_Parameter & param) : SVM_Kernel(param) {};
	SVM_Kernel_Gaussian(const node_matrix & x, const SVM_Parameter & param);
	~SVM_Kernel_Gaussian();

	Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
        //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(1) << exp(-param.gamma * dist_2_sqr(x, y)) << " ";
        return exp(-param.gamma * dist_2_sqr(x, y));
    }
};

#endif	/* SVM_KERNEL_GAUSSIAN_H */

