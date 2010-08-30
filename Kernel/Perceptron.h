/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEL_PERCEPTRON_H
#define	SVM_KERNEL_PERCEPTRON_H

#include "../SVM_Kernel.h"

class SVM_Kernel_Perceptron : public SVM_Kernel {
public:
	SVM_Kernel_Perceptron(const SVM_Parameter & param) : SVM_Kernel(param) {};
	SVM_Kernel_Perceptron(const node_matrix & x, const SVM_Parameter & param);
	~SVM_Kernel_Perceptron();

	Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const;
    std::vector<double> get_diag(double c = 0) const;
};

#endif	/* SVM_KERNEL_PERCEPTRON_H */

