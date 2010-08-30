/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEL_LINEAR_H
#define	SVM_KERNEL_LINEAR_H

#include "../SVM_Kernel.h"

class SVM_Kernel_Linear : public SVM_Kernel {
public:
	SVM_Kernel_Linear(const SVM_Parameter & param) : SVM_Kernel(param) {};
	SVM_Kernel_Linear(const node_matrix & x, const SVM_Parameter & param);
	~SVM_Kernel_Linear();

	Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const;
    std::vector<double> get_diag(double c = 0) const;
};

#endif	/* SVM_KERNEL_LINEAR_H */

