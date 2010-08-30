/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEL_EXPONENTIAL_H
#define	SVM_KERNEL_EXPONENTIAL_H

#include "../SVM_Kernel.h"

class SVM_Kernel_Exponential : public SVM_Kernel {
public:
    SVM_Kernel_Exponential(const SVM_Parameter & param) : SVM_Kernel(param) {};
	SVM_Kernel_Exponential(const node_matrix & x, const SVM_Parameter & param);
	~SVM_Kernel_Exponential();

	Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const;
};

#endif	/* SVM_KERNEL_EXPONENTIAL_H */

