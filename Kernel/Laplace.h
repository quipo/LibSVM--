/* 
 * File:   Laplace.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEL_LAPLACE_H
#define	SVM_KERNEL_LAPLACE_H

#include "../SVM_Kernel.h"

class SVM_Kernel_Laplace : public SVM_Kernel {
public:
	SVM_Kernel_Laplace(const SVM_Parameter & param) : SVM_Kernel(param) {};
	SVM_Kernel_Laplace(const node_matrix & x, const SVM_Parameter & param);
	~SVM_Kernel_Laplace();

    Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const;
};

#endif	/* SVM_KERNEL_LAPLACE_H */

