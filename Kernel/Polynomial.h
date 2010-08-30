/* 
 * File:   Polynomial.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEN_POLYNOMIAL_H
#define	SVM_KERNEN_POLYNOMIAL_H

#include "../SVM_Kernel.h"

static inline double powi(double base, int times) {
    double tmp = base, ret = 1.0;
    for (int t=times; t>0; t/=2) {
        if (t%2 == 1) {
            ret *= tmp;
        }
        tmp *= tmp;
    }
    return ret;
}

class SVM_Kernel_Polynomial : public SVM_Kernel {
  public:
    SVM_Kernel_Polynomial(const SVM_Parameter & param) : SVM_Kernel(param) {};
    SVM_Kernel_Polynomial(const node_matrix & x, const SVM_Parameter & param);
    ~SVM_Kernel_Polynomial();

    Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const;
    std::vector<double> get_diag(double c = 0) const;
};

#endif	/* SVM_KERNEN_POLYNOMIAL_H */

