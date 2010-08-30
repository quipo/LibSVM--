/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_SVDD_L2_H
#define	SVM_SVMType_SVDD_L2_H

#include "./ONE_CLASS_L2.h"

class SVM_SVMType_SVDD_L2 : public SVM_SVMType_ONE_CLASS_L2 {
public:
    SVM_SVMType_SVDD_L2(const SVM_Parameter & param) : SVM_SVMType_ONE_CLASS_L2(param) { }
	SVM_SVMType_SVDD_L2(const SVM_Problem & prob, const SVM_Parameter & param) : SVM_SVMType_ONE_CLASS_L2(prob, param) { }
	~SVM_SVMType_SVDD_L2();

    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;
};

#endif	/* SVM_SVMType_SVDD_L2_H */

