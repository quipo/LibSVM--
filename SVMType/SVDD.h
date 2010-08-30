/**
 * Support Vector Domain Description (SVDD) 
 *  
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_SVDD_H
#define	SVM_SVMType_SVDD_H

#include "./ONE_CLASS.h"

class SVM_SVMType_SVDD : public SVM_SVMType_ONE_CLASS {
public:
    SVM_SVMType_SVDD(const SVM_Parameter & param) : SVM_SVMType_ONE_CLASS(param) { };
	SVM_SVMType_SVDD(const SVM_Problem & prob, const SVM_Parameter & param) : SVM_SVMType_ONE_CLASS(prob, param) { }
	~SVM_SVMType_SVDD();

    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;
};

#endif	/* SVM_SVMType_SVDD_H */

