/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_NU_SVR_H
#define	SVM_SVMType_NU_SVR_H

#include "./SVR.h"
#include "../Solver/NU.h"

class SVM_SVMType_NU_SVR : public SVM_SVMType_SVR {
public:
    SVM_SVMType_NU_SVR(const SVM_Parameter & param) : SVM_SVMType_SVR(param) { }
    SVM_SVMType_NU_SVR(const SVM_Problem & prob, const SVM_Parameter & param) : SVM_SVMType_SVR(prob, param) { }
    ~SVM_SVMType_NU_SVR();

    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;
};

#endif	/* SVM_SVMType_NU_SVR_H */
