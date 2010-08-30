/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_NU_SVC_H
#define	SVM_SVMType_NU_SVC_H

#include "./SVC.h"
#include "../Solver/NU.h"

class SVM_SVMType_NU_SVC : public SVM_SVMType_SVC {
public:
	SVM_SVMType_NU_SVC(const SVM_Parameter & param) : SVM_SVMType_SVC(param) { };
	SVM_SVMType_NU_SVC(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<schar> & y) : SVM_SVMType_SVC(prob, param, y) { };
	~SVM_SVMType_NU_SVC();

    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;
};

#endif	/* SVM_SVMType_NU_SVC_H */

