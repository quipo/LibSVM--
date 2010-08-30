/* 
 * File:   SVC_L2.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_SVC_L2_H
#define	SVM_SVMType_SVC_L2_H

#include "./SVC.h"

class SVM_SVMType_SVC_L2 : public SVM_SVMType_SVC {
public:
	SVM_SVMType_SVC_L2(const SVM_Parameter & param) : SVM_SVMType_SVC(param) { };
	SVM_SVMType_SVC_L2(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<schar> & y_) : SVM_SVMType_SVC(prob, param, y_) {
        C = param.C;
    }
	~SVM_SVMType_SVC_L2();

	std::vector<Qfloat> get_Q(unsigned int column, unsigned int len) const;
    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;
    
private:
	std::vector<schar> y;
	Cache<unsigned int, std::vector<Qfloat> > * cache;
	double C; /* cost for C_SVC, SVDD, EPSILON_SVR and NU_SVR */
};

#endif	/* SVM_SVMType_SVC_L2_H */
