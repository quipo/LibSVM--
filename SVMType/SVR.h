/* 
 * File:   SVR.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_SVR_H
#define	SVM_SVMType_SVR_H

#include "../SVM_SVMType.h"
#include "../SVM_Solver.h"
#include "../interface_functions.h"

class SVM_SVMType_SVR : public SVM_SVMType {
public:
	SVM_SVMType_SVR(const SVM_Parameter & param) : SVM_SVMType(param) { };
	SVM_SVMType_SVR(const SVM_Problem & prob, const SVM_Parameter & param);
	~SVM_SVMType_SVR();

	std::vector<Qfloat> get_Q(unsigned int column, unsigned int len) const;
    void swap_index(unsigned int i, unsigned int j);
    double svm_predict(const SVM_Model & model, const std::vector<SVM_Node> & x) const;
    double svm_predict_probability(const SVM_Model & model, const std::vector<SVM_Node> & x, std::vector<double> & prob_estimates) const {
        if (prob_estimates.size() > 0) {
            std::cerr << "wrong SVM type for probability estimates";
        }
        return svm_predict(model, x);
    }
    SVM_Model svm_train(const SVM_Problem & prob, const SVM_Parameter & param) const;
    double svm_svr_probability(const SVM_Problem & prob, const SVM_Parameter & param) const;
    double svm_get_svr_probability(const SVM_Model & model) const;

    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;

protected:
	unsigned int l;
	Cache<unsigned int, std::vector<Qfloat> > * cache;
	std::vector<schar> sign;
	std::vector<int> index;
	mutable int next_buffer;
	std::vector<std::vector<Qfloat> > buffer; //dim = 2
};

#endif	/* SVM_SVMType_SVR_H */

