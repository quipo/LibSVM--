/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_ONE_CLASS_H
#define	SVM_SVMType_ONE_CLASS_H

#include "../SVM_SVMType.h"
#include "../SVM_Solver.h"

class SVM_SVMType_ONE_CLASS : public SVM_SVMType {
public:
	SVM_SVMType_ONE_CLASS(const SVM_Parameter & param) : SVM_SVMType(param) { };
	SVM_SVMType_ONE_CLASS(const SVM_Problem & prob, const SVM_Parameter & param);
	~SVM_SVMType_ONE_CLASS();

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

    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;

protected:
	Cache<unsigned int, std::vector<Qfloat> > * cache;
};

#endif	/* SVM_SVMType_ONE_CLASS_H */

