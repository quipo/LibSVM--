/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_ONE_CLASS_L2_H
#define	SVM_SVMType_ONE_CLASS_L2_H

#include "./ONE_CLASS.h"

class SVM_SVMType_ONE_CLASS_L2 : public SVM_SVMType_ONE_CLASS {
public:
    SVM_SVMType_ONE_CLASS_L2(const SVM_Parameter & param) : SVM_SVMType_ONE_CLASS(param) { };
	SVM_SVMType_ONE_CLASS_L2(const SVM_Problem & prob, const SVM_Parameter & param) : SVM_SVMType_ONE_CLASS(prob, param) {
        C = param.C;
    }
	~SVM_SVMType_ONE_CLASS_L2();

	std::vector<Qfloat> get_Q(unsigned int column, unsigned int len) const;
	double svm_predict_probability(const SVM_Model & model, const std::vector<SVM_Node> & x, std::vector<double> & prob_estimates) const {
        if (prob_estimates.size() > 0) {
            std::cerr << "wrong SVM type for probability estimates";
        }
        return svm_predict(model, x);
    }

protected:
	double C;
};

#endif	/* SVM_SVMType_ONE_CLASS_L2_H */

