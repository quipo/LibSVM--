/* 
 * File:   SVC.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_SVC_H
#define	SVM_SVMType_SVC_H

#include "../SVM_SVMType.h"
#include "../SVM_Solver.h"
#include "../interface_functions.h"

class SVM_SVMType_SVC : public SVM_SVMType {
public:
    SVM_SVMType_SVC(const SVM_Parameter & param) : SVM_SVMType(param) { };
    SVM_SVMType_SVC(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<schar> & y);
    ~SVM_SVMType_SVC();

    std::vector<Qfloat> get_Q(unsigned int column, unsigned int len) const;
    void swap_index(unsigned int i, unsigned int j);
    std::vector<double> svm_predict_values(const SVM_Model & model, const std::vector<SVM_Node> & x) const;
    double svm_predict(const SVM_Model & model, const std::vector<SVM_Node> & x) const;
    double svm_predict_probability(const SVM_Model & model, const std::vector<SVM_Node> & x, std::vector<double> & prob_estimates) const;
    void multiclass_probability(int k, std::vector<std::vector<double> > & pairwise_prob, std::vector<double> & prob_estimates) const;
    double sigmoid_predict(double decision_value, double A, double B) const;
    SVM_Model svm_train(const SVM_Problem & prob, const SVM_Parameter & param) const;
    void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const;

    /**
     * Cross-validation decision values for probability estimates
     */
    void svm_binary_svc_probability(const SVM_Problem & prob, const SVM_Parameter & param, double Cp, double Cn, double& probA, double& probB) const;

    /**
     * Platt's binary SVM Probablistic Output: an improvement from Lin et al.
     */
    void sigmoid_train(int l, const std::vector<double> & dec_values, const std::vector<double> & labels, double & A, double & B) const;

protected:
    std::vector<schar> y;
    Cache<unsigned int, std::vector<Qfloat> > * cache;
};

#endif	/* SVM_SVMType_SVC_H */

