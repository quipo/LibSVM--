/* 
 * File:   ONE_CLASS.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "ONE_CLASS.h"

SVM_SVMType_ONE_CLASS::SVM_SVMType_ONE_CLASS(const SVM_Problem & prob, const SVM_Parameter & param)
  : SVM_SVMType(prob, param)
{
    cache = new LRUCache<unsigned int, std::vector<Qfloat> >(static_cast<size_t>(param.cache_size*(1<<20)));
    QD.reserve(prob.dimension);
    for (unsigned int i=0; i<prob.dimension; i++) {
        QD.push_back(kernel->kernel(i, i));
    }
}

std::vector<Qfloat> SVM_SVMType_ONE_CLASS::get_Q(unsigned int column, unsigned int len) const {
    std::vector<Qfloat> data;
    data.assign(len, 0);
    //unsigned int start = (cache->find(column, data)) ? data.size() : 0;
    unsigned int start = 0;
    if (start < len) {
        data.resize(len, 0);
        for (unsigned int j=start; j<len; ++j) {
            data[j] = kernel->kernel(column, j);
        }
        //cache->insert(column, data);
    }
    return data;
}

void SVM_SVMType_ONE_CLASS::swap_index(unsigned int i, unsigned int j) {
    cache->swap_index(i, j);
    kernel->swap_index(i, j);
    std::swap(QD[i], QD[j]);
}

double SVM_SVMType_ONE_CLASS::svm_predict(const SVM_Model & model, const std::vector<SVM_Node> & x) const {
    double sum = 0;
    for (unsigned int i = 0; i < model.dimension; ++i) {
        sum += model.sv_coef[0][i] * kernel->kernel(x, model.SV[i], model.param);
    }
    sum -= model.rho[0];
    return (sum > 0) ? 1 : -1;
}

SVM_Model SVM_SVMType_ONE_CLASS::svm_train(const SVM_Problem & prob, const SVM_Parameter & param) const {
    SVM_Model model(param, 1u);
    model.free_sv = 0; // XXX
    // regression or one-class-svm or SVDD
    model.nr_class = 2;
    model.label.clear();
    model.nSV.clear();
    model.probA.clear();
    model.probB.clear();
    //model.sv_coef = Malloc(double *, 1);

    //std::cout << "\n\nTraining one...\n\n";
    decision_function f = svm_train_one(prob, param, 0.0, 0.0);
    //std::cout << "\n\nTrained One.\n\n";
    //model.rho = Malloc(double, 1);
    model.rho.push_back(f.rho[0]);
    //model.rho[0] = f.rho;

    model.dimension = f.dimension; //nSV
    model.lbsv = f.lbsv; //nBSV

    model.SV_idx  = f.SV_idx;
    model.BSV_idx = f.BSV_idx;

    //model.SV = Malloc(SVM_Node , f.l);
    //model.sv_coef[0] = Malloc(double, f.l);
    std::vector<double> coeff;
    for (unsigned int i = 0; i < prob.dimension; ++i) {
        if (fabs(f.alpha[i]) > 0) {
            model.SV.push_back(prob.x[i]);
            coeff.push_back(f.alpha[i]);
        }
    }
    model.sv_coef.assign(1, coeff);

    f.alpha.clear();
    return model;
}

void SVM_SVMType_ONE_CLASS::solve(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    std::vector<double> & alpha,
    SVM_SolutionInfo & si,
    double Cp,
    double Cn
) const {
    std::vector<double> zeros(prob.dimension, 0);
    std::vector<schar> ones(prob.dimension, 1);
    
    if (Cp > Cn) {
        // ugly hack to avoid unused variable compiler warnings
    }

    unsigned int n = static_cast<unsigned int>(param.nu * prob.dimension); // # of alpha's at upper bound

    //std::cout << "\nparam->nu, prob->l, n = " << param.nu << ", " << prob.dimension << ", " << n;

    alpha.assign(n, 1);
    if (n < prob.dimension) {
        alpha.push_back(param.nu * prob.dimension - n); // alpha[n]
    }
    alpha.resize(prob.dimension, 0);

    SVM_Solver s;
    SVM_SVMType * svmtype = new SVM_SVMType_ONE_CLASS(prob, param);
    s.solve(prob.dimension, svmtype, zeros, ones,
        alpha, 1.0, 1.0, param.eps, si, param.shrinking);
    delete svmtype;

    zeros.clear();
    ones.clear();
}

SVM_SVMType_ONE_CLASS::~SVM_SVMType_ONE_CLASS() {
    delete cache;
    QD.clear();
}
