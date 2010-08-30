/* 
 * File:   SVR.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "SVR.h"

SVM_SVMType_SVR::SVM_SVMType_SVR(const SVM_Problem & prob, const SVM_Parameter & param)
  : SVM_SVMType(prob, param), l(prob.dimension)
{
    cache = new LRUCache<unsigned int, std::vector<Qfloat> >(static_cast<size_t>(param.cache_size*(1<<20)));
    QD.assign(2 * l, 0);
    sign.assign(2 * l, 0);
    index.assign(2 * l, 0);
    for (unsigned int k=0; k<l; ++k) {
        sign[k]    = 1;
        sign[k+l]  = -1;
        index[k]   = k;
        index[k+l] = k;
        QD[k]   = kernel->kernel(k, k);
        QD[k+l] = QD[k];
    }
    std::vector<Qfloat> a(2 * l);
    std::vector<Qfloat> b(2 * l);
    buffer.push_back(a);
    buffer.push_back(b);
    next_buffer = 0;
}

std::vector<Qfloat> SVM_SVMType_SVR::get_Q(unsigned int column, unsigned int len) const {
    std::vector<Qfloat> data;
    unsigned int real_i = index[column];
    //unsigned int start = (cache->find(real_i, data)) ? data.size() : 0;
    unsigned int start = 0;
    if (start < len) {
        data.resize(len, 0);
        for (unsigned int j=start; j<len; ++j) {
            data[j] = kernel->kernel(real_i, j);
        }
        //cache->insert(column, data);
    }

    // reorder and copy
    std::vector<Qfloat> buf(buffer[next_buffer]);
    next_buffer = 1 - next_buffer;
    schar si = sign[column];
    for (unsigned int j=0; j<len; j++) {
        buf.push_back(static_cast<Qfloat>(si) * static_cast<Qfloat>(sign[j]) * data[index[j]]);
    }
    return buf;
}

void SVM_SVMType_SVR::swap_index(unsigned int i, unsigned int j) {
    std::swap(sign[i],  sign[j]);
    std::swap(index[i], index[j]);
    std::swap(QD[i], QD[j]);
}

double SVM_SVMType_SVR::svm_predict(const SVM_Model & model, const std::vector<SVM_Node> & x) const {
    double sum = 0;
    for (unsigned int i = 0; i < model.dimension; ++i) {
        sum += model.sv_coef[0][i] * kernel->kernel(x, model.SV[i], model.param);
    }
    sum -= model.rho[0];
    return sum;
}

SVM_Model SVM_SVMType_SVR::svm_train(const SVM_Problem & prob, const SVM_Parameter & param) const {
    SVM_Model model(param, 1u);
    model.free_sv = 0; // XXX
    // regression or one-class-svm or SVDD
    model.nr_class = 2;
    model.label.clear();
    model.nSV.clear();
    model.probA.clear();
    model.probB.clear();
    //model.sv_coef = Malloc(double *, 1);

    if (param.probability) {
        //model.probA = Malloc(double, 1);
        model.probA.push_back(svm_svr_probability(prob, param));
    }
    //std::cout << "\n\nTraining one...\n\n";
    decision_function f = svm_train_one(prob, param, 0, 0);
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

double SVM_SVMType_SVR::svm_get_svr_probability(const SVM_Model & model) const {
    if (model.probA.size() > 0) {
        return model.probA[0];
    }
    std::cerr << "Model doesn't contain information for SVR probability inference\n";
    return 0;
}

/**
 * Return parameter of a Laplace distribution
 */
double SVM_SVMType_SVR::svm_svr_probability(const SVM_Problem & prob, const SVM_Parameter & param) const {
    unsigned int i;
    unsigned int nr_fold = 5;
    std::vector<double> ymv(prob.dimension);
    double mae = 0;

    SVM_Parameter newparam = param;
    newparam.probability = 0;
    svm_cross_validation(prob, newparam, nr_fold, ymv);
    for (i = 0; i < prob.dimension; i++) {
        ymv[i] = prob.labels[i] - ymv[i];
        mae += fabs(ymv[i]);
    }
    mae /= prob.dimension;
    double std = sqrt(2 * mae * mae);
    int count = 0;
    mae = 0;
    double threshold = 5 * std;
    for (unsigned i = 0; i < prob.dimension; ++i) {
        double fabsi = fabs(ymv[i]);
        if (fabsi > threshold) {
            ++count;
        } else {
            mae += fabsi;
        }
    }
    mae /= (prob.dimension - count);
    std::cout << "Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma= " << mae << "\n";
    ymv.clear();
    return mae;
}

void SVM_SVMType_SVR::solve(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    std::vector<double> & alpha,
    SVM_SolutionInfo & si,
    double Cp,
    double Cn
) const {
    int l = prob.dimension;
    std::vector<double> alpha2(2 * l, 0);
    std::vector<double> linear_term(2 * l, param.p);
    std::vector<schar> y(2 * l, 1);
    int i;

    if (Cp < Cn) {
        // ugly hack to avoid unused variable warning
    }

    for (i = 0; i < l; i++) {
        linear_term[i]     -= prob.labels[i];
        linear_term[i + l] += prob.labels[i];
        y[i + l] = -1;
    }

    SVM_Solver s;
    SVM_SVMType * svmtype = new SVM_SVMType_SVR(prob, param);
    s.solve(2 * l, svmtype, linear_term, y,
        alpha2, param.C, param.C, param.eps, si, param.shrinking);
    delete svmtype;

    double sum_alpha = 0, a;
    for (i = 0; i < l; i++) {
        a = alpha2[i] - alpha2[i + l];
        sum_alpha += fabs(a);
        alpha.push_back(a);
    }
    //std::cout << "nu = " << (sum_alpha / (param.C * l)) << "\n";

    alpha2.clear();
    linear_term.clear();
    y.clear();
}

SVM_SVMType_SVR::~SVM_SVMType_SVR() {
    delete cache;
    sign.clear();
    index.clear();
    buffer.clear();
    QD.clear();
}
