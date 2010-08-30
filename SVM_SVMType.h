/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_SVMType_H
#define	SVM_SVMType_H

typedef signed char schar;

#include "SVM_Node.h"
#include "SVM_Kernel.h"
#include "SVM_Model.h"
#include "SVM_Problem.h"
#include "SVM_Parameter.h"
#include "SVM_Solutioninfo.h"
#include "Cache.h"
#include "LRUCache.h"
#include "Kernel/Linear.h"
#include "Kernel/Polynomial.h"
#include "Kernel/Gaussian.h"
#include "Kernel/Sigmoid.h"
#include "Kernel/Stump.h"
#include "Kernel/Perceptron.h"
#include "Kernel/Laplace.h"
#include "Kernel/Exponential.h"
#include "Kernel/Precomputed.h"
#include <vector>

struct SVM_Exception {
    const char * p;
    SVM_Exception(const char * q) {
        p = q;
    }
};


class SVM_SVMType {
  public:
    SVM_SVMType(const SVM_Parameter & param) {
        switch (param.kernel_type) {
            case LINEAR:
                kernel = new SVM_Kernel_Linear(param);
                break;
            case POLYNOMIAL:
                kernel = new SVM_Kernel_Polynomial(param);
                break;
            case GAUSSIAN:
                kernel = new SVM_Kernel_Gaussian(param);
                break;
            case SIGMOID:
                kernel = new SVM_Kernel_Sigmoid(param);
                break;
            case STUMP:
                kernel = new SVM_Kernel_Stump(param);
                break;
            case PERCEPTRON:
                kernel = new SVM_Kernel_Perceptron(param);
                break;
            case LAPLACE:
                kernel = new SVM_Kernel_Laplace(param);
                break;
            case EXPONENTIAL:
                kernel = new SVM_Kernel_Exponential(param);
                break;
            case PRECOMPUTED:
                kernel = new SVM_Kernel_Precomputed(param);
                break;
        }
    };

    SVM_SVMType(const SVM_Problem & prob, const SVM_Parameter & param) {
    std::cout << "\nKERNEL: " << param.kernel_type;
    switch (param.kernel_type) {
        case LINEAR:
            std::cout << " (linear)";
            kernel = new SVM_Kernel_Linear(prob.x, param);
            break;
        case POLYNOMIAL:
            std::cout << " (polynomial)";
            kernel = new SVM_Kernel_Polynomial(prob.x, param);
            break;
        case GAUSSIAN:
            std::cout << " (gaussian)";
            kernel = new SVM_Kernel_Gaussian(prob.x, param);
            break;
        case SIGMOID:
            std::cout << " (sigmoid)";
            kernel = new SVM_Kernel_Sigmoid(prob.x, param);
            break;
        case STUMP:
            std::cout << " (stump)";
            kernel = new SVM_Kernel_Stump(prob.x, param);
            break;
        case PERCEPTRON:
            std::cout << " (perc)";
            kernel = new SVM_Kernel_Perceptron(prob.x, param);
            break;
        case LAPLACE:
            std::cout << " (laplace)";
            kernel = new SVM_Kernel_Laplace(prob.x, param);
            break;
        case EXPONENTIAL:
            std::cout << " (exponential)";
            kernel = new SVM_Kernel_Exponential(prob.x, param);
            break;
        case PRECOMPUTED:
            std::cout << " (precomputed)";
            kernel = new SVM_Kernel_Precomputed(prob.x, param);
            break;
        }
    };
    SVM_SVMType() { }
    
    virtual ~SVM_SVMType() {
        QD.clear();
        delete kernel;
    }
    
    virtual std::vector<Qfloat> get_Q(const unsigned int column, const unsigned int len) const = 0;
    std::vector<Qfloat> get_QD() const { return QD; }
    virtual double svm_predict(const SVM_Model & model, const std::vector<SVM_Node> & x) const = 0;
    virtual double svm_predict_probability(const SVM_Model & model, const std::vector<SVM_Node> & x, std::vector<double> & prob_estimates) const = 0;
    virtual SVM_Model svm_train(const SVM_Problem & prob, const SVM_Parameter & param) const = 0;
    virtual void solve(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<double> & alpha, SVM_SolutionInfo & si, double Cp, double Cn) const = 0;

    decision_function svm_train_one(
        const SVM_Problem & prob,
        const SVM_Parameter & param,
        double Cp,
        double Cn
    ) const {
        std::vector<double> alpha(prob.dimension, 0.0);
        SVM_SolutionInfo si;
        this->solve(prob, param, alpha, si, Cp, Cn);
        std::cout << "obj = " << si.obj << ", rho = " << si.rho << "\n";

        // output SVs

        unsigned int nSV = 0;
        unsigned int nBSV = 0;
        for (unsigned int i = 0; i < prob.dimension; i++) {
            if (fabs(alpha[i]) > 0) {
                ++nSV;
                // TODO check if the control "if (param->svm_type != C_SVC_L2 && param->svm_type != SVDD_L2) {"
                // is really useful
                //if (param->svm_type != C_SVC_L2 && param->svm_type != SVDD_L2) {
                if (prob.labels[i] > 0) {
                    if (fabs(alpha[i]) >= si.upper_bound_p) {
                        ++nBSV;
                    }
                } else {
                    if (fabs(alpha[i]) >= si.upper_bound_n) {
                        ++nBSV;
                    }
                }
                //}
            }
        }

        std::cout << "\nnSV = " << nSV << ", nBSV = " << nBSV;

        decision_function f;
        f.alpha.assign(alpha.begin(), alpha.end());
        f.rho.assign(1, si.rho); //only 1 element

        std::cout << "\nf.rho[0] = " << f.rho[0];
        
        f.dimension = nSV;
        f.lbsv      = nBSV;
        //f.BSV_idx = Malloc(unsigned int, nBSV);
        //f.SV_idx  = Malloc(unsigned int, nSV - nBSV);

        // TODO check if the control "if (param->svm_type != C_SVC_L2 && param->svm_type != SVDD_L2) {"
        // is needed also below

        for (unsigned int i = 0; i < prob.dimension; ++i) {
            if (fabs(alpha[i]) > 0) {
                if (prob.labels[i] > 0) {
                    if (fabs(alpha[i]) >= si.upper_bound_p) {
                        f.BSV_idx.push_back(i);
                    } else {
                        f.SV_idx.push_back(i);
                    }
                } else {
                    if (fabs(alpha[i]) >= si.upper_bound_n) {
                        f.BSV_idx.push_back(i);
                    } else {
                        f.SV_idx.push_back(i);
                    }
                }
            }
        }
        return f;
    }

    virtual void swap_index(unsigned int i, unsigned int j) = 0;
    
    std::vector<double> get_diag(double c = 0) const {
        return kernel->get_diag(c);
    }

    SVM_Kernel * kernel;

  protected:
    std::vector<Qfloat> QD;
    /* kernel_type */
    enum { LINEAR, POLYNOMIAL, GAUSSIAN, SIGMOID, STUMP, PERCEPTRON, LAPLACE, EXPONENTIAL, PRECOMPUTED };
};

#endif	/* SVM_SVMType_H */
