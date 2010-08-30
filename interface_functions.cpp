/* 
 * File:   interface_functions.cpp
 * Author: lorenzo
 * 
 * Created on July 5, 2010, 2:00 PM
 */

#include <map>

#include "interface_functions.h"
#include "SVMType/NU_SVC.h"
#include "SVMType/NU_SVR.h"
#include "SVMType/ONE_CLASS_L2.h"
#include "SVMType/ONE_CLASS.h"
#include "SVMType/SVC_L2.h"
#include "SVMType/SVC.h"
#include "SVMType/SVDD.h"
#include "SVMType/SVDD_L2.h"
#include "SVMType/SVR.h"

SVM_SVMType * svmtype_factory(const SVM_Parameter & param) {
     switch (param.svm_type) {
        case C_SVC:
        case NU_SVC:
            return new SVM_SVMType_SVC(param);
        case C_SVC_L2:
            return new SVM_SVMType_SVC_L2(param);
        case ONE_CLASS:
        case SVDD:
        case SVDD_L2:
            return new SVM_SVMType_ONE_CLASS(param);
        case EPSILON_SVR:
        case NU_SVR:
            return new SVM_SVMType_SVR(param);
        default:
            throw SVM_Exception("Unknown svm_type param");
    }
}

SVM_SVMType * svmtype_factory(const SVM_Problem & prob, const SVM_Parameter & param) {
    std::vector<schar> ones(prob.dimension, 1);

    switch (param.svm_type) {
        case C_SVC:
        case NU_SVC:
            std::cout << "\nSVM_SVMType_SVC";
            return new SVM_SVMType_SVC(prob, param, ones);
        case C_SVC_L2:
            std::cout << "\nSVM_SVMType_SVC_L2";
            return new SVM_SVMType_SVC_L2(prob, param, ones);
        case ONE_CLASS:
        case SVDD:
        case SVDD_L2:
            std::cout << "\nSVM_SVMType_ONE_CLASS";
            return new SVM_SVMType_ONE_CLASS(prob, param);
        case EPSILON_SVR:
        case NU_SVR:
            std::cout << "\nSVM_SVMType_SVR";
            return new SVM_SVMType_SVR(prob, param);
        default:
            throw SVM_Exception("Unknown svm_type param");
    }
}

//
// Interface functions
//

/*
SVM_Model svm_train(const SVM_Problem & prob, const SVM_Parameter & param) {
    SVM_Model model(param, 1u);
    model.free_sv = 0; // XXX

    if (param.svm_type == ONE_CLASS ||
        param.svm_type == SVDD ||
        param.svm_type == SVDD_L2 ||
        param.svm_type == EPSILON_SVR ||
        param.svm_type == NU_SVR
    ) {
        // regression or one-class-svm or SVDD
        model.nr_class = 2;
        model.label.clear();
        model.nSV.clear();
        model.probA.clear();
        model.probB.clear();
        //model.sv_coef = Malloc(double *, 1);

        if (param.probability &&
            (param.svm_type == EPSILON_SVR || param.svm_type == NU_SVR)
        ) {
            //model.probA = Malloc(double, 1);
            model.probA.push_back(svm_svr_probability(prob, param));
        }
        std::cout << "\n\nTraining one...\n\n";
        decision_function f = svm_train_one(prob, param, 0, 0);
        std::cout << "\n\nTrained One.\n\n";
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
        model.sv_coef.push_back(coeff);

        f.alpha.clear();
    } else {
        std::cout << "\n\nClassification... (prob.dimension=" << prob.dimension << ")\n\n";
        // classification
        //int l = prob.dimension;
        //int nr_class;

        // index of labels after sorting by label
        std::vector<int> perm(prob.dimension);

        // group training data of the same class
        //std::vector<int> label;
        //std::vector<int> count;
        //std::vector<int> start;
        //svm_group_classes(prob, &nr_class, label, start, count, perm);

        // copy the problem's vectors into a new structure, with the vectors
        // sorted/grouped by label
        node_matrix x(prob.dimension);
        for (unsigned int i = 0; i < prob.dimension; ++i) {
            x.push_back(prob.x[perm[i]]);
        }
        perm.clear();

        std::map<int, SVM_Problem> subproblems = split_problem_by_class(prob);
        int nr_class = subproblems.size();
        std::vector<int> labels(nr_class);
        std::map<int, SVM_Problem>::const_iterator it;
        for (it = subproblems.begin(); it != subproblems.end(); ++it) {
            labels.push_back(it->first);
        }

        // calculate weighted C

        int i;
        std::map<int, double> weighted_C;
        for (i = 0; i < param.nr_weight; ++i) {
            it = subproblems.find(param.weight_label[i]);
            if (it == subproblems.end()) {
                std::cerr << "warning: class label " << param.weight_label[i] << " specified in weight is not found\n";
            } else {
                weighted_C[it->first] = param.C * param.weight[i];
            }
        }

        // train k*(k-1)/2 models
        int n_models = nr_class * (nr_class - 1) / 2;



        // ???
        std::vector<decision_function> f(n_models);

        std::vector<double> probA, probB;
        if (param.probability) {
            probA.assign(n_models, 0.0);
            probB.assign(n_models, 0.0);
        }

        std::map<int, std::vector<bool> > nonzero;
        for (it = subproblems.begin(); it != subproblems.end(); ++it) {
            std::vector<bool> nz(it->second.dimension, false);
            nonzero[it->first] = nz;
        }
        int p = 0, c = 0;
        for (c = 0; c < nr_class; ++c) {
            for (int j = c + 1; j < nr_class; ++j) {
                SVM_Problem sub_prob;
                sub_prob.add(subproblems[c]);
                sub_prob.add(subproblems[j]);

                if (param.probability) {
                    double probabA, probabB;
                    svm_binary_svc_probability(sub_prob, param, weighted_C[c], weighted_C[j], probabA, probabB);
                    probA[p] = probabA;
                    probB[p] = probabB;
                }

                f.push_back(svm_train_one(sub_prob, param, weighted_C[c], weighted_C[j]));
                int k, end;
                for (k = 0, end = subproblems[c].dimension; k < end; ++k) {
                    if (!nonzero[c][k] && fabs(f[p].alpha[k]) > 0) {
                        nonzero[c][k] = true;
                    }
                }
                int offset = subproblems[c].dimension;
                for (k = 0, end = subproblems[j].dimension; k < end; ++k) {
                    if (!nonzero[j][k] && fabs(f[p].alpha[k + offset]) > 0) {
                        nonzero[j][k] = true;
                    }
                }
                sub_prob.x.clear();
                sub_prob.labels.clear();
                ++p;
            }
        }

        std::cout << "\n\nBuilding Output...(nr_class=" << nr_class << ", n_models=" << n_models << ")\n\n";

        // build output
        // http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f402

        model.nr_class = nr_class;

        model.label.reserve(nr_class);
        for (it = subproblems.begin(); it != subproblems.end(); ++it) {
            model.label.push_back(it->first);
        }

        model.rho.reserve(n_models);
        for (i = 0; i < n_models; ++i) {
            model.rho.push_back(f[i].rho[0]);
        }

        if (param.probability) {
            model.probA.reserve(n_models);
            model.probB.reserve(n_models);
            for (i = 0; i < n_models; ++i) {
                model.probA.push_back(probA[i]);
                model.probB.push_back(probB[i]);
            }
        } else {
            model.probA.clear();
            model.probB.clear();
        }

        unsigned int total_sv = 0;
        std::vector<unsigned int> nz_count(nr_class);
        model.nSV.reserve(nr_class);
        for (c = 0; c < nr_class; ++c) {
            unsigned int nSV = 0;
            int end = 0;
            for (i = 0, end = nonzero[c].size(); i < end; ++i) {
                if (nonzero[c][i]) {
                    ++nSV;
                    model.SV.push_back(subproblems[c].x[i]);
                }
            }
            //unsigned int nSV = std::count(nonzero[c].begin(), nonzero[c].end(), true);
            total_sv += nSV;
            model.nSV.push_back(nSV);
            nz_count.push_back(nSV);
        }

        std::cout << "Total nSV = " << total_sv << "\n";

        model.dimension = total_sv;
//        model.SV.reserve(total_sv);
//        p = 0;
//        for (i = 0; i < l; ++i) {
//            if (nonzero[i]) {
//                model.SV.push_back(x[i]);
//            }
//        }

        std::vector<int> nz_start(nr_class);
        nz_start.push_back(0);
        for (c = 1; c < nr_class; ++c) {
            nz_start.push_back(nz_start[c - 1] + nz_count[c - 1]);
        }

        model.sv_coef.reserve(nr_class - 1);
        for (c = 0; c < nr_class - 1; ++c) {
            std::vector<double> sv(total_sv);
            model.sv_coef.push_back(sv);
        }

        p = 0;
        for (c = 0; c < nr_class; c++) {
            for (int j = c + 1; j < nr_class; ++j) {
                // classifier (i,j): coefficients with
                // i are in sv_coef[j-1][nz_start[i]...],
                // j are in sv_coef[i][nz_start[j]...]

                int q = nz_start[c];
                int k, end;
                for (k = 0, end = subproblems[c].dimension; k < end; ++k) {
                    if (nonzero[c][k]) {
                        model.sv_coef[j - 1][q++] = f[p].alpha[k];
                    }
                }
                int offset = subproblems[c].dimension;
                q = nz_start[j];
                for (k = 0, end = subproblems[j].dimension; k < end; ++k) {
                    if (nonzero[j][k]) {
                        model.sv_coef[c][q++] = f[p].alpha[k + offset];
                    }
                }
                ++p;
            }
        }

        //free(label);
        probA.clear();
        probB.clear();
        //free(count);
        //start.clear();
        x.clear();
        weighted_C.clear();
        nonzero.clear();
        f.clear();
        nz_count.clear();
        nz_start.clear();
    }
    return model;
}
*/

// Stratified cross validation

void svm_cross_validation(const SVM_Problem & prob, const SVM_Parameter & param, unsigned int nr_fold, std::vector<double> & target) {
    unsigned int i;
    unsigned int dim = prob.dimension;
    unsigned int nr_class;
    std::vector<int> fold_start(nr_fold + 1);
    std::vector<int> perm(dim);
    for (unsigned int i=0; i<dim; ++i) {
        perm.push_back(i);
    }
    //std::generate(perm.begin(), perm.end(), SequenceNumber);

    SVM_SVMType * svmtype = svmtype_factory(prob, param);
    
    // stratified cv may not give leave-one-out rate
    // Each class to l folds -> some folds may have zero elements
    if ((param.svm_type == C_SVC ||
        param.svm_type == NU_SVC) && nr_fold < dim
    ) {
        std::vector<int> count;
        std::vector<int> label;
        std::vector<int> start;
        svm_group_classes(prob, & nr_class, label, start, count, perm);

        // random shuffle and then data grouped by fold using the array perm
        std::vector<int> index(perm); //copy perm vector into index vector
        unsigned int c;
        for (c = 0; c < nr_class; c++) {
            std::vector<int>::iterator from = index.begin() + start[c];
            std::vector<int>::iterator to   = from + count[c];
            std::random_shuffle(from, to);
        }
        std::vector<int> fold_count(nr_fold);
        for (i = 0; i < nr_fold; ++i) {
            fold_count.push_back(0);
            for (c = 0; c < nr_class; c++) {
                //fold_count[i] += (i + 1) * count[c] / nr_fold - i * count[c] / nr_fold;
                fold_count[i] += count[c] / nr_fold;
            }
        }
        fold_start.push_back(0);
        for (i = 1; i <= nr_fold; ++i) {
            fold_start.push_back(fold_start[i - 1] + fold_count[i - 1]);
        }
        for (c = 0; c < nr_class; c++) {
            int dimc = count[c] / nr_fold;
            for (i = 0; i < nr_fold; ++i) {
                int begin = start[c] + i * dimc;
                int end   = begin + dimc;
                for (int j = begin; j < end; ++j) {
                    perm[fold_start[i]] = index[j];
                    fold_start[i]++;
                }
            }
        }
        fold_start[0] = 0;
        for (i = 1; i <= nr_fold; ++i) {
            fold_start[i] = fold_start[i - 1] + fold_count[i - 1];
        }
        start.clear();
        label.clear();
        count.clear();
        index.clear();
        fold_count.clear();
    } else {
        std::random_shuffle(perm.begin(), perm.end());
        
        for (i = 0; i <= nr_fold; ++i) {
            fold_start.push_back(i * dim / nr_fold);
        }
    }

    target.assign(nr_fold, 0);

    for (i = 0; i < nr_fold; ++i) {
        unsigned int begin = fold_start[i];
        unsigned int end   = fold_start[i + 1];
        unsigned int j;

        SVM_Problem subprob(dim - (end - begin));

        for (j = 0; j < begin; ++j) {
            subprob.x.push_back(prob.x[perm[j]]);
            subprob.labels.push_back(prob.labels[perm[j]]);
        }
        for (j = end; j < dim; ++j) {
            subprob.x.push_back(prob.x[perm[j]]);
            subprob.labels.push_back(prob.labels[perm[j]]);
        }
        SVM_Model submodel = svmtype->svm_train(subprob, param);
        if (param.probability &&
            (param.svm_type == C_SVC || param.svm_type == NU_SVC)) {
            std::vector<double> prob_estimates(submodel.nr_class);
            for (j = begin; j < end; ++j) {
                target[perm[j]] = svmtype->svm_predict_probability(submodel, prob.x[perm[j]], prob_estimates);
            }
            prob_estimates.clear();
        } else {
            for (j = begin; j < end; ++j) {
                target[perm[j]] = svmtype->svm_predict(submodel, prob.x[perm[j]]);
            }
        }
        //delete submodel;
        //delete subprob.x;
        //delete subprob.y;
    }
    fold_start.clear();
    perm.clear();
    delete svmtype;
}


/**
 *
 * @param is
 * @return
 */
int count_elements(std::ifstream & is) {
    std::streamoff original_pos = is.tellg();

    int elements = 0;
    std::string line;
    while (!is.eof()) {
        std::getline(is, line);
        std::string::size_type pos = line.find(':');
        while (pos != std::string::npos) {
            ++elements;
            pos = line.find(':', pos + 1);
        }
    }

    // restore cursor position
    is.seekg(original_pos);
    return elements;
}

const char *svm_check_parameter(const SVM_Problem & prob, const SVM_Parameter & param) {
    // svm_type

    int svm_type = param.svm_type;
    if (svm_type != C_SVC &&
        svm_type != C_SVC_L2 &&
        svm_type != NU_SVC &&
        svm_type != ONE_CLASS &&
        svm_type != EPSILON_SVR &&
        svm_type != NU_SVR &&
        svm_type != SVDD &&
        svm_type != SVDD_L2
    ) {
        return "unknown svm type";
    }

    // kernel_type, degree

    int kernel_type = param.kernel_type;
    if (kernel_type != LINEAR &&
        kernel_type != POLYNOMIAL &&
        kernel_type != GAUSSIAN &&
        kernel_type != SIGMOID &&
        kernel_type != STUMP &&
        kernel_type != PERCEPTRON &&
        kernel_type != LAPLACE &&
        kernel_type != EXPONENTIAL &&
        kernel_type != PRECOMPUTED
    ) {
        return "unknown kernel type";
    }

    if (param.gamma < 0) {
        return "gamma < 0";
    }

    if (param.degree < 0) {
        return "degree of polynomial kernel < 0";
    }

    // cache_size,eps,C,nu,p,shrinking

    if (param.cache_size <= 0) {
        return "cache_size <= 0";
    }

    if (param.eps <= 0) {
        return "eps <= 0";
    }

    if (   svm_type == C_SVC
        || svm_type == EPSILON_SVR
        || svm_type == NU_SVR
    ) {
        if (param.C <= 0) {
            return "C <= 0";
        }
    }

    if (   svm_type == NU_SVC
        || svm_type == ONE_CLASS
        || svm_type == NU_SVR
    ) {
        if (param.nu <= 0 || param.nu > 1) {
            return "nu <= 0 or nu > 1";
        }
    }

    if (svm_type == EPSILON_SVR) {
        if (param.p < 0) {
            return "p < 0";
        }
    }

    if (param.shrinking != 0 &&
        param.shrinking != 1
    ) {
        return "shrinking != 0 and shrinking != 1";
    }

    if (param.probability != 0 &&
        param.probability != 1
    ) {
        return "probability != 0 and probability != 1";
    }

    if (param.probability == 1 &&
        (svm_type == ONE_CLASS || svm_type == SVDD || svm_type == SVDD_L2)
    ) {
        return "one-class SVM probability output not supported yet";
    }


    // check whether nu-svc is feasible
    if (svm_type == NU_SVC) {
        int dim = prob.dimension;
        int nr_class = 0;
        std::vector<int> label;
        std::vector<int> count;

        int i;
        for (i = 0; i < dim; ++i) {
            int this_label = static_cast<int>(prob.labels[i]);
            int c;
            for (c = 0; c < nr_class; ++c) {
                if (this_label == label[c]) {
                    ++count[c];
                    break;
                }
            }
            if (c == nr_class) {
                label.push_back(this_label);
                count.push_back(1);
                ++nr_class;
            }
        }

        for (i = 0; i < nr_class; ++i) {
            int n1 = count[i];
            for (int j = i + 1; j < nr_class; ++j) {
                int n2 = count[j];
                if (param.nu * (n1 + n2) / 2 > min(n1, n2)) {
                    label.clear();
                    count.clear();
                    return "specified nu is infeasible";
                }
            }
        }
        label.clear();
        count.clear();
    }

    return NULL;
}

int svm_check_probability_model(const SVM_Model & model) {
    return (
           (model.param.svm_type == C_SVC || model.param.svm_type == NU_SVC)
        && (model.probA.size() > 0)
        && (model.probB.size() > 0)
    ) || (
           (model.param.svm_type == EPSILON_SVR || model.param.svm_type == NU_SVR)
        && (model.probA.size() > 0)
    );
}
