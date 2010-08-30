/* 
 * File:   Kernel.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "SVC.h"
//#include <typeinfo>

SVM_SVMType_SVC::SVM_SVMType_SVC(const SVM_Problem & prob, const SVM_Parameter & param, std::vector<schar> & y_)
  : SVM_SVMType(prob, param)
{
    //x.assign(prob.x.begin(), prob.x.end());
    //clone(x_, x);
    y.assign(y_.begin(), y_.end());

    cache = new LRUCache<unsigned int, std::vector<Qfloat> >(static_cast<size_t>(param.cache_size*(1<<20)));

    //std::cout << "\n\nSVC prob dim = " << this->kernel->dimension << " (svm_kernel = " << param.kernel_type << ")\n\n";

    QD.reserve(prob.dimension);
    for (unsigned int i=0; i<prob.dimension; i++) {
        //QD.push_back(kernel->kernel(i, i));
        QD.push_back(kernel->kernel(prob.x[i], prob.x[i], param));
    }

    //std::cout << "\n\nEnd SVC constructor\n\n";
}

std::vector<Qfloat> SVM_SVMType_SVC::get_Q(unsigned int column, unsigned int len) const {
    //std::cout << "\nSVC get_Q(" << column << ")";
    std::vector<Qfloat> data;
    //unsigned int start = (cache->find(column, data)) ? data.size() : 0;
    unsigned int start = 0;
    if (start < len) {
        data.resize(len, 0.0);
      //  for (unsigned int j=start; j<len; ++j) {
        for (unsigned int j=0; j<len; ++j) {
            //if (j < 30) {
            //    std::cout << "\n    " << (int)y[column] << " * " << (int)y[j] << " * " << std::setiosflags(std::ios::fixed) << std::setprecision(6) << kernel->kernel(column, j) << " " << typeid(kernel).name();
            //}
            data[j] = static_cast<Qfloat>(static_cast<int>(y[column]) * static_cast<int>(y[j]) * kernel->kernel(column, j));
        }
        //cache->insert(column, data);
    }
    return data;
}

void SVM_SVMType_SVC::swap_index(unsigned int i, unsigned int j) {
    //std::cout << "\nSVC::swap(" << i << ", " << j << ")";
    cache->swap_index(i,j);
    kernel->swap_index(i, j);
    std::swap(y[i], y[j]);
    std::swap(QD[i], QD[j]);
}

std::vector<double> SVM_SVMType_SVC::svm_predict_values(const SVM_Model & model, const std::vector<SVM_Node> & x) const {
    unsigned int i;
    unsigned int nr_class = model.nr_class;
    std::vector<double> dec_values;
    dec_values.reserve(model.dimension);

    std::vector<double> kvalue;
    kvalue.reserve(model.dimension);
    for (i = 0; i < model.dimension; ++i) {
        //kvalue.push_back(Kernel::k_function(x, model.SV[i], model.param));
        //kvalue.push_back(model.kernel->kernel(x, model.SV[i], model.param));
        kvalue.push_back(kernel->kernel(x, model.SV[i], model.param));
    }

    std::vector<int> start;
    start.reserve(nr_class);
    start.push_back(0);
    for (i = 1; i < nr_class; ++i) {
        start.push_back(start[i - 1] + model.nSV[i - 1]);
    }

    int p = 0;
    for (i = 0; i < nr_class; ++i) {
        for (unsigned int j = i + 1; j < nr_class; ++j) {
            double sum = 0;
            int si = start[i];
            int sj = start[j];
            int ci = model.nSV[i];
            int cj = model.nSV[j];

            int k;
            std::vector<double> coef1 = model.sv_coef[j - 1];
            std::vector<double> coef2 = model.sv_coef[i];
            for (k = 0; k < ci; ++k) {
                sum += coef1[si + k] * kvalue[si + k];
            }
            for (k = 0; k < cj; ++k) {
                sum += coef2[sj + k] * kvalue[sj + k];
            }
            sum -= model.rho[p];
            dec_values.push_back(sum);
            ++p;
        }
    }

    kvalue.clear();
    start.clear();

    return dec_values;
}

double SVM_SVMType_SVC::svm_predict(const SVM_Model & model, const std::vector<SVM_Node> & x) const {
    int i;
    int nr_class = model.nr_class;
    std::vector<double> dec_values = svm_predict_values(model, x);
    std::vector<int> vote(nr_class, 0);

    int pos = 0;
    for (i = 0; i < nr_class; ++i) {
        for (int j = i + 1; j < nr_class; ++j) {
            if (dec_values[pos++] > 0) {
                ++vote[i];
            } else {
                ++vote[j];
            }
        }
    }

    int vote_max_idx = 0;
    for (i = 1; i < nr_class; ++i) {
        if (vote[i] > vote[vote_max_idx]) {
            vote_max_idx = i;
        }
    }
    vote.clear();
    dec_values.clear();
    return model.label[vote_max_idx];
}

double SVM_SVMType_SVC::svm_predict_probability(const SVM_Model & model, const std::vector<SVM_Node> & x, std::vector<double> & prob_estimates) const {
    if ((model.probA.size() > 0) && (model.probB.size() > 0)) {
        int i;
        int nr_class = model.nr_class;
        std::vector<double> dec_values = svm_predict_values(model, x);

        double min_prob = 1e-7;
        std::vector<std::vector<double> > pairwise_prob(nr_class);
        for (i = 0; i < nr_class; ++i) {
            std::vector<double> pp(nr_class, 0.0);
            pairwise_prob.push_back(pp);
        }
        int k = 0;
        for (i = 0; i < nr_class; ++i) {
            for (int j = i + 1; j < nr_class; ++j) {
                pairwise_prob[i][j] = min(max(sigmoid_predict(dec_values[k], model.probA[k], model.probB[k]), min_prob), 1 - min_prob);
                pairwise_prob[j][i] = 1 - pairwise_prob[i][j];
                ++k;
            }
        }
        multiclass_probability(nr_class, pairwise_prob, prob_estimates);

        int prob_max_idx = 0;
        for (i = 1; i < nr_class; ++i) {
            if (prob_estimates[i] > prob_estimates[prob_max_idx]) {
                prob_max_idx = i;
            }
        }
        pairwise_prob.clear();
        dec_values.clear();
        return model.label[prob_max_idx];
    }
    return svm_predict(model, x);
}

// Method 2 from the multiclass_prob paper by Wu, Lin, and Weng
/**
 *
 * @param k
 * @param r
 * @param p prob_estimates
 */
void SVM_SVMType_SVC::multiclass_probability(int k, std::vector<std::vector<double> > & pairwise_prob, std::vector<double> & prob_estimates) const {
    int t, j;
    int iter = 0, max_iter = max(100, k);
    double pQp, eps = 0.005 / k;

    std::vector<std::vector<double> > Q(k);
    for (int i = 0; i < k; ++i) {
        std::vector<double> pp(k, 0.0);
        Q.push_back(pp);
    }

    prob_estimates.assign(k, 1.0 / k); // Valid if k = 1

    for (t = 0; t < k; t++) {
        Q[t][t] = 0;
        for (j = 0; j < t; j++) {
            Q[t][t] += pairwise_prob[j][t] * pairwise_prob[j][t];
            Q[t][j] = Q[j][t];
        }
        for (j = t + 1; j < k; j++) {
            Q[t][t] += pairwise_prob[j][t] * pairwise_prob[j][t];
            Q[t][j] = -pairwise_prob[j][t] * pairwise_prob[t][j];
        }
    }
    for (iter = 0; iter < max_iter; iter++) {
        std::vector<double> Qp(k, 0.0);
        // stopping condition, recalculate QP,pQP for numerical accuracy
        pQp = 0;
        for (t = 0; t < k; t++) {
            Qp[t] = 0;
            for (j = 0; j < k; j++) {
                Qp[t] += Q[t][j] * prob_estimates[j];
            }
            pQp += prob_estimates[t] * Qp[t];
        }
        double max_error = 0;
        for (t = 0; t < k; t++) {
            double error = fabs(Qp[t] - pQp);
            if (error > max_error) {
                max_error = error;
            }
        }
        if (max_error < eps) {
            break;
        }

        for (t = 0; t < k; t++) {
            double diff = (-Qp[t] + pQp) / Q[t][t];
            prob_estimates[t] += diff;
            pQp = (pQp + diff * (diff * Q[t][t] + 2 * Qp[t])) / (1 + diff) / (1 + diff);
            for (j = 0; j < k; j++) {
                Qp[j] = (Qp[j] + diff * Q[t][j]) / (1 + diff);
                prob_estimates[j] /= (1 + diff);
            }
        }
    }
    if (iter >= max_iter) {
        std::cout << "Exceeds max_iter in multiclass_prob\n";
    }
    Q.clear();
}

double SVM_SVMType_SVC::sigmoid_predict(double decision_value, double A, double B) const {
    double fApB = decision_value * A + B;
    if (fApB >= 0) {
        return exp(-fApB) / (1.0 + exp(-fApB));
    }
    return 1.0 / (1 + exp(fApB));
}

/**
 * Cross-validation decision values for probability estimates
 *
 * @param prob
 * @param param
 * @param Cp
 * @param Cn
 * @param probA
 * @param probB
 */ 
void SVM_SVMType_SVC::svm_binary_svc_probability(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    double Cp,
    double Cn,
    double& probA,
    double& probB
) const {
    unsigned int i;
    unsigned int nr_fold = 5;
    std::vector<int> perm(prob.dimension);
    std::vector<double> dec_values; //decision values
    dec_values.assign(prob.dimension, 0);

    // random shuffle
    //std::generate(perm.begin(), perm.end(), SequenceNumber);
    for (unsigned int i=0; i<prob.dimension; ++i) {
        perm.push_back(i);
    }
    std::random_shuffle(perm.begin(), perm.end());

    for (i = 0; i < nr_fold; i++) {
        unsigned int dim = prob.dimension / nr_fold;
        unsigned int begin = i * dim;
        unsigned int end = begin * dim;
        unsigned int j;
        
        SVM_Problem subprob(prob.dimension - (end - begin));

        for (j = 0; j < begin; j++) {
            subprob.x.push_back(prob.x[perm[j]]);
            subprob.labels.push_back(prob.labels[perm[j]]);
        }
        for (j = end; j < prob.dimension; j++) {
            subprob.x.push_back(prob.x[perm[j]]);
            subprob.labels.push_back(prob.labels[perm[j]]);
        }
        unsigned int p_count = 0, n_count = 0, size = 0;
        for (j = 0, size = subprob.x.size(); j < size; j++) {
            if (subprob.labels[j] > 0) {
                ++p_count;
            } else {
                ++n_count;
            }
        }

        if (p_count == 0 && n_count == 0) {
            for (j = begin; j < end; j++) {
                dec_values[perm[j]] = 0;
            }
        } else if (p_count > 0 && n_count == 0) {
            for (j = begin; j < end; j++) {
                dec_values[perm[j]] = 1;
            }
        } else if (p_count == 0 && n_count > 0) {
            for (j = begin; j < end; j++) {
                dec_values[perm[j]] = -1;
            }
        } else {
            SVM_Parameter subparam = param;
            subparam.probability = 0;
            subparam.C = 1.0;
            subparam.nr_weight = 2;
            subparam.weight_label.push_back(+1);
            subparam.weight_label.push_back(-1);
            subparam.weight.push_back(Cp);
            subparam.weight.push_back(Cn);
            SVM_Model submodel = svm_train(subprob, subparam);
            for (j = begin; j < end; j++) {
                std::vector<double> dec_values_sub = svm_predict_values(submodel, prob.x[perm[j]]);
                // copy the decision values for this subproblem to dec_values starting from dec_values[perm[j]]
                std::copy(dec_values_sub.begin(), dec_values_sub.end(), dec_values.begin() + perm[j]);
                // ensure +1 -1 order; reason not using CV subroutine
                dec_values[perm[j]] *= submodel.label[0];
            }
            //delete submodel;
            //delete subparam;
        }
        subprob.x.clear();
        subprob.labels.clear();
    }
    sigmoid_train(prob.dimension, dec_values, prob.labels, probA, probB);
    dec_values.clear();
    perm.clear();
}

/**
 * Platt's binary SVM Probablistic Output: an improvement from Lin et al.
 *
 * @param l
 * @param dec_values
 * @param labels
 * @param A
 * @param B
 */
void SVM_SVMType_SVC::sigmoid_train(
    int l,
    const std::vector<double> & dec_values,
    const std::vector<double> & labels,
    double& A,
    double& B
) const {
    double prior1 = 0, prior0 = 0;
    int i;

    for (i = 0; i < l; i++) {
        if (labels[i] > 0) {
            prior1 += 1;
        } else {
            prior0 += 1;
        }
    }

    int max_iter = 100; // Maximal number of iterations
    double min_step = 1e-10; // Minimal step taken in line search
    double sigma = 1e-12; // For numerically strict PD of Hessian
    double eps = 1e-5;
    double hiTarget = (prior1 + 1.0) / (prior1 + 2.0);
    double loTarget = 1 / (prior0 + 2.0);
    std::vector<double> t(l);
    double fApB, p, q, h11, h22, h21, g1, g2, det, dA, dB, gd, stepsize;
    double newA, newB, newf, d1, d2;
    int iter;

    // Initial Point and Initial Fun Value
    A = 0.0;
    B = log((prior0 + 1.0) / (prior1 + 1.0));
    double fval = 0.0;

    for (i = 0; i < l; i++) {
        if (labels[i] > 0) {
            t.push_back(hiTarget);
        } else {
            t.push_back(loTarget);
        }
        fApB = dec_values[i] * A + B;
        if (fApB >= 0) {
            fval += t[i] * fApB + log(1 + exp(-fApB));
        } else {
            fval += (t[i] - 1) * fApB + log(1 + exp(fApB));
        }
    }
    for (iter = 0; iter < max_iter; iter++) {
        // Update Gradient and Hessian (use H' = H + sigma I)
        h11 = sigma; // numerically ensures strict PD
        h22 = sigma;
        h21 = 0.0;
        g1 = 0.0;
        g2 = 0.0;
        for (i = 0; i < l; i++) {
            fApB = dec_values[i] * A + B;
            if (fApB >= 0) {
                p = exp(-fApB) / (1.0 + exp(-fApB));
                q = 1.0 / (1.0 + exp(-fApB));
            } else {
                p = 1.0 / (1.0 + exp(fApB));
                q = exp(fApB) / (1.0 + exp(fApB));
            }
            d2 = p*q;
            h11 += dec_values[i] * dec_values[i] * d2;
            h22 += d2;
            h21 += dec_values[i] * d2;
            d1 = t[i] - p;
            g1 += dec_values[i] * d1;
            g2 += d1;
        }

        // Stopping Criteria
        if (fabs(g1) < eps && fabs(g2) < eps) {
            break;
        }

        // Finding Newton direction: -inv(H') * g
        det = h11 * h22 - h21*h21;
        dA = -(h22 * g1 - h21 * g2) / det;
        dB = -(-h21 * g1 + h11 * g2) / det;
        gd = g1 * dA + g2*dB;


        stepsize = 1; // Line Search
        while (stepsize >= min_step) {
            newA = A + stepsize * dA;
            newB = B + stepsize * dB;

            // New function value
            newf = 0.0;
            for (i = 0; i < l; i++) {
                fApB = dec_values[i] * newA + newB;
                if (fApB >= 0) {
                    newf += t[i] * fApB + log(1 + exp(-fApB));
                } else {
                    newf += (t[i] - 1) * fApB + log(1 + exp(fApB));
                }
            }
            // Check sufficient decrease
            if (newf < fval + 0.0001 * stepsize * gd) {
                A = newA;
                B = newB;
                fval = newf;
                break;
            } else {
                stepsize = stepsize / 2.0;
            }
        }

        if (stepsize < min_step) {
            std::cout << "Line search fails in two-class probability estimates\n";
            break;
        }
    }

    if (iter >= max_iter) {
        std::cout << "Reaching maximal iterations in two-class probability estimates\n";
    }
    t.clear();
}


SVM_Model SVM_SVMType_SVC::svm_train(const SVM_Problem & prob, const SVM_Parameter & param) const {
    SVM_Model model(param, 1u);
    model.free_sv = 0; // XXX

    //std::cout << "\n\nClassification...\n\n";

    //std::cout << "PROB[0] SIZE =" << prob.x[0].size() << " (dim=" << prob.dimension << ", prob.x.size=" << prob.x.size() << ")\n\n";
    // classification
    //int l = prob.dimension;
    //int nr_class;

    // index of labels after sorting by label
    std::vector<int> perm(prob.dimension);

    // group training data of the same class
    unsigned int nr_class;
    std::vector<int> label;
    std::vector<int> count;
    std::vector<int> start;
    svm_group_classes(prob, &nr_class, label, start, count, perm);

    std::cout << "\n\nGrouped classes: " << nr_class << "\n\n";

    // group training data of the same class
    node_matrix x;
    unsigned int i;

    x.reserve(prob.dimension);
    for (i = 0; i < prob.dimension; ++i) {
        x.push_back(prob.x[perm[i]]);
    }

    //std::cout << "x[" << x.size() << "] SIZE =" << x[0].size() << ", cnt=" << cnt << ", prob.dimension=" << prob.dimension << ", perm.size=" << perm.size() << "\n\n";

    // calculate weighted C

    std::vector<double> weighted_C(nr_class, param.C);
    for (i = 0; i < param.nr_weight; ++i) {
        unsigned int j;
        for (j = 0; j < nr_class; ++j) {
            if (param.weight_label[i] == label[j]) {
                break;
            }
        }
        if (j == nr_class) {
            std::cerr << "warning: class label " << param.weight_label[i] << " specified in weight is not found\n";
        } else {
            weighted_C[j] *= param.weight[i];
        }
    }

//        std::cout << "\nweigthed_C[] = \n";
//        for (i=0; i<nr_class; ++i) {
//            std::cout << weighted_C[i] << " ";
//        }
//        std::cout << "\nparam.weight[] = \n";
//        for (i=0; i<param.nr_weight; ++i) {
//            std::cout << param.weight[i] << " ";
//        }

    // train k*(k-1)/2 models

    std::vector<bool> nonzero(prob.dimension, false);
    //std::vector<decision_function> f(nr_class * (nr_class - 1) / 2);
    std::vector<decision_function> f;

    std::vector<double> probA, probB;
    if (param.probability) {
        probA.reserve(nr_class * (nr_class - 1) / 2);
        probB.reserve(nr_class * (nr_class - 1) / 2);
    }

    int p = 0;
    for (i = 0; i < nr_class; ++i) {
        for (unsigned int j = i + 1; j < nr_class; ++j) {
            std::cout << "\n[sub_prob " << i << " " << j << "]";
            int si = start[i], sj = start[j];
            int ci = count[i], cj = count[j];

            SVM_Problem sub_prob(ci + cj);
            //std::cout << "(si, sj, ci, cj, dim) = (" << si << ", " << sj << ", " << ci << ", " << cj << ", " << sub_prob.dimension <<")\n";

            int k;
            for (k = 0; k < ci; ++k) {
                sub_prob.x.push_back(x[si + k]);
                sub_prob.labels.push_back(+1);
            }
            for (k = 0; k < cj; ++k) {
                sub_prob.x.push_back(x[sj + k]);
                sub_prob.labels.push_back(-1);
            }

            //std::cout << "SUBPROB[0] SIZE =" << sub_prob.x[0].size() << "\n\n";
            //exit(0);

            if (param.probability) {
                //std::cout << "\n\n param.probability = " << param.probability << "\n\n";
                svm_binary_svc_probability(sub_prob, param, weighted_C[i], weighted_C[j], probA[p], probB[p]);
            }

            //std::cout << "\n\nsvm_train_one()\n\n";
            decision_function df = svm_train_one(sub_prob, param, weighted_C[i], weighted_C[j]);
            
            f.push_back(df);

//            std::cout << "\n(prob.dimension, ci, si, cj, sj) = (" << prob.dimension << ", " << ci << ", " << si << ", " << cj << ", " << sj << ")";
//            std::cout << "\nnonzero[" << nonzero.size() << ", " << "df.alpha[" << df.alpha.size() << "], accessing si+k=" << (si + ci -1) << ", sj+k=" << (sj + cj -1);
//
//            std::cout.flush();
            //exit(0);

            for (k = 0; k < ci; ++k) {
                if (!nonzero[si + k] && fabs(df.alpha[k]) > 0) {
                    nonzero[si + k] = true;
                }
            }
            for (k = 0; k < cj; ++k) {
                if (!nonzero[sj + k] && fabs(df.alpha[ci + k]) > 0) {
                    nonzero[sj + k] = true;
                }
            }
//            std::cout.flush();
//            std::cout << "\n\nSo far so good\n\n";
//            std::cout.flush();
//            exit(0);
            sub_prob.x.clear();
            sub_prob.labels.clear();
            ++p;
        }
    }

//    std::cout << "\n\nBuild output\n\n";
//    std::cout.flush();

//    for (std::vector<decision_function>::const_iterator it = f.begin(); it != f.end(); ++it) {
//        std::cout << "\n - Rho size: " << it->rho.size();
//        std::cout.flush();
//    }
//
//    std::cout << "\nLABELS: ";
//    for (std::vector<int>::const_iterator it = label.begin(); it != label.end(); ++it) {
//        std::cout << *it;
//    }


    model.nr_class = nr_class;
    model.label.assign(label.begin(), label.end());

//    std::cout << "\nLabels assigned";
//    std::cout.flush();
//    std::cout << "\nf.size = " << f.size();
//    std::cout.flush();

    model.rho.assign(nr_class * (nr_class - 1) / 2, 0.0);
    for (i = 0; i < nr_class * (nr_class - 1) / 2; ++i) {
        model.rho[i] = f[i].rho[0];
    }
//    std::cout << "\nRho assigned";
//    std::cout.flush();

    if (param.probability) {
        model.probA.assign(probA.begin(), probA.end());
        model.probB.assign(probB.begin(), probB.end());
//        model.probA.reserve(nr_class * (nr_class - 1) / 2);
//        model.probB.reserve(nr_class * (nr_class - 1) / 2);
//        for (i = 0; i < nr_class * (nr_class - 1) / 2; ++i) {
//            model.probA.push_back(probA[i]);
//            model.probB.push_back(probB[i]);
//        }
    }


    // build output

    unsigned int total_sv = 0;
    std::vector<int> nz_count;
    nz_count.reserve(nr_class);
    model.nSV.reserve(nr_class);
    for (i = 0; i < nr_class; ++i) {
        int nSV = 0;
        for (int j = 0; j < count[i]; ++j) {
            if (nonzero[start[i] + j]) {
                ++nSV;
                ++total_sv;
            }
        }
        model.nSV.push_back(nSV);
        nz_count.push_back(nSV);
    }

    std::cout << "\nTotal nSV = " << total_sv;

    model.dimension = total_sv;
    model.SV.reserve(total_sv);
    //std::cout << "\nnonzero[] = ";
    for (i = 0; i < prob.dimension; ++i) {
        //std::cout << " " << static_cast<int>(nonzero[i]);
        if (nonzero[i]) {
            model.SV.push_back(x[i]);
        }
    }

    std::vector<int> nz_start;
    nz_start.reserve(nr_class);
    nz_start.push_back(0);
    for (i = 1; i < nr_class; ++i) {
        nz_start.push_back(nz_start[i - 1] + nz_count[i - 1]);
    }
    //std::cout << "\nnz_start=[" << nz_start[0] << "," << nz_start[1] << "], nz_count=[" << nz_count[0] << "," << nz_count[1] << "]";

    model.sv_coef.clear();
    model.sv_coef.reserve(nr_class - 1);
    for (i = 0; i < nr_class - 1; ++i) {
        std::vector<double> vd(total_sv, 0.0);
        model.sv_coef.push_back(vd);
    }

    p = 0;
    for (i = 0; i < nr_class; ++i) {
        for (unsigned int j = i + 1; j < nr_class; ++j) {
            // classifier (i,j): coefficients with
            // i are in sv_coef[j-1][nz_start[i]...],
            // j are in sv_coef[i][nz_start[j]...]

            int si = start[i];
            int sj = start[j];
            int ci = count[i];
            int cj = count[j];

            //std::cout << "\n(si, sj, ci, cj) = (" << si << ", " << sj << ", " << ci << ", " << cj << ")";

            int q = nz_start[i];
            int k;
            //std::cout << "\n--[" << i << "][" << j << "][" << ci << "]-----------------";
            for (k = 0; k < ci; ++k) {
                //std::cout << "\n";
                if (nonzero[si + k]) {
                    model.sv_coef[j - 1][q++] = f[p].alpha[k];
                    //std::cout << f[p].alpha[k] << " ";
                }
            }
            q = nz_start[j];
            //std::cout << "\n--[" << i << "][" << j << "][" << cj << "]-----------------";
            for (k = 0; k < cj; ++k) {
                //std::cout << "\n";
                if (nonzero[sj + k]) {
                    model.sv_coef[i][q++] = f[p].alpha[ci + k];
                    //std::cout << f[p].alpha[ci + k] << " ";
                }
            }
            ++p;
        }
    }

//    std::cout << "\n-------------------------------------";
//    for (i = 0; i < total_sv; ++i) {
//        std::cout << model.sv_coef[0][i] << " ";
//        if (i % 20 == 19) std::cout << "\n";
//    }

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

    return model;
}

void SVM_SVMType_SVC::solve(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    std::vector<double> & alpha,
    SVM_SolutionInfo & si,
    double Cp,
    double Cn
) const {
    unsigned int l = prob.dimension;
    std::vector<double> minus_ones(l, -1);
    std::vector<schar> labels(l, 1);
    unsigned int i;

    alpha.assign(l, 0);
    for (i = 0; i < l; i++) {
        if (prob.labels[i] > 0) {
            labels[i] = +1;
        } else {
            labels[i] = -1;
        }
    }

    //std::cout << "\n\nSolving SVC(prob,param,labels) with DIM=" << prob.dimension << "(" << prob.x.size() << ")\n\n";

    SVM_Solver s;
    SVM_SVMType * svmtype = new SVM_SVMType_SVC(prob, param, labels);
    s.solve(l, svmtype, minus_ones, labels,
        alpha, Cp, Cn, param.eps, si, param.shrinking);
    delete svmtype;

    double sum_alpha = 0;
    for (i = 0; i < l; i++) {
        sum_alpha += alpha[i];
    }

    if (Cp == Cn) {
        std::cout << "nu = " << (sum_alpha / (Cp * prob.dimension)) << "\n";
    }
    for (i = 0; i < l; i++) {
        alpha[i] *= labels[i];
    }

    minus_ones.clear();
    labels.clear();
}

SVM_SVMType_SVC::~SVM_SVMType_SVC() {
    y.clear();
    delete cache;
    QD.clear();
}
