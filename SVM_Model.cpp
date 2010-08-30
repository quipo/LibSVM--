/* 
 * File:   SVM_Model.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "SVM_Model.h"

static const char *svm_type_table[] = {
    "c_svc", "c_svc_l2", "nu_svc", "one_class", "epsilon_svr", "nu_svr", "svdd", "svdd_l2", NULL
};

static const char *kernel_type_table[] = {
    "linear", "polynomial", "gaussian", "sigmoid", "stump", "perc", "laplace", "exponential", "precomputed", "neural", NULL
};

/**
 *
 * @param SVM_Parameter param
 * @param unsigned int  dim
 */
SVM_Model::SVM_Model(const SVM_Parameter & param, unsigned int dim) : param(param), dimension(dim) {
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
        case NEURAL:
            kernel = new SVM_Kernel_Neural(param);
            break;
    }
}

SVM_Model::~SVM_Model() {
    SV.clear();
    sv_coef.clear();
    rho.clear();
    label.clear();
    probA.clear();
    probB.clear();
    nSV.clear();
    //delete kernel;
}

/**
 * Load model from file
 *
 * @param model_file_name
 */
SVM_Model::SVM_Model(const char *model_file_name) {
    std::ifstream is(model_file_name, std::ios::binary);
    if (!is.is_open()) {
        return;
    }

    std::string line, cmd;
    while (1) {
        std::getline(is, line);
        std::istringstream ln(line);
        ln >> cmd;
        if (cmd == "svm_type") {
            ln >> cmd;
            int i;
            for (i = 0; svm_type_table[i]; ++i) {
                if (cmd == svm_type_table[i]) {
                    param.svm_type = i;
                    break;
                }
            }
            if (svm_type_table[i] == NULL) {
                std::cerr << "unknown svm type.\n";
                return;
            }
        } else if (cmd == "kernel_type") {
            ln >> cmd;
            int i;
            for (i = 0; kernel_type_table[i]; ++i) {
                if (cmd == kernel_type_table[i]) {
                    param.kernel_type = i;
                    break;
                }
            }
            if (kernel_type_table[i] == NULL) {
                std::cerr << "unknown kernel function.\n";
                return;
            }
        } else if (cmd == "degree") {
            ln >> param.degree;
        } else if (cmd == "gamma") {
            ln >> param.gamma;
        } else if (cmd == "coef0") {
            ln >> param.coef0;
        } else if (cmd == "nr_class") {
            ln >> nr_class;
        } else if (cmd == "total_sv") {
            ln >> dimension;
        } else if (cmd == "rho") {
            unsigned int n = nr_class * (nr_class - 1) / 2;
            rho.reserve(n);
            double rhoEl;
            for (unsigned int i = 0; i < n; ++i) {
                ln >> rhoEl;
                rho.push_back(rhoEl);
            }
        } else if (cmd == "label") {
            unsigned int n = nr_class;
            label.reserve(n);
            int c;
            for (unsigned int i = 0; i < n; ++i) {
                ln >> c;
                label.push_back(c);
            }
        } else if (cmd == "probA") {
            unsigned int n = nr_class * (nr_class - 1) / 2;
            probA.reserve(n);
            double prob;
            for (unsigned int i = 0; i < n; ++i) {
                ln >> prob;
                probA.push_back(prob);
            }
        } else if (cmd == "probB") {
            int n = nr_class * (nr_class - 1) / 2;
            probB.reserve(n);
            double prob;
            for (int i = 0; i < n; ++i) {
                ln >> prob;
                probB.push_back(prob);
            }
        } else if (cmd == "nr_sv") {
            unsigned int n = nr_class;
            nSV.reserve(n);
            unsigned int u;
            for (unsigned int i = 0; i < n; ++i) {
                ln >> u;
                nSV.push_back(u);
            }
        } else if (cmd == "SV") {
            break;
        } else {
            std::cerr << "unknown text in model file: [" << cmd << "]\n";
            return;
        }
    }

    // read sv_coef and SV

    //int elements = count_elements(fp) + dimension;

    unsigned int m = nr_class - 1;
    sv_coef.clear();
    sv_coef.reserve(m);
    unsigned int i;
    for (i = 0; i < m; ++i) {
        std::vector<double> vd(dimension, 0.0);
        sv_coef.push_back(vd);
    }

    for (i = 0; i < dimension; ++i) {
        std::getline(is, line);
        std::vector<SVM_Node> x_space;

        std::istringstream p(line);
        double d;
        for (unsigned int k = 0; k < m; ++k) {
            p >> d;
            sv_coef[k][i] = d;
        }

        std::string vpair;
        while (p >> vpair) {
            std::string index;

            SVM_Node node;
            std::stringstream ss(vpair);
            //idx = strtok(NULL, ":");
            //val = strtok(NULL, " \t");
            std::getline(ss, index, ':');
            std::istringstream idx(index);
            idx >> node.index;
            ss >> node.value;
            x_space.push_back(node);
        }
        SV.push_back(x_space);
    }

    is.close();

    free_sv = 1; // XXX
}

/**
 *
 * @param model_file_name
 * @return int
 */
int SVM_Model::save(const char *model_file_name) const {
    std::ofstream fp(model_file_name);
    if (!fp.is_open()) {
        return -1;
    }

    fp << "svm_type "    << svm_type_table[param.svm_type]       << "\n";
    fp << "kernel_type " << kernel_type_table[param.kernel_type] << "\n";

    if (param.kernel_type == POLYNOMIAL) {
        fp << "degree " << param.degree << "\n";
    }

    if (   param.kernel_type == POLYNOMIAL
        || param.kernel_type == GAUSSIAN
        || param.kernel_type == SIGMOID
        || param.kernel_type == LAPLACE
        || param.kernel_type == EXPONENTIAL
        || param.kernel_type == NEURAL
    ) {
        fp << "gamma " << param.gamma << "\n";
    }

    if (   param.kernel_type == POLYNOMIAL
        || param.kernel_type == SIGMOID
        || param.kernel_type == STUMP
        || param.kernel_type == PERCEPTRON
        || param.kernel_type == NEURAL
    ) {
        fp << "coef0 " << param.coef0 << "\n";
    }

    fp << "nr_class " << nr_class << "\n";
    fp << "total_sv " << dimension << "\n";

    {
        fp << "rho";
        for (unsigned int i = 0; i < nr_class * (nr_class - 1) / 2; ++i) {
            fp << " " << rho[i];
        }
        fp << "\n";
    }

    if (label.size() > 0) {
        fp << "label";
        for (std::vector<int>::const_iterator it = label.begin(); it < label.end(); ++it) {
            fp << " " << static_cast<int>(*it);
        }
        fp << "\n";
    }

    unsigned int i = 0;
    if (probA.size() > 0) {
        // regression has probA only
        fp << "probA";
        for (i = 0; i < nr_class * (nr_class - 1) / 2; ++i) {
            fp << " " << probA[i];
        }
        fp << "\n";
    }
    if (probB.size() > 0) {
        fp << "probB";
        for (i = 0; i < nr_class * (nr_class - 1) / 2; ++i) {
            fp << " " << probB[i];
        }
        fp << "\n";
    }

    if (nSV.size() > 0) {
        fp << "nr_sv";
        for (i = 0; i < nr_class; ++i) {
            fp << " " << nSV[i];
        }
        fp << "\n";
    }

    fp << "SV\n";

    for (i = 0; i < dimension; ++i) {
        for (unsigned int j = 0; j < nr_class - 1; ++j) {
            //fp << std::setiosflags(std::ios::fixed) << std::setprecision(16) << sv_coef[j][i] << " ";
            fp << std::setprecision(16) << sv_coef[j][i] << " ";
        }

        //const SVM_Node p = SV[i];

        if (param.kernel_type == PRECOMPUTED) {
            fp << "0:" << static_cast<int>(SV[i][0].value) << " ";
        } else {
            for (std::vector<SVM_Node>::const_iterator p = SV[i].begin(); p < SV[i].end(); ++p) {
                fp << *p;
            }
        }
        fp << "\n";
    }
    if (!fp.good()) {
        return -1;
    }
    fp.close();
    return 0;
}
