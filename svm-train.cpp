#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "static_trainers.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void print_null(const char *s);
void exit_with_help();
void exit_input_error(int line_num);
SVM_Parameter parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name);
void do_cross_validation(const SVM_Problem & prob, const SVM_Parameter & param, unsigned int nr_fold);
int main(int argc, char **argv);

void print_null(const char *s) {
    if (s) {
        
    }
}

void exit_with_help() {
    printf(
        "Usage: svm-train [options] training_set_file [model_file]\n"
        "options:\n"
        "-s svm_type : set type of SVM (default 0)\n"
        "\t%d -- C-SVC (L1SVM)\n"
        "\t%d -- C-SVC (L2SVM)\n"
        "\t%d -- nu-SVC\n"
        "\t%d -- one-class SVM\n"
        "\t%d -- epsilon-SVR\n"
        "\t%d -- nu-SVR\n"
        "\t%d -- SVDD (L1SVM)\n"
        "\t%d -- SVDD (L2SVM)\n"
        "-t kernel_type : set type of kernel function (default 2)\n"
        "\t%d -- linear: u'*v\n"
        "\t%d -- polynomial: (gamma*u'*v + coef0)^degree\n"
        "\t%d -- gaussian: exp(-gamma*||u-v||^2)\n"
        "\t%d -- sigmoid: tanh(gamma*u'*v + coef0)\n"
        "\t%d -- stump: -|u-v| + coef0\n"
        "\t%d -- perceptron: -||u-v|| + coef0\n"
        "\t%d -- laplacian: exp(-gamma*|u-v|)\n"
        "\t%d -- exponential: exp(-gamma*||u-v||)\n"
        "\t%d -- precomputed kernel (kernel values in training_set_file)\n"
        "-d degree: set degree in kernel function (default 3)\n"
        "-g gamma: set gamma in kernel function (default 1/num_features)\n"
        "-r coef0: set coef0 in kernel function (default 0)\n"
        "-c cost: set the parameter C of C-SVC (L1/L2-SVM), SVDD (L1/L2-SVM), epsilon-SVR, and nu-SVR (default 1)\n"
        "-n nu: set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
        "-p epsilon: set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
        "-m cachesize: set cache memory size in MB (default 100)\n"
        "-e epsilon: set tolerance of termination criterion (default 0.001)\n"
        "-h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
        "-b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
        "-wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
        "-v n: n-fold cross validation mode\n"
        "-q: quiet mode (no outputs)\n", C_SVC, C_SVC_L2, NU_SVC, ONE_CLASS,
        EPSILON_SVR, NU_SVR, SVDD, SVDD_L2,
        LINEAR, POLYNOMIAL, GAUSSIAN, SIGMOID, STUMP,
        PERCEPTRON, LAPLACE, EXPONENTIAL, PRECOMPUTED);
    exit(1);
}

void exit_input_error(int line_num) {
    std::cerr << "Wrong input format at line " << line_num << "\n";
    exit(1);
}

int cross_validation;
int nr_fold;

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
    char input_file_name[1024];
    char model_file_name[1024];
    const char *error_msg;

    SVM_Parameter param = parse_command_line(argc, argv, input_file_name, model_file_name);
    SVM_Problem prob(input_file_name, param);
    
    error_msg = svm_check_parameter(prob, param);
    
    if (error_msg) {
        std::cerr << "Error: " << error_msg << "\n";
        exit(1);
    }

    std::cout << param;

    if (cross_validation) {
        std::cout << "\n\nDo cross validation...\n\n";
        do_cross_validation(prob, param, nr_fold);
    } else {
        std::cout << "\n\nRead Problem. Training...\n\n";
        SVM_SVMType * svmtype = svmtype_factory(prob, param);
        SVM_Model model = svmtype->svm_train(prob, param);
        std::cout << "\n\ntrained. Saving model...\n\n";
        model.save(model_file_name);
        std::cout << "\n\nSAVED\n\n";
        //delete model;
    }
    //delete param;
    //delete prob;
    //free(line);

    return 0;
}

void do_cross_validation(const SVM_Problem & prob, const SVM_Parameter & param, unsigned int nr_fold) {
    unsigned int i;
    int total_correct = 0;
    double total_error = 0;
    double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
    std::vector<double> target(prob.dimension);

    svm_cross_validation(prob, param, nr_fold, target);
    if (param.svm_type == EPSILON_SVR ||
        param.svm_type == NU_SVR) {
        for (i = 0; i < prob.dimension; i++) {
            double y = prob.labels[i];
            double v = target[i];
            total_error += (v - y)*(v - y);
            sumv += v;
            sumy += y;
            sumvv += v*v;
            sumyy += y*y;
            sumvy += v*y;
        }
        std::cout << "Cross Validation Mean squared error = " << std::setprecision(8) << (total_error / prob.dimension) << "\n";
        std::cout << "Cross Validation Squared correlation coefficient = " << std::setprecision(8) <<
            (((prob.dimension * sumvy - sumv * sumy)*(prob.dimension * sumvy - sumv * sumy)) /
            ((prob.dimension * sumvv - sumv * sumv)*(prob.dimension * sumyy - sumy * sumy)))
            << "\n";
    } else {
        for (i = 0; i < prob.dimension; i++) {
            if (target[i] == prob.labels[i]) {
                ++total_correct;
            }
        }
        std::cout << "Cross Validation Accuracy = " << std::setprecision(8) << (100.0 * total_correct / prob.dimension) << "%\n";
    }
    target.clear();
}

SVM_Parameter parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name) {
    int i;

    SVM_Parameter param;
    // default values
    param.svm_type = C_SVC;
    param.kernel_type = GAUSSIAN;
    param.degree = 3;
    param.gamma = 0; // 1/num_features
    param.coef0 = 0;
    param.nu = 0.5;
    param.cache_size = 100;
    param.C = 1;
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label.clear();
    param.weight.clear();
    cross_validation = 0;

    // parse options
    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-') break;
        if (++i >= argc)
            exit_with_help();
        switch (argv[i - 1][1]) {
            case 's':
                param.svm_type = atoi(argv[i]);
                break;
            case 't':
                param.kernel_type = atoi(argv[i]);
                break;
            case 'd':
                param.degree = atoi(argv[i]);
                break;
            case 'g':
                param.gamma = atof(argv[i]);
                break;
            case 'r':
                param.coef0 = atof(argv[i]);
                break;
            case 'n':
                param.nu = atof(argv[i]);
                break;
            case 'm':
                param.cache_size = atof(argv[i]);
                break;
            case 'c':
                param.C = atof(argv[i]);
                break;
            case 'e':
                param.eps = atof(argv[i]);
                break;
            case 'p':
                param.p = atof(argv[i]);
                break;
            case 'h':
                param.shrinking = atoi(argv[i]);
                break;
            case 'b':
                param.probability = atoi(argv[i]);
                break;
            case 'q':
                //svm_print_string = &print_null;
                i--;
                break;
            case 'v':
                cross_validation = 1;
                nr_fold = atoi(argv[i]);
                if (nr_fold < 2) {
                    std::cerr << "n-fold cross validation: n must >= 2\n";
                    exit_with_help();
                }
                break;
            case 'w':
                ++param.nr_weight;
                param.weight_label.reserve(param.nr_weight);
                param.weight.reserve(param.nr_weight);
                param.weight_label[param.nr_weight - 1] = atoi(&argv[i - 1][2]);
                param.weight[param.nr_weight - 1] = atof(argv[i]);
                break;
            default:
                std::cerr << "Unknown option: -" << argv[i - 1][1] << "\n";
                exit_with_help();
        }
    }

    // determine filenames

    if (i >= argc) {
        exit_with_help();
    }

    strcpy(input_file_name, argv[i]);

    if (i < argc - 1) {
        strcpy(model_file_name, argv[i + 1]);
    } else {
        char *p = strrchr(argv[i], '/');
        if (p == NULL) {
            p = argv[i];
        } else {
            ++p;
        }
        sprintf(model_file_name, "%s.model", p);
    }

    return param;
}
