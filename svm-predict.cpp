#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "static_trainers.h"
#include "interface_functions.h"


void exit_input_error(int line_num) {
    std::cerr << "Wrong input format at line " << line_num << "\n";
    exit(1);
}

void predict(SVM_Model & model, int predict_probability, std::istream & input, std::ostream & output) {
    int correct = 0;
    int total = 0;
    double error = 0;
    double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

    int svm_type = model.param.svm_type;
    int nr_class = model.nr_class;
    std::vector<double> prob_estimates;
    int j;

    SVM_SVMType * svmtype = svmtype_factory(model.param);
    
    std::cout << std::setprecision(4);
    output << std::setprecision(16);

    if (predict_probability) {
        if (svm_type == NU_SVR || svm_type == EPSILON_SVR) {
            std::cout << "Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=" << model.get_svr_probability() << "\n";
        } else {
            prob_estimates.reserve(nr_class);
            output << "labels";
            for (std::vector<int>::const_iterator it = model.label.begin(); it != model.label.end(); ++it) {
                output << " " << *it;
            }
            output << "\n";
        }
    }
    
    std::string line;
    while (std::getline(input, line)) {
        std::vector<SVM_Node> x;
        double target_label, predict_label;
        int inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0

        std::istringstream is(line);
        is >> target_label;
        //std::cout << "\n";
        SVM_Node node;
        while (is.good()) {
            is >> node;
            if (!is.good()) {
                break;
            }
          //  std::cout << node; std::cout.flush();
            if (node.index <= inst_max_index) {
                exit_input_error(total + 1);
            } else {
                inst_max_index = node.index;
            }
            x.push_back(node);
        }

        if (predict_probability && (svm_type == C_SVC || svm_type == NU_SVC)) {
            predict_label = svmtype->svm_predict_probability(model, x, prob_estimates);
            output << predict_label;
            for (j = 0; j < nr_class; j++) {
                output << " " << prob_estimates[j];
            }
            output <<"\n";
        } else {
            predict_label = svmtype->svm_predict(model, x);
            output << predict_label << "\n";
        }
        x.clear();
        
        if (predict_label == target_label) {
            ++correct;
        }
        error += (predict_label - target_label) * (predict_label - target_label);
        sump  += predict_label;
        sumt  += target_label;
        sumpp += predict_label * predict_label;
        sumtt += target_label * target_label;
        sumpt += predict_label * target_label;
        ++total;
 
    }
    if (svm_type == NU_SVR || svm_type == EPSILON_SVR) {
        std::cout << "Mean squared error = " << (error / total) << " (regression)\n";
        std::cout << "Squared correlation coefficient = "
            << (((total * sumpt - sump * sumt)*(total * sumpt - sump * sumt)) /
                ((total * sumpp - sump * sump)*(total * sumtt - sumt * sumt)))
            << "(regression)\n";
    } else {
        std::cout << "Accuracy = " << ((double) correct / total * 100) << "% (" << correct << "/" << total << ") (classification)\n";
    }
    if (predict_probability) {
        prob_estimates.clear();
    }
}

void exit_with_help() {
    std::cout <<
        "Usage: svm-predict [options] test_file model_file output_file\n"
        "options:\n"
        "-b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); for one-class SVM only 0 is supported\n"
    ;
    exit(1);
}

int main(int argc, char **argv) {
    int i;
    int predict_probability = 0;
    
    // parse options
    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-') break;
        ++i;
        switch (argv[i - 1][1]) {
            case 'b':
                predict_probability = atoi(argv[i]);
                break;
            default:
                std::cerr << "Unknown option: -" << argv[i - 1][1] << "\n";
                exit_with_help();
        }
    }
    if (i >= argc - 2) {
        exit_with_help();
    }

    std::ifstream input;
    std::ofstream output;
    input.open(argv[i]);
    if (!input.is_open()) {
        std::cerr << "can't open input file " << argv[i] << "\n";
        exit(1);
    }

    output.open(argv[i + 2]);
    if (!output.is_open()) {
        std::cerr << "can't open output file " << argv[i + 2] << "\n";
        exit(1);
    }

    //TODO add exceptions to this constructor in case of failure
    SVM_Model model(argv[i + 1]);
//    std::cout << "KERNEL: " << model.param.kernel_type;
//    std::cout.flush();

    
//    if ((model = svm_load_model(argv[i + 1])) == 0) {
//        fprintf(stderr, "can't open model file %s\n", argv[i + 1]);
//        exit(1);
//    }

    if (predict_probability) {
        if (svm_check_probability_model(model) == 0) {
            std::cerr << "Model does not support probabiliy estimates\n";
            exit(1);
        }
    } else {
        if (svm_check_probability_model(model) != 0) {
            std::cout << "Model supports probability estimates, but disabled in prediction.\n";
        }
    }
    predict(model, predict_probability, input, output);
    input.close();
    output.close();
    return 0;
}
