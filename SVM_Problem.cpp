/* 
 * File:   SVM_Problem.cpp
 * Author: lorenzo
 * 
 * Created on July 5, 2010, 2:00 PM
 */

#include "SVM_Problem.h"

/**
 * Append the labels and support vectors from another problem
 * @param p
 */
void SVM_Problem::add(const SVM_Problem & p) {
    x.insert(x.end(), p.x.begin(), p.x.end());
    labels.insert(labels.end(), p.labels.begin(), p.labels.end());
    dimension = labels.size();
}

SVM_Problem::~SVM_Problem() { }

/**
 * Read in a problem (in svmlight format)
 *
 * @param filename
 * @param param
 */
SVM_Problem::SVM_Problem(const char *filename, SVM_Parameter & param) {
    std::ifstream is(filename, std::ios::binary);
    if (!is.is_open()) {
        std::cerr << "can't open input file " << filename << "\n";
        exit(1);
    }

    dimension = 0;
    int max_index = 0;
    unsigned int i = 0;
    int inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
    std::string line;
    while (std::getline(is, line)) {
        ++i;
        std::istringstream p(line);

        int label;
        p >> label;
        labels.push_back(label);
//        if (endptr == label) {
//            exit_input_error(i);
//        }

        std::vector<SVM_Node> x_space;
        x_space.reserve(dimension);

//        std::string vpair;
//        while (p >> vpair) {
//            std::string index;
//
//            SVM_Node node;
//            std::stringstream ss(vpair);
//            //idx = strtok(NULL, ":");
//            //val = strtok(NULL, " \t");
//            std::getline(ss, index, ':');
//            std::istringstream idx(index);
//            idx >> node.index;
//            ss >> node.value;
        SVM_Node node;
        while (p >> node) {
            if (node.index > inst_max_index) {
                inst_max_index = node.index;
            }
            //std::cout << node << ' ';
            x_space.push_back(node);
        }
        x.push_back(x_space);

        //std::cout << prob.x;

        if (inst_max_index > max_index) {
            max_index = inst_max_index;
        }
    }

    dimension = x.size();

//    for (node_matrix::const_iterator it = prob.x.begin(); it != prob.x.end(); ++it) {
//        for (std::vector<SVM_Node>::const_iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
//            std::cout << *it2;
//        }
//        std::cout << "\n";
//    }
//    exit(0);

    if (param.gamma == 0 && max_index > 0) {
        param.gamma = 1.0 / max_index;
    }

    if (param.kernel_type == PRECOMPUTED) {
        for (i = 0; i < dimension; i++) {
            if (x[i][0].index != 0) {
                std::cerr << "Wrong input format: first column must be 0:sample_serial_number\n";
                exit(1);
            }
            if (static_cast<int>(x[i][0].value) <= 0 || static_cast<int>(x[i][0].value) > max_index) {
                std::cerr << "Wrong input format: sample_serial_number out of range\n";
                exit(1);
            }
        }
    }

    is.close();
}
