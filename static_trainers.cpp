/* 
 * File:   static_trainers.cpp
 * Author: lorenzo
 * 
 * Created on July 5, 2010, 2:00 PM
 */

#include <map>

#include "static_trainers.h"



// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
/**
 * 
 * @param SVM_Problem prob
 * @param nr_class_ret
 * @param label_ret        label name
 * @param start_ret        begin of each class
 * @param count_ret        #data of classes
 * @param perm             indices to the original data
 */
void svm_group_classes(const SVM_Problem & problem, unsigned int *nr_class_ret, std::vector<int> & label, std::vector<int> & start, std::vector<int> & count, std::vector<int> & perm) {
    unsigned int dimension = problem.dimension;
    unsigned int nr_class = 0;
    std::vector<int> data_label(dimension, 0);
    unsigned int i;

//    label.reserve(max_nr_class);
//    count.reserve(max_nr_class);
    label.clear();
    count.clear();
    
    for (i = 0; i < dimension; i++) {
        int this_label = static_cast<int>(problem.labels[i]);
        unsigned int c;
        for (c = 0; c < nr_class; c++) {
            if (this_label == label[c]) {
                ++count[c];
                break;
            }
        }
        data_label[i] = c;
        if (c == nr_class) {
            label.push_back(this_label);
            count.push_back(1);
            ++nr_class;
        }
    }

    start.reserve(nr_class);
    start.push_back(0);
    for (i = 1; i < nr_class; i++) {
        start.push_back(start[i - 1] + count[i - 1]);
    }
    for (i = 0; i < dimension; i++) {
        perm[start[data_label[i]]] = i;
        ++start[data_label[i]];
    }
    start[0] = 0;
    for (i = 1; i < nr_class; i++) {
        start[i] = start[i - 1] + count[i - 1];
    }

    *nr_class_ret = nr_class;
    data_label.clear();
}

/**
 * Groups the vectors in the original problem into smaller problems having the same class
 * 
 * @param problem
 *
 * @return std::map<label, problem>
std::map<int, SVM_Problem> split_problem_by_class(const SVM_Problem & problem) {
    std::map<int, SVM_Problem> m;

    for (unsigned int i=0; i<problem.dimension; ++i) {
        int label = static_cast<int>(problem.labels[i]);
        //std::cout << " label[" << i << "]=" << label;
        if (m.find(label) == m.end()) {
            SVM_Problem subproblem(0);
            m[label] = subproblem;
        }
        m[label].dimension++;
        m[label].labels.push_back(label);
        m[label].x.push_back(problem.x[i]);
    }

    return m;
}
 */
