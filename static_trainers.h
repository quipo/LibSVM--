/* 
 * File:   SVM_Solver.h
 * Author: lorenzo
 *
 * Created on July 5, 2010, 2:00 PM
 */

#ifndef STATIC_TRAINERS_H
#define	STATIC_TRAINERS_H

#include <map>
#include <vector>
#include <algorithm>
#include <functional>

#include "SVM_Kernel.h"
#include "SVM_Model.h"
#include "SVM_Parameter.h"
#include "SVM_Problem.h"
#include "SVM_Solver.h"
#include "interface_functions.h"

decision_function svm_train_one(
    const SVM_Problem & prob,
    const SVM_Parameter & param,
    double Cp,
    double Cn
);

void svm_group_classes(const SVM_Problem & problem, unsigned int *nr_class_ret, std::vector<int> & label, std::vector<int> & start, std::vector<int> & count, std::vector<int> & perm);
//std::map<int, SVM_Problem> split_problem_by_class(const SVM_Problem & problem);

#endif	/* STATIC_TRAINERS_H */

