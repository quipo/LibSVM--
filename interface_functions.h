/* 
 * File:   interface_functions.h
 * Author: Lorenzo Alberton
 *
 * Created on July 5, 2010, 2:00 PM
 */

#ifndef INTERFACE_FUNCTIONS_H
#define	INTERFACE_FUNCTIONS_H

#include <math.h>
#include <vector>
#include <string>
#include <algorithm>

#include "static_trainers.h"

//template <class T> struct sequence : public std::unary_function<T, void> {
//    void operator()(T & x) { x++; }
//};


SVM_SVMType * svmtype_factory(const SVM_Parameter & param);
SVM_SVMType * svmtype_factory(const SVM_Problem & prob, const SVM_Parameter & param);

void svm_cross_validation(const SVM_Problem & prob, const SVM_Parameter & param, unsigned int nr_fold, std::vector<double> & target);
const char * svm_check_parameter(const SVM_Problem & prob, const SVM_Parameter & param);
int svm_check_probability_model(const SVM_Model & model);
int count_elements(std::ifstream & is);

#endif	/* INTERFACE_FUNCTIONS_H */

