/* 
 * File:   Matrix.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 5:30 PM
 */

#ifndef SVM_PROBLEM_H
#define	SVM_PROBLEM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "SVM_Node.h"
#include "SVM_Parameter.h"



/**
 * SVM Type: Classification, Classification with L2SVM, nu-Classification,
 * One Class Classification, epsilon Regression, nu Regression,
 * Support Vector Domain Description (SVDD), SVDD with L2SVM
 *
 * C_SVC_L2, SVDD and SVDD_L2 were backported from the code available at
 *
 *    http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/#18
 *
 * written by Leland Wang and extended by Holger Froehlich.
 *
 * --
 * Vincenzo Russo (vincenzo.russo@neminis.org)
 */
enum { C_SVC, C_SVC_L2, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR, SVDD, SVDD_L2 };

/**
 * Kernel Type: Linear, Polynomial, Gaussian, Sigmoid, Stump, Perceptron, Laplace, Exponential, Precomputed
 *
 * Stump, Perceptron, Laplace, Exponential were backported from the LIBSVM extension by Hsuan-Tien Lin
 *
 * For his version of the software, see http://www.work.caltech.edu/~htlin/program/libsvm/#infensemble
 * For explanation about the four above kernels, see his Master's Thesis:
 *   "H.-T. Lin. Infinite Ensemble Learning with Support Vector Machines.
 *   Master's Thesis, California Institute of Technology, May 2005."
 *   http://www.work.caltech.edu/~htlin/publication/
 *
 * --
 * Vincenzo Russo (vincenzo.russo@neminis.org)
 */
enum KernelTypes { LINEAR, POLYNOMIAL, GAUSSIAN, SIGMOID, STUMP, PERCEPTRON, LAPLACE, EXPONENTIAL, PRECOMPUTED };	/* kernel_type */


/*
// class generator:
struct c_sequence {
  int current;
  c_sequence() {
      current = 0;
  }
  int operator()() {
      return ++current;
  }
} SequenceNumber;
*/

/**
 * 	The problem: the dimension (l), the labels (*y) and the patterns (**x)
 */
class SVM_Problem
{
  public:
    unsigned int dimension;
	std::vector<Qfloat> labels; 
	node_matrix x; //patterns

    SVM_Problem(const unsigned int dim) : dimension(dim) {
        x.reserve(dim);
        labels.reserve(dim);
    }
    
    SVM_Problem(const char *filename, SVM_Parameter & param);
    ~SVM_Problem();
	void add(const SVM_Problem & p);

    friend std::ostream & operator<<(std::ostream & os, const SVM_Problem & prob) {
        os << std::setprecision(3);
        os << "\n----[PROBLEM]-------------------";
        os << "\ndimension \t= "    << prob.dimension;
//        os << "\nlabels[] \t= ";
//        for (std::vector<Qfloat>::const_iterator l_it = prob.labels.begin(), l_end = prob.labels.end(); l_it != l_end; ++l_it) {
//            os << " " << *l_it;
//        }
//        os << "\nx[] \t= ";
//        for (std::vector<std::vector<SVM_Node> >::const_iterator m_it = prob.x.begin(), m_end = prob.x.end(); m_it != m_end; ++m_it) {
//            os << "\n\t";
//            for (std::vector<SVM_Node>::const_iterator n_it = m_it->begin(), n_end = m_it->end(); n_it != n_end; ++n_it) {
//                os << *n_it;
//            }
//        }
        os << "\nx[] \t= ";
        std::vector<Qfloat>::const_iterator l_it = prob.labels.begin(), l_end = prob.labels.end();
        for (std::vector<std::vector<SVM_Node> >::const_iterator m_it = prob.x.begin(), m_end = prob.x.end(); m_it != m_end && l_it != l_end; ++m_it, ++l_it) {
            os.setf(std::ios::showpos);
            os << "\n" << *l_it << " ";
            os.unsetf(std::ios::showpos);
            for (std::vector<SVM_Node>::const_iterator n_it = m_it->begin(), n_end = m_it->end(); n_it != n_end; ++n_it) {
                os << *n_it;
            }
        }
        os << "\n----[END PROBLEM]-------------------";
        return os;
    }
};


#endif	/* SVM_PROBLEM_H */

