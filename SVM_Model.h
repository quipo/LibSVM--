/* 
 * File:   Kernel.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_MODEL_H
#define	SVM_MODEL_H

#include "SVM_Node.h"
#include "SVM_Parameter.h"
#include "SVM_Kernel.h"
#include "Cache.h"
#include "LRUCache.h"
#include "Kernel/Linear.h"
#include "Kernel/Polynomial.h"
#include "Kernel/Gaussian.h"
#include "Kernel/Sigmoid.h"
#include "Kernel/Stump.h"
#include "Kernel/Perceptron.h"
#include "Kernel/Laplace.h"
#include "Kernel/Exponential.h"
#include "Kernel/Precomputed.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

/**
 * The decision function
 */
class decision_function
{
public:
	//double *alpha;	// Lagrangian multipliers
    std::vector<double> alpha;  // Lagrangian multipliers
	//double rho;
    std::vector<double> rho; //constants in decision functions

	//unsigned int *BSV_idx;	// indices of the BSVs in the original dataset
	//unsigned int *SV_idx;	// indices of the SVs in the original dataset

    std::vector<unsigned int> BSV_idx;  // indices of the BSVs in the original dataset
    std::vector<unsigned int> SV_idx;   // indices of the SVs in the original dataset

	unsigned int dimension;		// number of total SVs (SVs + BSVs)
	unsigned int lbsv;	// number of BSVs
};

/**
 * The model yielded by the training of an SVM
 * */
class SVM_Model : public decision_function {
public:
    enum { LINEAR, POLYNOMIAL, GAUSSIAN, SIGMOID, STUMP, PERCEPTRON, LAPLACE, EXPONENTIAL, PRECOMPUTED };	/* kernel_type */

public:
    SVM_Model(unsigned int dim) : dimension(dim) { };
    SVM_Model(const SVM_Parameter & param, unsigned int dim);
    SVM_Model(const char *model_file_name);
    ~SVM_Model();

    int save(const char *model_file_name) const;

    double get_svr_probability() {
        if (probA.size() > 0) {
            return probA[0];
        }
        std::cerr << "Model doesn't contain information for SVR probability inference\n";
        return 0.0;
    }

	SVM_Parameter param;	// parameter
	unsigned int nr_class;	// number of classes, = 2 in regression/one class svm
	unsigned int dimension;	// total #SV
	node_matrix SV;	// SVs (SV[l])

    SVM_Kernel * kernel;

	//std::vector<unsigned int> BSV_idx;	// indices of the BSVs in the original dataset
	//std::vector<unsigned int> SV_idx;	// indices of the SVs in the original dataset

	//unsigned int lbsv;		// total #BSV

	//double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[k-1][l])
    std::vector<std::vector<double> > sv_coef; // sv_coef[class][coefficient]
	//double *rho;		// constants in decision functions (rho[k*(k-1)/2])
    //std::vector<double> rho;
	
    std::vector<double> probA;  // pairwise probability information
    std::vector<double> probB;

	// for classification only
	std::vector<int> label;         // label of each class (label[k])
	std::vector<unsigned int> nSV;	// number of SVs for each class (nSV[k])
                        // nSV[0] + nSV[1] + ... + nSV[k-1] = l
	// XXX
	unsigned int free_sv;	// 1 if svm_model is created by svm_load_model
                            // 0 if svm_model is created by svm_train
};

#endif	/* SVM_MODEL_H */

