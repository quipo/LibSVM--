/* 
 * File:   SVM_Parameter.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 5:30 PM
 */

#ifndef SVM_PARAMETER_H
#define	SVM_PARAMETER_H

/**
 * The parameters of a task
 */
class SVM_Parameter {
  public:
	int svm_type;	/* SVM Type */
	int kernel_type;/* Kernel Type */
	int degree;		/* for poly */
	double gamma;	/* for poly/gaussian/sigmoid */
	double coef0;	/* for poly/sigmoid */

	/* these are for training only */
	double cache_size; 	/* in MB */
	double eps;			/* stopping criteria (the lower eps, the better the approximation) */
	double C;			/* cost for C_SVC, SVDD, EPSILON_SVR and NU_SVR */
	unsigned int nr_weight;		/* for C_SVC */
	std::vector<int> weight_label;	/* for C_SVC */
	std::vector<double> weight;		/* for C_SVC */
	double nu;			/* for NU_SVC, ONE_CLASS, and NU_SVR */
	double p;			/* for EPSILON_SVR */
	int shrinking;		/* use the shrinking heuristics */
	int probability; 	/* do probability estimates */

    friend std::ostream & operator<<(std::ostream & os, const SVM_Parameter & param) {
        os << std::setprecision(3);
        os << "\n----[PARAMETER]-------------------";
        os << "\nsvm_type \t= "    << param.svm_type;    //<< " (" << svm_type_table[param.svm_type] << ")";
        os << "\nkernel_type \t= " << param.kernel_type; //<< " (" << kernel_type_table[param.kernel_type] << ")";
        os << "\ndegree \t\t= "    << param.degree;
        os << "\ngamma \t\t= "     << param.gamma;
        os << "\ncoef0 \t\t= "     << param.coef0;
        os << "\neps \t\t= "       << param.eps;
        os << "\nC \t\t= "         << param.C;
        os << "\nnr_weight \t= "   << param.nr_weight;
        os << "\nnu \t\t= "        << param.nu;
        os << "\np \t\t= "         << param.p;
        os << "\nshrinking \t= "   << param.shrinking;
        os << "\nprobability \t= " << param.probability;
        os << "\nweight_label[] \t= ";
        for (std::vector<int>::const_iterator it = param.weight_label.begin(), end = param.weight_label.end(); it != end; ++it) {
            os << " " << *it;
        }
        os << "\nweight[] \t= ";
        for (std::vector<double>::const_iterator it = param.weight.begin(), end = param.weight.end(); it != end; ++it) {
            os << " " << *it;
        }
        os << "\n----[END PARAMETER]-------------------";
        return os;
    }
};

#endif	/* SVM_PARAMETER_H */

