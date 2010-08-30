/* 
 * File:   SVM_Kernel.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 4:54 PM
 */

#ifndef SVM_KERNEL_H
#define	SVM_KERNEL_H

#include <math.h>
#include "SVM_Matrix.h"
#include "SVM_Parameter.h"

#ifndef min
template <class T> static inline T min(T x, T y) { return (x<y)?x:y; }
#endif
#ifndef max
template <class T> static inline T max(T x, T y) { return (x>y)?x:y; }
#endif
//template <class T> static inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
/*
template <class S, class T> static inline void clone(T*& dst, S* src, int n) {
	dst = new T[n];
	memcpy((void *)dst,(void *)src, sizeof(T)*n);
}
*/

class SVM_Kernel : public SVM_Matrix {
  public:
    SVM_Kernel();
    SVM_Kernel(const SVM_Parameter & param);
	SVM_Kernel(const node_matrix & x, const SVM_Parameter & param);
	virtual ~SVM_Kernel();

	Qfloat kernel(const unsigned int i, const unsigned int j) const {
        return this->kernel(x[i], x[j], param);
    }
    virtual Qfloat kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const = 0;

    std::vector<double> get_diag(double c = 0) const;


	//virtual
    void swap_index(unsigned int i, unsigned int j) {
		std::swap(this->x[i], this->x[j]);
		if (this->x_square.size()) {
            std::swap(this->x_square[i], this->x_square[j]);
        }
    }

    unsigned int dimension;

  protected:
    node_matrix x;
    std::vector<Qfloat> x_square;

    // svm_parameter
    /* const */ int kernel_type;
    /* const */ int degree;
    /* const */ double gamma;
    /* const */ double coef0;
    /* const */ SVM_Parameter param;

    static double        dot(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y);
    static double     dist_1(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y);
    static double dist_2_sqr(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y);

    inline double dist_2_sqr(int i, int j) const {
            double sum = this->x_square[i] + this->x_square[j] - 2 * dot(this->x[i], this->x[j]);
            return (sum > 0.0 ? sum : 0.0);
    }
};

#endif	/* SVM_KERNEL_H */
