/* 
 * File:   SVM_Solver.h
 * Author: Lorenzo Alberton
 *
 * Created on July 5, 2010, 2:00 PM
 */

#ifndef SVM_SOLVER_H
#define	SVM_SOLVER_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "SVM_Solutioninfo.h"
#include "SVM_Matrix.h"
#include "SVM_SVMType.h"

#define INF HUGE_VAL
#define TAU 1e-12

//#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
/* static int sequence() {
    static int x = 0;
    return x++;
}
 */

/* static void info(const char *fmt,...) {
    std::cout << fmt;
//    char buf[BUFSIZ];
//	va_list ap;
//	va_start(ap, fmt);
//	vsprintf(buf, fmt, ap);
//	va_end(ap);
//	(*svm_print_string)(buf);
}
*/
class SVM_Solver {
public:
	SVM_Solver() {};
	virtual ~SVM_Solver() {};

	void solve(
        unsigned int l_,
        SVM_SVMType * Q,
        const std::vector<double> & p_,
        const std::vector<schar> & y_,
        std::vector<double> & alpha_,
        double Cp,
        double Cn,
        double eps,
        SVM_SolutionInfo & si,
        int shrinking
    );

protected:
	unsigned int active_size;
	std::vector<schar> y;
	std::vector<double> G;      // gradient of objective function
	std::vector<double> G_bar;	// gradient, if we treat free variables as 0
	enum AlphaStatusEnum { LOWER_BOUND, UPPER_BOUND, FREE };
	std::vector<AlphaStatusEnum> alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
	std::vector<double> alpha;
	//const
    SVM_SVMType * Q;
	//const
    std::vector<Qfloat> QD;
	SVM_SolutionInfo si;
	double eps;
	double Cp, Cn;
	std::vector<double> p;
	std::vector<int> active_set;
	unsigned int dimension;
	bool unshrink;	// XXX

	double get_C(int i) {
		return (static_cast<int>(y[i]) > 0) ? Cp : Cn;
	}

    void init_alpha_status(unsigned int size) {
        alpha_status.assign(size, FREE);
    }

	void update_alpha_status(unsigned int i) {
        if (alpha[i] >= get_C(i)) {
			alpha_status[i] = UPPER_BOUND;
        } else if (alpha[i] <= 0) {
			alpha_status[i] = LOWER_BOUND;
        } else {
            alpha_status[i] = FREE;
        }
	}

	bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
	bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
	bool is_free(int i)        { return alpha_status[i] == FREE; }
    
	void swap_index(unsigned int i, unsigned int j);
	void reconstruct_gradient();
	virtual int select_working_set(unsigned int &i, unsigned int &j);
	virtual double calculate_rho();
	virtual void do_shrinking();
private:
	bool be_shrunk(int i, double Gmax1, double Gmax2);
};


#endif	/* SVM_SOLVER_H */

