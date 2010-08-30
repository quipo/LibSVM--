/* 
 * File:   NU.h
 * Author: Lorenzo Alberton
 *
 * Created on July 5, 2010, 2:00 PM
 */

#ifndef SVM_SOLVER_NU_H
#define	SVM_SOLVER_NU_H

#include "../SVM_Solver.h"

/**
 *  Solver for nu-svm classification and regression
 * 
 * additional constraint: e^T \alpha = constant
 */
class SVM_Solver_NU : public SVM_Solver
{
  public:
    SVM_Solver_NU() {}
    
  protected:
    int select_working_set(int &i, int &j);
    double calculate_rho();
    bool be_shrunk(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4);
    void do_shrinking();
};

#endif	/* SVM_SOLVER_NU_H */
