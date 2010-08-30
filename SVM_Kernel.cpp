/* 
 * File:   Kernel.cpp
 * Author: lorenzo
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "SVM_Kernel.h"
#include "Kernel/Exponential.h"

SVM_Kernel::SVM_Kernel() : dimension(0), param() {
    
}

SVM_Kernel::SVM_Kernel(const SVM_Parameter & param) :
    kernel_type(param.kernel_type),
    degree(param.degree),
    gamma(param.gamma),
    coef0(param.coef0),
    param(param)
{
    //clone(x_, x);
}

SVM_Kernel::SVM_Kernel(const node_matrix & x_, const SVM_Parameter & param) :
    //dimension(x.size()), //disabled, for some reason not working correctly
    x(x_),
    kernel_type(param.kernel_type),
    degree(param.degree),
    gamma(param.gamma),
    coef0(param.coef0),
    param(param)
{
    dimension = x.size();
    //clone(x_, x);
}

SVM_Kernel::~SVM_Kernel() {
	x.clear();
	x_square.clear();
}

//Qfloat SVM_Kernel::kernel(const unsigned int i, const unsigned int j) const {
//    return this->kernel(x[i], x[j], param);
//}

/*
Qfloat SVM_Kernel::kernel(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y, const SVM_Parameter & param) const {
    //override this method
    std::cout << "BASE_kernel ";
    if (param.svm_type) {
        //ugly hack to avoid "unused parameter 'param'" warning
    }
    return dot(x, y);
}
*/

std::vector<double> SVM_Kernel::get_diag(double c) const {
    std::vector<double> diag(dimension, -0.5);
    return (c) ? diag : diag; //ugly hack to avoid "unused parameter 'c'" warning
}

double SVM_Kernel::dot(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y) {
    std::vector<SVM_Node>::const_iterator px = x.begin(), py = y.begin();
    double sum = 0;
    while (px != x.end() && py != y.end()) {
        if (px->index == py->index) {
            sum += px->value * py->value;
            ++px;
            ++py;
        } else {
            if (px->index > py->index) {
                ++py;
            } else {
                ++px;
            }
        }
    }
    return sum;
}

double SVM_Kernel::dist_1(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y) {
    std::vector<SVM_Node>::const_iterator px = x.begin(), py = y.begin();
    double sum = 0;
    while (px != x.end() && py != y.end()) {
        if (px->index == py->index) {
            sum += fabs(px->value - py->value);
            ++px;
            ++py;
        } else {
            if (px->index > py->index) {
                sum += fabs(py->value);
                ++py;
            } else {
                sum += fabs(px->value);
                ++px;
            }
        }
    }
    while (px != x.end()) {
        sum += fabs(px->value);
        ++px;
    }
    while (py != y.end()) {
        sum += fabs(py->value);
        ++py;
    }
    return sum;
}

double SVM_Kernel::dist_2_sqr(const std::vector<SVM_Node> & x, const std::vector<SVM_Node> & y) {
    std::vector<SVM_Node>::const_iterator px = x.begin(), py = y.begin();
    double sum = 0;
    while (px != x.end() && py != y.end()) {
        if (px->index == py->index) {
            double d = px->value - py->value;
            sum += d * d;
            ++px;
            ++py;
        } else {
            if (px->index > py->index) {
                sum += py->value * py->value;
                ++py;
            } else {
                sum += px->value * px->value;
                ++px;
            }
        }

    }
    while (px != x.end()) {
        sum += px->value * px->value;
        ++px;
    }
    while (py != y.end()) {
        sum += py->value * py->value;
        ++py;
    }

    //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(1) << sum << " ";
    return (sum > 0.0 ? sum : 0.0);
}
