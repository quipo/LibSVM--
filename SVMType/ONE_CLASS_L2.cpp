/** 
 * File:   ONE_CLASS_L2.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 4, 2010, 4:54 PM
 */

#include "ONE_CLASS_L2.h"

std::vector<Qfloat> SVM_SVMType_ONE_CLASS_L2::get_Q(unsigned int column, unsigned int len) const {
    std::vector<Qfloat> data(len, 0);
    //unsigned int start = (cache->find(column, data)) ? data.size() : 0;
    unsigned int start = 0;
    if (start < len) {
        data.resize(len, 0);
        for (unsigned int j=start; j<len; ++j) {
            data[j] = kernel->kernel(column, j);
        }
        if (column >= start && column < len) {
            data[column] += 1/C;
        }
        //cache->insert(column, data);
    }
    return data;
}

SVM_SVMType_ONE_CLASS_L2::~SVM_SVMType_ONE_CLASS_L2() {
    delete cache;
    QD.clear();
}
