/* 
 * File:   SVM_Matrix.h
 * Author: Lorenzo Alberton
 *
 * Created on July 4, 2010, 5:30 PM
 */

#ifndef SVM_MATRIX_H
#define	SVM_MATRIX_H

#include <vector>
#include "SVM_Node.h"

/**
 *
 * @param column
 * @param len
 * @return
 */
class SVM_Matrix {
  public:
	//virtual void swap(unsigned int i, unsigned int j) const = 0;
	virtual ~SVM_Matrix() {}
};

#endif	/* SVM_MATRIX_H */

