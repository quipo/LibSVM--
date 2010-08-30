/* 
 * File:   Matrix.h
 * Author: lorenzo
 *
 * Created on July 4, 2010, 5:30 PM
 */

#ifndef SVM_NODE_H
#define	SVM_NODE_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

typedef double Qfloat;

/**
 *
 */
class SVM_Node {
  public:
    SVM_Node() : index(0), value(0.0) { }
    SVM_Node(int idx, Qfloat val) : index(idx), value(val) { }
    ~SVM_Node() {}
    
    int index;
    Qfloat value;

    friend std::ostream & operator<<(std::ostream & os, const SVM_Node & p) {
        //return os << p.index << ":" << std::setiosflags(std::ios::fixed) << std::setprecision(8) << p.value << " ";
        return os << p.index << ":" << std::setprecision(8) << p.value << " ";
    }
    
    friend std::istream & operator>>(std::istream & in, SVM_Node & node) {
        std::string str;
        in >> str;
        std::stringstream ss(str);
        std::getline(ss, str, ':');
        std::istringstream (str) >> node.index;
        ss  >> node.value;
        return in;
    }
};

typedef std::vector<std::vector<SVM_Node> > node_matrix;

#endif	/* SVM_NODE_H */
