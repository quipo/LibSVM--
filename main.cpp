/* 
 * File:   main.cpp
 * Author: Lorenzo Alberton
 *
 * Created on July 3, 2010, 2:26 PM
 */

#include <cstdlib>
#include <iostream>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    std::cout << argc << "\n";
    for (int i=0; i<argc; ++i) {
        std::cout << " " << argv[i];
    }
    std::cout << "\n\n";
    return 0;
}

