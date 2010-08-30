#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "SVM_Node.h"

void exit_with_help() {
    std::cout <<
        "Usage: svm-scale [options] data_filename\n"
        "options:\n"
        "-l lower : x scaling lower limit (default -1)\n"
        "-u upper : x scaling upper limit (default +1)\n"
        "-y y_lower y_upper : y scaling limits (default: no y scaling)\n"
        "-s save_filename : save scaling parameters to save_filename\n"
        "-r restore_filename : restore scaling parameters from restore_filename\n"
        ;
    exit(1);
}

double lower = -1.0, upper = 1.0, y_lower, y_upper;
int y_scaling = 0;
double y_max = -DBL_MAX;
double y_min = DBL_MAX;
int max_index;
long int num_nonzeros = 0;
long int new_num_nonzeros = 0;
double eps = 1e-5;

void output_target(double value);
void output(const int index, double value, const std::vector<double> & feature_min, const std::vector<double> & feature_max);


int main(int argc, char **argv) {
    int i;
    char *save_filename = NULL;
    char *restore_filename = NULL;

    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-') break;
        ++i;
        switch (argv[i - 1][1]) {
            case 'l': lower = atof(argv[i]);
                break;
            case 'u': upper = atof(argv[i]);
                break;
            case 'y':
                y_lower = atof(argv[i]);
                ++i;
                y_upper = atof(argv[i]);
                y_scaling = 1;
                break;
            case 's': save_filename = argv[i];
                break;
            case 'r': restore_filename = argv[i];
                break;
            default:
                std::cerr << "unknown option '" << argv[i - 1][1] << "'\n";
                exit_with_help();
        }
    }

    if (!(upper > lower) || (y_scaling && !(y_upper > y_lower))) {
        std::cerr << "inconsistent lower/upper specification\n";
        exit(1);
    }

    if (restore_filename && save_filename) {
        std::cerr << "cannot use -r and -s simultaneously\n";
        exit(1);
    }

    if (argc != i + 1) {
        exit_with_help();
    }

/*
#define SKIP_TARGET\
	while(isspace(*p)) ++p;\
	while(!isspace(*p)) ++p;

#define SKIP_ELEMENT\
	while(*p!=':') ++p;\
	++p;\
	while(isspace(*p)) ++p;\
	while(*p && !isspace(*p)) ++p;
*/
    /* assumption: min index of attributes is 1 */
    /* pass 1: find out max index of attributes */
    max_index = 0;

    std::string line;

    std::ifstream fp_restore;
    if (restore_filename) {
        fp_restore.open(restore_filename);
        if (!fp_restore.is_open()) {
            std::cerr << "can't open file " << restore_filename << "\n";
            exit(1);
        }

        int c = fp_restore.peek();
        if (c == 'y') {
            std::getline(fp_restore, line);
            std::getline(fp_restore, line);
            std::getline(fp_restore, line);
        }
        std::getline(fp_restore, line);
        std::getline(fp_restore, line);

        int idx;
        while (std::getline(fp_restore, line)) {
            std::istringstream (line) >> idx;
            max_index = std::max(idx, max_index);
        }
        fp_restore.clear();
        fp_restore.seekg(0); //rewind
    }

    std::ifstream fp(argv[i]);
    if (!fp.is_open()) {
        std::cerr << "can't open file " << argv[i] << "\n";
        exit(1);
    }
    std::string dummy;
    while (std::getline(fp, line)) {
        std::istringstream is(line);
        is >> dummy; // SKIP_TARGET

        SVM_Node node;
        while (is >> node) {
            max_index = std::max(max_index, node.index);
            //is >> dummy; // SKIP_ELEMENT
            ++num_nonzeros;
        }
    }
    fp.clear();
    fp.seekg(0, std::ios_base::beg); //rewind

    std::vector<double> feature_max(max_index + 1, -DBL_MAX);
    std::vector<double> feature_min(max_index + 1, DBL_MAX);

    /* pass 2: find out min/max value */
    while (std::getline(fp, line)) {
        std::istringstream is(line);
        int next_index = 1;
        double target;

        is >> target;
        y_max = std::max(y_max, target);
        y_min = std::min(y_min, target);

        // is >> dummy; // SKIP_TARGET

        SVM_Node node;
        while (is >> node) {
            for (i = next_index; i < node.index; i++) {
                feature_max[i] = std::max(feature_max[i], 0.0);
                feature_min[i] = std::min(feature_min[i], 0.0);
            }

            feature_max[node.index] = std::max(feature_max[node.index], node.value);
            feature_min[node.index] = std::min(feature_min[node.index], node.value);

            //is >> dummy; // SKIP_ELEMENT
            next_index = node.index + 1;
        }

        for (i = next_index; i <= max_index; i++) {
            feature_max[i] = std::max(feature_max[i], 0.0);
            feature_min[i] = std::min(feature_min[i], 0.0);
        }
    }

//    std::cout << "\nFEATURES[]=\n";
//    for (int k=0; k<max_index; ++k) {
//        std::cout << std::setprecision(16) << feature_max[k] << " ";
//        if (k % 20 == 19) std::cout << "\n";
//    }

    fp.clear();
    fp.seekg(0); //rewind

    /* pass 2.5: save/restore feature_min/feature_max */

    if (restore_filename) {
        /* fp_restore rewinded in finding max_index */
        int idx;
        double fmin, fmax;

        int c = fp_restore.peek();
        if (c == 'y') {
            std::getline(fp, line);
            std::istringstream(line) >> y_lower >> y_upper;
            std::getline(fp, line);
            std::istringstream(line) >> y_min >> y_max;
            y_scaling = 1;
        }

        c = fp_restore.peek();
        if (c == 'x') {
            std::getline(fp, line);
            std::istringstream(line) >> lower >> upper;
            while (std::getline(fp, line)) {
                std::istringstream(line) >> idx >> fmin >> fmax;
                if (idx <= max_index) {
                    feature_min[idx] = fmin;
                    feature_max[idx] = fmax;
                }
            }
        }
        fp_restore.close();
    }

    if (save_filename) {
        std::ofstream fp_save(save_filename);
        if (!fp_save.is_open()) {
            std::cerr << "can't open file " << save_filename << "\n";
            exit(1);
        }
        fp_save << std::setprecision(16);

        if (y_scaling) {
            fp_save << "y\n";
            fp_save << y_lower << " " << y_upper << "\n";
            fp_save << y_min << " " << y_max << "\n";
        }
        fp_save << "x\n";
        fp_save << lower << " " << upper << "\n";
        for (i = 1; i <= max_index; i++) {
            if (feature_min[i] != feature_max[i]) {
                fp_save << i << " " << feature_min[i] << " " << feature_max[i] << "\n";
            }
        }
        fp_save.close();
    }

    /* pass 3: scale */
    while (std::getline(fp, line)) {
        std::istringstream is(line);
        int next_index = 1;
        double target;

        is >> target;
        output_target(target);

        //is >> dummy; // SKIP_TARGET

        SVM_Node node;
        while (is >> node) {
            for (i = next_index; i < node.index; ++i) {
                output(i, 0, feature_min, feature_max);
            }

            output(node.index, node.value, feature_min, feature_max);

            //is >> dummy; // SKIP_ELEMENT
            next_index = node.index + 1;
        }

        for (i = next_index; i <= max_index; ++i) {
            output(i, 0, feature_min, feature_max);
        }

        std::cout << "\n";
    }

    if (new_num_nonzeros > num_nonzeros) {
        std::cerr <<
            "Warning: original #nonzeros " << num_nonzeros << "\n"
            "         new      #nonzeros " << new_num_nonzeros << "\n"
            "Use -l 0 if many original feature values are zeros\n";
    }

    //free(line);
    feature_max.clear();
    feature_min.clear();
    fp.close();
    return 0;
}

void output_target(double value) {
    if (y_scaling) {
        if (fabs(value - y_min) < eps) {
            value = y_lower;
        } else if (fabs(value - y_max) < eps) {
            value = y_upper;
        } else {
            value = y_lower + (y_upper - y_lower) *
            (value - y_min) / (y_max - y_min);
        }
    }
    std::cout.setf(std::ios::showpos);
    std::cout << std::setprecision(16) << value << " ";
    std::cout.unsetf(std::ios::showpos);
}

void output(const int index, double value, const std::vector<double> & feature_min, const std::vector<double> & feature_max) {
    /* skip single-valued attribute */
    if (feature_max[index] == feature_min[index]) {
        return;
    }

    if (value == feature_min[index]) {
        value = lower;
    } else if (value == feature_max[index]) {
        value = upper;
    } else {
        value = lower + (upper - lower) *
            (value - feature_min[index]) /
            (feature_max[index] - feature_min[index]);
    }

    if (value != 0.0) {
        std::cout << index << ":" << std::setprecision(16) << value << " ";
        ++new_num_nonzeros;
    }
}
