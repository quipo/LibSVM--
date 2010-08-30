/* 
 * File:   SVM_Solver.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 5, 2010, 2:00 PM
 */

#include "SVM_Solver.h"

void SVM_Solver::swap_index(unsigned int i, unsigned int j) {
    Q->swap_index(i, j);
    std::swap(y[i], y[j]);
    std::swap(G[i], G[j]);
    std::swap(alpha_status[i], alpha_status[j]);
    std::swap(alpha[i], alpha[j]);
    std::swap(p[i], p[j]);
    std::swap(active_set[i], active_set[j]);
    std::swap(G_bar[i], G_bar[j]);
}

/**
 * reconstruct inactive elements of G from G_bar and free variables
 */
void SVM_Solver::reconstruct_gradient() {
    //std::cout << "\n==[Reconstructing Gradient]==";
    if (active_size == dimension) {
        //std::cout << "\nactive_size == dimension == " << dimension << "\n\n";
        return;
    }
    //std::cout << "\nactive_size = " << active_size << ", dimension = " << dimension;

    unsigned int i, j;
    unsigned int nr_free = 0;
    
    for (j = active_size; j < dimension; ++j) {
        G[j] = G_bar[j] + p[j];
    }

    for (j = 0; j < active_size; ++j) {
        if (is_free(j)) {
            ++nr_free;
        }
    }

    if (2 * nr_free < active_size) {
        std::cout << "\nWarning: using -h 0 may be faster\n";
    }

    if (nr_free * dimension > 2 * active_size * (dimension - active_size)) {
        //std::cout << "\nBranch 1";
        for (i = active_size; i < dimension; ++i) {
            const std::vector<Qfloat> Q_i = Q->get_Q(i, active_size);
            for (j = 0; j < active_size; ++j) {
                if (is_free(j)) {
                    G[i] += alpha[j] * Q_i[j];
                }
            }
        }
    } else {
        //std::cout << "\nBranch 2";
        for (i = 0; i < active_size; ++i) {
            if (is_free(i)) {
                const std::vector<Qfloat> Q_i = Q->get_Q(i, dimension);
                double alpha_i = alpha[i];
                for (j = active_size; j < dimension; ++j) {
                    G[j] += alpha_i * Q_i[j];
                }
            }
        }
    }
}

void SVM_Solver::solve(
    unsigned int dim,
    SVM_SVMType * Q_,
    const std::vector<double> & p_,
    const std::vector<schar> & y_,
    std::vector<double> & alpha_,
    double Cp_,
    double Cn_,
    double eps_,
    SVM_SolutionInfo & si_,
    int shrinking
) {
    //std::cout << "\n\nsolve(" << dim << ")\n\n";
    dimension = dim;
    Q = Q_;
    QD = Q->get_QD();
    p.assign(p_.begin(), p_.end());
    y.assign(y_.begin(), y_.end());
    alpha.assign(alpha_.begin(), alpha_.end());
    Cp = Cp_;
    Cn = Cn_;
    eps = eps_;
    si = si_;
    unshrink = false;

    //std::cout << "Cp, Cn = " << Cp << ", " << Cn << "\n";

//    for (unsigned int t = 0; t < dim; t++) {
//        std::cout << (int)y[t] << " ";
//    }

    // initialize alpha_status
    // initialize active set (for shrinking)
    {
        init_alpha_status(dim);
        active_set.reserve(dim);
        G.reserve(dim);
        //std::generate(active_set.begin(), active_set.end(), SequenceNumber);
        for (unsigned int i = 0; i < dim; ++i) {
            update_alpha_status(i);
            active_set.push_back(i);
            G.push_back(p[i]);
        }
        active_size = dim;
    }



//    std::cout << "\nG[" << G.size() << "] = ";
//    for (std::vector<double>::const_iterator it = G.begin(); it != G.end(); ++it) {
//        std::cout << *it << " ";
//    }

//    std::cout << "\n\nalpha[] = \n" << std::setprecision(3);
//    for (unsigned int t = 0; t < active_size; t++) {
//        std::cout << alpha[t] << " ";
//        if (t % 20 == 19) std::cout << "\n";
//    }
//    std::cout << "\n\nalpha_status[] = \n";
//    for (unsigned int t = 0; t < active_size; t++) {
//        std::cout << (int)alpha_status[t] << " ";
//        if (t % 20 == 19) std::cout << "\n";
//    }
//    exit(0);


    // initialize gradient
    {
        G_bar.assign(dim, 0.0);
        for (unsigned int i = 0; i < dim; ++i) {
            if (!is_lower_bound(i)) {
                const std::vector<Qfloat> Q_i = Q->get_Q(i, dim);
                double alpha_i = alpha[i];
                unsigned int j =0;

//                if (alpha_i !=0 && Q_i[j] != 0) {
//                    std::cout << "(" << i << "," << j << ") ";
//                    std::cout << "(" << alpha_i << "," << Q_i[j] << ") ";
//                    exit(0);
//                }
                for (j = 0; j < dim; ++j) {
                    G[j] += alpha_i * Q_i[j];
                }
//                if (alpha_i !=0 && Q_i[j] != 0) {
//                    std::cout << "(" << i << "," << j << ") ";
//                    std::cout << "(" << alpha_i << "," << Q_i[j] << ") ";
//                    exit(0);
//                }
                if (is_upper_bound(i)) {
                    for (j = 0; j < dim; ++j) {
                        G_bar[j] += get_C(i) * Q_i[j];
                    }
                }
            }
        }
    }

//    std::cout << "\nG[" << active_size << "] = ";
//    for (std::vector<double>::const_iterator it = G.begin(); it != G.end(); ++it) {
//        std::cout << *it << " ";
//    }
//
//    std::cout << "\nG_bar[" << active_size << "] = ";
//    for (std::vector<double>::const_iterator it = G_bar.begin(); it != G_bar.end(); ++it) {
//        std::cout << *it << " ";
//    }
//    exit(0);

//    std::cout << "\nG_bar[" << active_size << "] = ";
//    for (std::vector<double>::const_iterator it = G_bar.begin(); it != G_bar.end(); ++it) {
//        std::cout << *it << " ";
//    }

    // optimization step

    int iter = 0;
    unsigned int counter = std::min(dim, 1000u) + 1;
    //std::cout << "\ncounter=" << counter << ", shrinking=" << shrinking << "\n\n";
//    bool calledDoShrinking = false;
//    int loopAfterDoShrinking = 0;
    while (true) {
        // show progress and do shrinking
        if (--counter == 0) {
            counter = std::min(dim, 1000u);
//            std::cout << " new counter = " << counter;
            if (shrinking) {
//                std::cout << "\nCALLING do_shrinking()\n";
//                std::cout << "\n\nG[" << active_size << "][" << dim << "]";
//                std::cout.flush();
                do_shrinking();
            }
            std::cout << ".";
//            std::cout.flush();
//            std::cout << "\nNew counter = " << counter;
//            calledDoShrinking = true;
        }

//        if (calledDoShrinking) {
//            ++loopAfterDoShrinking;
//            std::cout << "\n---[loop " << loopAfterDoShrinking << " after do_shrinking";
//            if (4 == loopAfterDoShrinking) {
//                std::cout << "\n\nG[" << active_size << "][" << dim << "] = ";
//                std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3);
//                //for (std::vector<double>::const_iterator it = G.begin(); it != G.end(); ++it) {
//                for (unsigned int ite=0; ite<active_size; ++ite) {
//                    //std::cout << *it << " ";
//                    std::cout << G[ite] << " ";
//                    if (ite % 20 == 19) std::cout << "\n";
//                }
//                exit(0);
//            }
//        }

        unsigned int i, j;
//        std::cout << "select_working_set(UNDEF, UNDEF)\n\n";
        if (select_working_set(i, j) != 0) {
//            std::cout << "\nselect_working_set(" << i << ", " << j << ") != 0 (a)";
            // reconstruct the whole gradient
            reconstruct_gradient();
            // reset active set size and check
//            unsigned int prev_active_size = active_size;
            active_size = dim;
            std::cout << "*";
            if (select_working_set(i, j) != 0) {
//                std::cout << "\nselect_working_set(" << i << ", " << j << ") != 0 -> break";

//                std::cout << "\n\nG[" << active_size << "][" << dim << "][" << prev_active_size << "] = ";
//                for (unsigned int ite=0; ite<active_size; ++ite) {
//                    if (ite == prev_active_size) std::cout << "\n==========[active_size]================\n";
//                    std::cout << G[ite] << " ";
//                    if (ite % 20 == 19) std::cout << "\n";
//                }
//                exit(0);

                break;
            } else {
                counter = 1; // do shrinking next iteration
            }
        }

        ++iter;

        // update alpha[i] and alpha[j], handle bounds carefully

//        std::cout << "\n\nupdate alpha[" << i << "] and alpha[" << j << "], handle bounds carefully (iter " << iter << ")\n\n";

        const std::vector<Qfloat> Q_i = Q->get_Q(i, active_size);  //some ELEMENTS have DIFFERENT signs!!!!
        const std::vector<Qfloat> Q_j = Q->get_Q(j, active_size);

        double C_i = get_C(i);
        double C_j = get_C(j);

//        if (calledDoShrinking) {
////            std::cout << "\n\ny[" << i << "] = ";
////            for (unsigned int q = 0; q < active_size; ++q) {
////                std::cout << static_cast<int>(y[q]) << " ";
////                if (q % 20 == 19) std::cout << "\n";
////            }
//            std::cout << "\n\nQ_ij[" << i << "][" << j << "] = ";
//            for (unsigned int q = 0; q < active_size; ++q) {
//                std::cout << Q_i[q] << "," << Q_j[q] << "   ";
//                if (q % 20 == 19) std::cout << "\n";
//            }
//            exit(0);
//        }

        //std::cout << "\nC_i = " << C_i << ", C_j = " << C_j;
        

        double old_alpha_i = alpha[i];
        double old_alpha_j = alpha[j];

        if (y[i] != y[j]) {
            double quad_coef = Q_i[i] + Q_j[j] + 2 * Q_i[j];
            if (quad_coef <= 0) {
                quad_coef = TAU;
            }
            double delta = (-G[i] - G[j]) / quad_coef;
            double diff = alpha[i] - alpha[j];
            alpha[i] += delta;
            alpha[j] += delta;

            if (diff > 0) {
                if (alpha[j] < 0) {
                    alpha[j] = 0;
                    alpha[i] = diff;
                }
            } else {
                if (alpha[i] < 0) {
                    alpha[i] = 0;
                    alpha[j] = -diff;
                }
            }
            if (diff > C_i - C_j) {
                if (alpha[i] > C_i) {
                    alpha[i] = C_i;
                    alpha[j] = C_i - diff;
                }
            } else {
                if (alpha[j] > C_j) {
                    alpha[j] = C_j;
                    alpha[i] = C_j + diff;
                }
            }
        } else {
            double quad_coef = Q_i[i] + Q_j[j] - 2 * Q_i[j];
            if (quad_coef <= 0) {
                quad_coef = TAU;
            }
            double delta = (G[i] - G[j]) / quad_coef;
            double sum = alpha[i] + alpha[j];
            alpha[i] -= delta;
            alpha[j] += delta;

            if (sum > C_i) {
                if (alpha[i] > C_i) {
                    alpha[i] = C_i;
                    alpha[j] = sum - C_i;
                }
            } else {
                if (alpha[j] < 0) {
                    alpha[j] = 0;
                    alpha[i] = sum;
                }
            }
            if (sum > C_j) {
                if (alpha[j] > C_j) {
                    alpha[j] = C_j;
                    alpha[i] = sum - C_j;
                }
            } else {
                if (alpha[i] < 0) {
                    alpha[i] = 0;
                    alpha[j] = sum;
                }
            }
        }

        // update G


        double delta_alpha_i = alpha[i] - old_alpha_i;
        double delta_alpha_j = alpha[j] - old_alpha_j;

        //std::cout << "\ndelta_alpha[i,j] = " << std::setiosflags(std::ios::fixed) << std::setprecision(3) << delta_alpha_i << ", " << delta_alpha_j;

        for (unsigned int k = 0; k < active_size; ++k) {
            G[k] += Q_i[k] * delta_alpha_i + Q_j[k] * delta_alpha_j;
        }

//        std::cout << "\n==[Updated G]==\n";

//        for (unsigned int t = 0; t < 10; t++) {
//            std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) << "\ndelta_alpha_i=" << delta_alpha_i << ", delta_alpha_j=" << delta_alpha_j << ", alpha[t]=" << alpha[t] << ", Q_i[t]=" << Q_i[t] << ", Q_j[" << j << "]=" << Q_j[t] << " ";
//        }

        //exit(0);

        // update alpha_status and G_bar

        {
            bool ui = is_upper_bound(i);
            bool uj = is_upper_bound(j);
            update_alpha_status(i);
            update_alpha_status(j);
            unsigned int k;
            if (ui != is_upper_bound(i)) {
                const std::vector<Qfloat> Q_i = Q->get_Q(i, dim);
                if (ui) {
                    for (k = 0; k < dim; ++k) {
                        G_bar[k] -= C_i * Q_i[k];
                    }
                } else {
                    for (k = 0; k < dim; ++k) {
                        G_bar[k] += C_i * Q_i[k];
                    }
                }
            }

            if (uj != is_upper_bound(j)) {
                const std::vector<Qfloat> Q_j = Q->get_Q(j, dim);
                if (uj) {
                    for (k = 0; k < dim; ++k) {
                        G_bar[k] -= C_j * Q_j[k];
                    }
                } else {
                    for (k = 0; k < dim; ++k) {
                        G_bar[k] += C_j * Q_j[k];
                    }
                }
            }
        }
    }

    // calculate rho

    si.rho = calculate_rho();

    //std::cout << "\nRho = " << si.rho;

    // calculate objective value
    {
        //std::cout << "\n\nQ_[i,j] = \n";
        double v = 0;
        for (unsigned int i = 0; i < dim; ++i) {
         //   std::cout << "(" << Q_i[i] << "," << Q_j[i] << ") ";
         //   if (i % 20 == 19) std::cout << "\n";
            v += alpha[i] * (G[i] + p[i]);
        }
//std::cout << "\n\n";
        si.obj = v / 2;
    }
    //std::cout << "\nOBJ = " << si.obj;

    // put back the solution
    {
        for (unsigned int i = 0; i < dim; ++i) {
            alpha_[active_set[i]] = alpha[i];
        }
    }

    // juggle everything back
    /*{
        for(int i=0;i<l;i++)
            while(active_set[i] != i)
                swap_index(i,active_set[i]);
                // or Q.swap_index(i,active_set[i]);
    }*/

    si.upper_bound_p = Cp;
    si.upper_bound_n = Cn;

    //std::cout << "\nUpper_Bound Cp, Cn = " << Cp << ", " << Cn;

    std::cout << "\noptimization finished, #iter = " << iter << "\n";

    si_ = si;

    p.clear();
    y.clear();
    alpha.clear();
    alpha_status.clear();
    active_set.clear();
    G.clear();
    G_bar.clear();
}

// return 1 if already optimal, return 0 otherwise

int SVM_Solver::select_working_set(unsigned int &out_i, unsigned int &out_j) {
    // return i,j such that
    // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
    // j: minimizes the decrease of obj value
    //    (if quadratic coefficeint <= 0, replace it with tau)
    //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
    
    double Gmax  = -INF;
    double Gmax2 = -INF;
    int Gmax_idx = -1;
    int Gmin_idx = -1;
    double obj_diff_min = INF;

    for (unsigned int t = 0; t < active_size; t++) {
//        std::cout << (int)y[t] << " ";
//        if (t % 20 == 19) std::cout << "\n";
        if (static_cast<int>(y[t]) == +1) {
            if (!is_upper_bound(t)) {
                if (-G[t] >= Gmax) {
                    Gmax = -G[t];
                    Gmax_idx = t;
                }
            }
        } else {
            if (!is_lower_bound(t)) {
                if (G[t] >= Gmax) {
                    Gmax = G[t];
                    Gmax_idx = t;
                }
            }
        }
    }
//    exit(0);

//    for (unsigned int t = 0; t < active_size; t++) {
//        std::cout << G[t] << " ";
//        if (t % 20 == 19) std::cout << "\n";
//    }
//    std::cout << "\nGmax_idx=" << Gmax_idx << ", Gmax=" << Gmax;
//    exit(0);

    int i = Gmax_idx;
    //const std::vector<Qfloat> Q_i;
    //if (i != -1) { // NULL Q_i not accessed: Gmax=-INF if i=-1
        const std::vector<Qfloat> Q_i = Q->get_Q((i > 0 ? i : 0), active_size);
    //}

//        std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3);
        
//        int counter = 0;
//    std::cout << "\nCalculating j. Gmax2 = ";
    for (unsigned int j = 0; j < active_size; ++j) {
        if (static_cast<int>(y[j]) == +1) {
            if (!is_lower_bound(j)) {
                //std::cout << "O ";
                double grad_diff = Gmax + G[j];
                if (G[j] >= Gmax2) {
                    Gmax2 = G[j];
                    //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3) << "(" << j << ",+)" << Gmax2 << " ";
                }
                if (grad_diff > 0 && Gmax_idx > -1) {
                    double obj_diff;
                    double quad_coef = Q_i[i] + QD[j] - 2.0 * static_cast<int>(y[i]) * Q_i[j];
                    if (quad_coef > 0) {
                        obj_diff = -(grad_diff * grad_diff) / quad_coef;
                    } else {
                        obj_diff = -(grad_diff * grad_diff) / TAU;
                    }
                    if (obj_diff <= obj_diff_min) {
//                        std::cout << "\nGmin_idx[a]=" << j << " obj_diff=" << obj_diff;
                        Gmin_idx = j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        } else {
            if (!is_upper_bound(j)) {
                double grad_diff = Gmax - G[j];
                if (-G[j] >= Gmax2) {
                    //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3) << "(Gmax=" << Gmax << ", -G[" << j << "]=Gmax2=" << (-G[j]) << ", grad=" << grad_diff << ", Gmax2_was=" << Gmax2 << ", -) ";
                    Gmax2 = -G[j];
                }
                if (grad_diff > 0 && Gmax_idx > -1) {
                    double obj_diff;
                    double quad_coef = Q_i[i] + QD[j] + 2.0 * static_cast<int>(y[i]) * Q_i[j];
//                    if (counter < 40) {
//                        ++counter;
//                    std::cout << "\n" << std::setiosflags(std::ios::fixed) << std::setprecision(2) << quad_coef << " " <<
//                        Q_i[i] << " " <<
//                        QD[j] << " " <<
//                        (int)y[i] << " " <<
//                        Q_i[j] << " ";
//                    }
                    if (quad_coef > 0) {
                        obj_diff = -(grad_diff * grad_diff) / quad_coef;
                    } else {
                        obj_diff = -(grad_diff * grad_diff) / TAU;
                    }
                    //std::cout <<  std::setiosflags(std::ios::fixed) << std::setprecision(2) << obj_diff << " ";
                    if (obj_diff <= obj_diff_min) {
//                        std::cout << "\nGmin_idx[b]=" << j << " obj_diff=" << obj_diff;
                        Gmin_idx = j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
    }


//    std::cout  << std::setiosflags(std::ios::fixed) << std::setprecision(3) << "\n[2] Gmax_idx=" << Gmax_idx << ", Gmin_idx=" << Gmin_idx << ", Gmax=" << Gmax << ", Gmax2=" << Gmax2 << ", (Gmax + Gmax2)=" << (Gmax + Gmax2) << ", eps=" << eps;

    if (Gmax + Gmax2 < eps) {
        return 1;
    }

    out_i = Gmax_idx;
    out_j = Gmin_idx;
    return 0;
}

bool SVM_Solver::be_shrunk(int i, double Gmax1, double Gmax2) {
    if (is_upper_bound(i)) {
        if (static_cast<int>(y[i]) == +1) {
            return (-G[i] > Gmax1);
        }
        return (-G[i] > Gmax2);
    }
    if (is_lower_bound(i)) {
        if (static_cast<int>(y[i]) == +1) {
            return (G[i] > Gmax2);
        }
        return (G[i] > Gmax1);
    }
    return false;
}

void SVM_Solver::do_shrinking() {
    unsigned int i;
    double Gmax1 = -INF; // max { -y_i * grad(f)_i | i in I_up(\alpha) }
    double Gmax2 = -INF; // max { y_i * grad(f)_i | i in I_low(\alpha) }
    //std::cout << "\n\ndo_shrinking(START) active_size=" << active_size << "\n\n";
    // find maximal violating pair first

    for (i = 0; i < active_size; ++i) {
        //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(8) << (double)G[i] << " ";
        if (static_cast<int>(y[i]) == +1) {
            if (!is_upper_bound(i)) {
                if (-G[i] >= Gmax1) {
                    Gmax1 = -G[i];
                }
            }
            if (!is_lower_bound(i)) {
                if (G[i] >= Gmax2) {
                    Gmax2 = G[i];
                }
            }
        } else {
            if (!is_upper_bound(i)) {
                if (-G[i] >= Gmax2) {
                    Gmax2 = -G[i];
                }
            }
            if (!is_lower_bound(i)) {
                if (G[i] >= Gmax1) {
                    Gmax1 = G[i];
                }
            }
        }
    }

    //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6) << "\n\nunshrink=" << (int)unshrink << ", (Gmax1 + Gmax2 <= eps * 10) = ("<< (double)Gmax1 << " + " << (double)Gmax2 << " <= " << (eps * 10) << ")\n\n";

    if (unshrink == false && Gmax1 + Gmax2 <= eps * 10) {
        unshrink = true;
        reconstruct_gradient();
        active_size = dimension;
        std::cout << "*";
    }

    //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3) << "\n[3] Gmax1=" << Gmax1 << ", Gmax2=" << Gmax2;
    
    for (i = 0; i < active_size; ++i) {
        //std::cout << "\nbe_shrunk(" << i << ") = ";
        if (be_shrunk(i, Gmax1, Gmax2)) {
            //std::cout << "true\n";
            //std::cout << "\nbe_shrunk(" << i << ") = true";
            --active_size;
            while (active_size > i) {
                if (!be_shrunk(active_size, Gmax1, Gmax2)) {
                    //std::cout << "\n!be_shrunk(" << active_size << ", " << Gmax1 << ", " << Gmax2 << ") -> swap(" << i << ", " << active_size << ")";
                    swap_index(i, active_size);
                    break;
                }
                --active_size;
            }
        } else {
            //std::cout << "false\n";
        }
    }
    //std::cout << "\n\ndo_shrinking(END) active_size = " << active_size << "\n\n";
}

double SVM_Solver::calculate_rho() {
    double r;
    int nr_free = 0;
    double ub = INF, lb = -INF, sum_free = 0;
    for (unsigned int i = 0; i < active_size; ++i) {
        double yG = static_cast<int>(y[i]) * G[i];

        if (is_upper_bound(i)) {
            if (static_cast<int>(y[i]) == -1) {
                ub = std::min(ub, yG);
            } else {
                lb = std::max(lb, yG);
            }
        } else if (is_lower_bound(i)) {
            if (static_cast<int>(y[i]) == +1) {
                ub = std::min(ub, yG);
            } else {
                lb = std::max(lb, yG);
            }
        } else {
            ++nr_free;
            sum_free += yG;
        }
    }

    if (nr_free > 0) {
        r = sum_free / nr_free;
    } else {
        r = (ub + lb) / 2;
    }
    return r;
}
