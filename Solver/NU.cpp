/* 
 * File:   NU.cpp
 * Author: Lorenzo Alberton
 * 
 * Created on July 5, 2010, 2:00 PM
 */

#include "NU.h"

/**
 *  return 1 if already optimal, return 0 otherwise
 */
int SVM_Solver_NU::select_working_set(int &out_i, int &out_j) {
    // return i,j such that y_i = y_j and
    // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
    // j: minimizes the decrease of obj value
    //    (if quadratic coefficeint <= 0, replace it with tau)
    //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)

    double Gmaxp  = -INF;
    double Gmaxp2 = -INF;
    int Gmaxp_idx = -1;

    double Gmaxn  = -INF;
    double Gmaxn2 = -INF;
    int Gmaxn_idx = -1;

    int Gmin_idx = -1;
    double obj_diff_min = INF;

    for (unsigned int t = 0; t < active_size; t++) {
        if (y[t] == +1) {
            if (!is_upper_bound(t)) {
                if (-G[t] >= Gmaxp) {
                    Gmaxp = -G[t];
                    Gmaxp_idx = t;
                }
            }
        } else {
            if (!is_lower_bound(t)) {
                if (G[t] >= Gmaxn) {
                    Gmaxn = G[t];
                    Gmaxn_idx = t;
                }
            }
        }
    }

    int idxpos = Gmaxp_idx;
    int idxneg = Gmaxn_idx;
    //const std::vector<Qfloat> Q_idxpos;
    //const std::vector<Qfloat> Q_idxneg;
    //if (idxpos != -1) {
        // NULL Q_ip not accessed: Gmaxp=-INF if ip=-1
        const std::vector<Qfloat> Q_idxpos = Q->get_Q(idxpos, active_size);
    //}
    //if (idxneg != -1) {
        const std::vector<Qfloat> Q_idxneg = Q->get_Q(idxneg, active_size);
    //}

    for (unsigned int j = 0; j < active_size; j++) {
        if (y[j] == +1) {
            if (!is_lower_bound(j)) {
                double grad_diff = Gmaxp + G[j];
                if (G[j] >= Gmaxp2) {
                    Gmaxp2 = G[j];
                }
                if (grad_diff > 0) {
                    double obj_diff;
                    double quad_coef = Q_idxpos[idxpos] + QD[j] - 2 * Q_idxpos[j];
                    if (quad_coef > 0) {
                        obj_diff = -(grad_diff * grad_diff) / quad_coef;
                    } else {
                        obj_diff = -(grad_diff * grad_diff) / TAU;
                    }
                    if (obj_diff <= obj_diff_min) {
                        Gmin_idx = j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        } else {
            if (!is_upper_bound(j)) {
                double grad_diff = Gmaxn - G[j];
                if (-G[j] >= Gmaxn2) {
                    Gmaxn2 = -G[j];
                }
                if (grad_diff > 0) {
                    double obj_diff;
                    double quad_coef = Q_idxneg[idxneg] + QD[j] - 2 * Q_idxneg[j];
                    if (quad_coef > 0) {
                        obj_diff = -(grad_diff * grad_diff) / quad_coef;
                    } else {
                        obj_diff = -(grad_diff * grad_diff) / TAU;
                    }
                    if (obj_diff <= obj_diff_min) {
                        Gmin_idx = j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
    }

    if (max(Gmaxp + Gmaxp2, Gmaxn + Gmaxn2) < eps) {
        return 1;
    }

    if (y[Gmin_idx] == +1) {
        out_i = Gmaxp_idx;
    } else {
        out_i = Gmaxn_idx;
    }
    out_j = Gmin_idx;

    return 0;
}

/**
 * 
 * @param i
 * @param Gmax1
 * @param Gmax2
 * @param Gmax3
 * @param Gmax4
 * @return
 */
bool SVM_Solver_NU::be_shrunk(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4) {
    if (is_upper_bound(i)) {
        if (y[i] == +1) {
            return (-G[i] > Gmax1);
        }
        return (-G[i] > Gmax4);
    }
    if (is_lower_bound(i)) {
        if (y[i] == +1) {
            return (G[i] > Gmax2);
        }
        return (G[i] > Gmax3);
    }
    return false;
}

void SVM_Solver_NU::do_shrinking() {
    double Gmax1 = -INF; // max { -y_i * grad(f)_i | y_i = +1, i in I_up(\alpha) }
    double Gmax2 = -INF; // max { y_i * grad(f)_i | y_i = +1, i in I_low(\alpha) }
    double Gmax3 = -INF; // max { -y_i * grad(f)_i | y_i = -1, i in I_up(\alpha) }
    double Gmax4 = -INF; // max { y_i * grad(f)_i | y_i = -1, i in I_low(\alpha) }

    // find maximal violating pair first
    unsigned int i;
    for (i = 0; i < active_size; i++) {
        if (!is_upper_bound(i)) {
            if (y[i] == +1) {
                if (-G[i] > Gmax1) {
                    Gmax1 = -G[i];
                }
            } else if (-G[i] > Gmax4) {
                Gmax4 = -G[i];
            }
        }
        if (!is_lower_bound(i)) {
            if (y[i] == +1) {
                if (G[i] > Gmax2) {
                    Gmax2 = G[i];
                }
            } else if (G[i] > Gmax3) {
                Gmax3 = G[i];
            }
        }
    }

    if (unshrink == false && max(Gmax1 + Gmax2, Gmax3 + Gmax4) <= eps * 10) {
        unshrink = true;
        reconstruct_gradient();
        active_size = dimension;
    }

    for (i = 0; i < active_size; i++)
        if (be_shrunk(i, Gmax1, Gmax2, Gmax3, Gmax4)) {
            active_size--;
            while (active_size > i) {
                if (!be_shrunk(active_size, Gmax1, Gmax2, Gmax3, Gmax4)) {
                    swap_index(i, active_size);
                    break;
                }
                --active_size;
            }
        }
}

double SVM_Solver_NU::calculate_rho() {
    int nr_free1     = 0, nr_free2  = 0;
    double sum_free1 = 0, sum_free2 = 0;
    double ub1 = INF,  ub2 = INF;
    double lb1 = -INF, lb2 = -INF;
    
    for (unsigned int i = 0; i < active_size; i++) {
        if (y[i] == +1) {
            if (is_upper_bound(i)) {
                lb1 = max(lb1, G[i]);
            } else if (is_lower_bound(i)) {
                ub1 = min(ub1, G[i]);
            } else {
                ++nr_free1;
                sum_free1 += G[i];
            }
        } else {
            if (is_upper_bound(i)) {
                lb2 = max(lb2, G[i]);
            } else if (is_lower_bound(i)) {
                ub2 = min(ub2, G[i]);
            } else {
                ++nr_free2;
                sum_free2 += G[i];
            }
        }
    }

    double r1, r2;
    if (nr_free1 > 0) {
        r1 = sum_free1 / nr_free1;
    } else {
        r1 = (ub1 + lb1) / 2;
    }

    if (nr_free2 > 0) {
        r2 = sum_free2 / nr_free2;
    } else {
        r2 = (ub2 + lb2) / 2;
    }
    si.r = (r1 + r2) / 2;
    return (r1 - r2) / 2;
}
