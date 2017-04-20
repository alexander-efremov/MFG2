#ifndef SOLVER1_H
#define SOLVER1_H

double *solve_1(int n, int n_1, double h, double h_sq, double h_2, double sigma_sq, double sigma, double a, double b,
                double a_coef, double tau, int time_step_cnt, double *exact_sol_to_fill);


#endif //SOLVER1_H