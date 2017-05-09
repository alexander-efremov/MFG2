#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "solver1.h"


void print_params(int n, int n_1, double a_coef, double sigma, double sigma_sq, double h, double h_sq, double tau,
                  int time_step_cnt) {
    printf("\nN = %d\n", n);
    //printf("N_1 = %d\n", n_1);
    //printf("A_COEF = %le\n", a_coef);
    //printf("SIGMA = %le\n", sigma);
    //printf("SIGMA_SQ = %le\n", sigma_sq);
    printf("H = %le\n", h);
    //printf("H_SQ = %le\n", h_sq);
    printf("TAU = %le\n", tau);
    printf("TIME_STEP_CNT = %d\n", time_step_cnt);
//    printf("%f ", A - H_2);
//    for (int i = 1; i < N_1; ++i) {
//        printf("%f ", i * H);
//    }
//    printf("%f\n", B + H_2);
    fflush(stdout);
}


void run_solver_1(int _n, double _tau, int _time_step_cnt, double _sigma) {
    double a = 0.;
    double b = 1.;
    double a_coef = 1.;
    int time_step_cnt = _time_step_cnt;
    double tau = _tau;
    double sigma = _sigma;
    double sigma_sq = sigma * sigma;
    int n_1 = _n + 1;
    double h = (b - a) / _n;
    double h_2 = 0.5 * h;
    double h_sq = h * h;
    print_params(n_1 - 1, n_1, a_coef, sigma, sigma_sq, h, h_sq, tau, time_step_cnt);
    solve_1(n_1, h, h_sq, h_2, sigma_sq, sigma, a, b, a_coef, tau, time_step_cnt);
}

void run_solver_1() {
    int time_step_cnt = 1;
    double tau = 1e-3;
    double sigma = 1.;
    int n = 20;
    run_solver_1(n, tau, time_step_cnt, sigma);
}

TEST_CASE("mfg2_solver_1", "[run_solver_1]") {
    run_solver_1();
}

TEST_CASE("mfg2_solver_1_many", "[run_solver_1_l1_diag]") {
    int exp_cnt = 5;
    for (int i = 0; i < exp_cnt; ++i) {
        printf("\n\n=============== EXPERIMENT %d ================ \n\n", i + 1);
        int nx = 0;
        double tau = 1e-3;
        double sigma = 1.;
        int tsc;
        switch (i) {
            case 0:
                nx = 50;
                tau = tau / 1.;
                tsc = 10 * (int) 1.;
                break;
            case 1:
                nx = 100;
                tau = tau / 4.;
                tsc = 10 * (int) 4;
                break;
            case 2:
                nx = 200;
                tau = tau / 16.;
                tsc = 10 * 16;
                break;
            case 3:
                nx = 400;
                tau = tau / 64.;
                tsc = 10 * 64;
                break;
            case 4:
                nx = 800;
                tau = tau / 256.;
                tsc = 10 * 256;
                break;
            case 5:
                nx = 1600;
                tau = tau / 512.;
                tsc = 10 * 512;
                break;
            default:
                return;
        }

        assert(tau * tsc == 1e-2);
        //printf("TAU * TIME_STEP_COUNT = %e", tau * tsc);
        run_solver_1(nx, tau, tsc, sigma);
    }
}

TEST_CASE("mfg2_solver_1_tsc_const", "[run_solver_1_many_tcs_const]") {
    int exp_cnt = 5;
    int col_cnt = 5;
    for (int j = 0; j < col_cnt; ++j) {
        for (int i = 0; i < exp_cnt; ++i) {
            printf("\n\n=============== COL %d EXPERIMENT %d ================ \n\n", j + 1, i + 1);
            int nx = 0;
            double tau_d = pow(4., j);
            double tau = 1e-3 / tau_d;
            double sigma = 1.;
            int tcs_m = (int) pow(4., j);
            int time_step_cnt = 10 * tcs_m;
            switch (i) {
                case 0:
                    nx = 50;
                    break;
                case 1:
                    nx = 100;
                    break;
                case 2:
                    nx = 200;
                    break;
                case 3:
                    nx = 400;
                    break;
                case 4:
                    nx = 800;
                    break;
                case 5:
                    nx = 1600;
                    break;
                default:
                    return;
            }
            printf("TAU * TIME_STEP_COUNT = %e", tau * time_step_cnt);
            run_solver_1(nx, tau, time_step_cnt, sigma);
        }
    }
}