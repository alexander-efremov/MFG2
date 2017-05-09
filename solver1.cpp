#include <cassert>
#include <cstdio>
#include <cstring>
#include "inc/malgo.h"
#include "inc/utils.h"

/**
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы a (нумеруется: [0;n-1])
	 * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
inline void fill_b(double *b, int n, double tau, double sigma_sq) {
    b[0] = 0.;
    for (int i = 1; i < n; ++i)
        b[i] = -tau * sigma_sq;
}

inline void fill_c(double *c, int n, double tau, double sigma_sq, double h_sq) {
    for (int i = 0; i < n; ++i)
        c[i] = h_sq + 2. * tau * sigma_sq;
}

inline void fill_d(double *d, int n, double tau, double sigma_sq) {
    for (int i = 0; i < n - 1; ++i)
        d[i] = -tau * sigma_sq;
    d[n - 1] = 0.;
}

inline double func_u(double a_coef, double t, double x) {
    return a_coef * t * x * (1. - x);
}

inline double get_f(double sigma_sq, double a_coef, double x, double t) {
    return a_coef * x * (1. - x) + func_u(a_coef, t, x) * a_coef * t * (1. - 2. * x) + 2. * sigma_sq * a_coef * t;
}

double analytical_solution_1(double a_coef, double t, double x) {
    return a_coef * t * x * (1. - x);
}

double
get_rp_value(int ii, const double *m_pr, double time, double a, double b, double h, double a_coef,
             double tau) {
    // определяем точку из которой будем опускать траекторию
    double x = a + ii * h;

    // опускаем траекторию из этой точки
    double alpha = func_u(a_coef, time, x);
    x -= tau * alpha;

    if (x < a || x > b)
        printf("Time value %.8le! ERROR INDEX i=%d : x=%.8le\n ", time, ii, x);

    auto i = (int) ((x - a) / h);
    double xi_l = a + i * h; // левая граница интервала
    double xi_r = a + (i + 1) * h; // правая граница интервала
    assert(xi_l <= x && x <= xi_r);

    // считаем интеграл
    // вычислим значение функции в нашей точке по формуле 9
    double m = m_pr[i] * (xi_r - x) / h + m_pr[i + 1] * (x - xi_l) / h;

    return m;
}

void fill_rp(double *rp, double *m_pr, double time, int n_1, double tau, double h, double h_sq, double sigma_sq,
             double a_coef, double a, double b) {
    for (int i = 1; i < n_1 - 1; ++i) {
        double f = get_f(sigma_sq, a_coef, i * h, time);
        double d = get_rp_value(i, m_pr, time, a, b, h, a_coef, tau);
        rp[i] = (tau * f + d) * h_sq;
    }
    rp[0] = 0.;
    rp[n_1 - 1] = 0.;
}

void assert_params(double h, double h_sq, double h_2, double sigma_sq, int n_1, double sigma,
                   double a, double tau, double a_coef, double time_step_cnt) {
    assert(h > 0.);
    assert(h_sq == h * h);
    assert(h_2 == h / 2.);
    assert(sigma_sq == sigma * sigma);
    assert(n_1 > 0);
    assert(a == 0.);
    assert(tau > 0.);
    assert(a_coef > 0.);
    assert(time_step_cnt >= 1);
    // (3.19)
//    printf("h * h = %e\n", h * h);
//    printf("8 * tau * sigma_sq = %e\n", 8. * tau * sigma_sq);
//    fflush(stdout);
    //assert(h_sq <= 8 * tau * sigma_sq);
    if (h_sq > 8. * tau * sigma_sq)
        printf("\n!!!!!!!!!!!!!!!!!! H*H<=8*tau*sigma_sq FAILED\n");
}

/**
    * n - число уравнений (строк матрицы)
    * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
    * c - главная диагональ матрицы a (нумеруется: [0;n-1])
    * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
    * f - правая часть (столбец)
    * x - решение, массив x будет содержать ответ
    */
void print_thomas_arrays(double *b, double *c, double *d, int n) {
    printf("b\n");
    for (int i = 0; i < n; ++i) {
        printf("%e ", b[i]);
    }
    printf("\nc\n");
    for (int i = 0; i < n; ++i) {
        printf("%e ", c[i]);
    }
    printf("\nd\n");
    for (int i = 0; i < n; ++i) {
        printf("%e ", d[i]);
    }
    printf("\n");
}

void fill_arr_by_ex_sol(double *arr, int n_1, double time, double a_coef, double h, double a) {
    for (int i = 0; i < n_1; ++i)
        arr[i] = analytical_solution_1(a_coef, time, a + i * h);
}

void solve_1(int n_1, double h, double h_sq, double h_2, double sigma_sq, double sigma, double a, double b,
             double a_coef, double tau, int time_step_cnt) {
    assert_params(h, h_sq, h_2, sigma_sq, n_1, sigma, a, tau, a_coef, time_step_cnt);

    double *m = (double *) malloc(n_1 * sizeof(double));
    double *m_pr = (double *) malloc(n_1 * sizeof(double));
    double *ex_m = (double *) malloc(n_1 * sizeof(double));
    double *rp = (double *) malloc(n_1 * sizeof(double));
    double *err = (double *) malloc(n_1 * sizeof(double));
    double *b_t = (double *) malloc(n_1 * sizeof(double));
    double *c_t = (double *) malloc(n_1 * sizeof(double));
    double *d_t = (double *) malloc(n_1 * sizeof(double));

    fill_b(b_t, n_1, tau, sigma_sq);
    fill_c(c_t, n_1, tau, sigma_sq, h_sq);
    fill_d(d_t, n_1, tau, sigma_sq);
    //print_thomas_arrays(b_t, c_t, d_t, n_1);

    for (int i = 0; i < n_1; ++i) m[i] = m_pr[i] = rp[i] = 0.;

    fill_arr_by_ex_sol(m_pr, n_1, 0., a_coef, h, a);

    //print_vector("M_PR", m_pr, n_1);

    for (int tl = 1; tl <= time_step_cnt; ++tl) {
        fill_rp(rp, m_pr, tau * tl, n_1, tau, h, h_sq, sigma_sq, a_coef, a, b);
        thomas_algo_verzh_modified(n_1, b_t, c_t, d_t, rp, m);
        memcpy(m_pr, m, n_1 * sizeof(double));
    }

    fill_arr_by_ex_sol(ex_m, n_1, tau * time_step_cnt, a_coef, h, a);
    fill_arr_diff(err, ex_m, m, n_1);
    print_XY("exact_numer", n_1, h, time_step_cnt, a, b, tau, ex_m, m);
    printf("uniform norm of l1 norm of error %22.14le\n", get_l1_norm(n_1, err));

    free(ex_m);
    free(m_pr);
    free(rp);
    free(err);
    free(b_t);
    free(c_t);
    free(d_t);
}