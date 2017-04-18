#include <malgo.h>
#include <utils.h>
#include "catch.hpp"


TEST_CASE("thomas_algo_test3") {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;

    //под
    a[0] = 0;
    a[1] = 5;
    a[2] = 1;

    //главная
    b[0] = 2;
    b[1] = 4;
    b[2] = -3;

    //над
    c[0] = -1;
    c[1] = 2;
    c[2] = 0;

    rp[0] = 3;
    rp[1] = 6;
    rp[2] = 2;

    printf("\n");
    print_matrix(rp, 1, n);

    thomas_algo(rp, n, a, b, c);

    print_matrix(rp, 1, n);

    REQUIRE(rp[0] == Approx(1.49).epsilon(1e-2));
    REQUIRE(rp[1] == Approx(-0.02).epsilon(1e-2));
    REQUIRE(rp[2] == Approx(-0.68).epsilon(1e-2));

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_CASE("thomas_algo_test4") {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i) {
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;
    }

    //под
    a[0] = 0;
    a[1] = 1;
    a[2] = 1;

    //главная
    b[0] = 1;
    b[1] = 2;
    b[2] = 3;

    //над
    c[0] = 1;
    c[1] = 1;
    c[2] = 0;

    rp[0] = 2;
    rp[1] = 4;
    rp[2] = 4;

    printf("\n");
    print_matrix(rp, 1, n);

    thomas_algo(rp, n, a, b, c);

    print_matrix(rp, 1, n);

    REQUIRE(rp[0] == Approx(1.));
    REQUIRE(rp[1] == Approx(1.));
    REQUIRE(rp[2] == Approx(1.));

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_CASE("thomas_algo_test5") {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i) {
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;
    }

    //под
    a[0] = 0;
    a[1] = 5;
    a[2] = 1;

    //главная
    c[0] = 2;
    c[1] = 4;
    c[2] = -3;

    //над
    b[0] = -1;
    b[1] = 2;
    b[2] = 0;

    rp[0] = 3;
    rp[1] = 6;
    rp[2] = 2;

    printf("\n");
    print_matrix(rp, 1, n);
/**
	 * n - число уравнений (строк матрицы)
	 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
    thomas_algo(n, a, c, b, rp, m);

    print_matrix(m, 1, n);

    REQUIRE(m[0] == Approx(1.49).epsilon(1e-2));
    REQUIRE(m[1] == Approx(-0.02).epsilon(1e-2));
    REQUIRE(m[2] == Approx(-0.68).epsilon(1e-2));

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_CASE("thomas_algo_test6") {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;

    //под
    a[0] = 0;
    a[1] = 1;
    a[2] = 1;

    //главная
    c[0] = 1;
    c[1] = 2;
    c[2] = 3;

    //над
    b[0] = 1;
    b[1] = 1;
    b[2] = 0;

    rp[0] = 2;
    rp[1] = 4;
    rp[2] = 4;

    printf("\n");
    print_matrix(rp, 1, n);

    thomas_algo(n, a, c, b, rp, m);

    print_matrix(m, 1, n);

    REQUIRE(m[0] == Approx(1.));
    REQUIRE(m[1] == Approx(1.));
    REQUIRE(m[2] == Approx(1.));

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_CASE("thomas_algo1_test1") {
    const int n = 5;

    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *d = (double *) malloc(n * sizeof(double));
    double *r = (double *) malloc(n * sizeof(double));
    double *x = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        b[i] = b[i] = c[i] = d[i] = x[i] = 0.;

    //под
    b[0] = 0.;
    b[1] = 2.;
    b[2] = 4.;
    b[3] = 4.;
    b[4] = 2.;

    //главная
    c[0] = 2.;
    c[1] = 9.;
    c[2] = 17.;
    c[3] = 15.;
    c[4] = 3.;

    //над
    d[0] = 1.;
    d[1] = 2.;
    d[2] = -4.;
    d[3] = -8.;
    d[4] = 0.;

    r[0] = -10.;
    r[1] = -26.;
    r[2] = -16.;
    r[3] = -2.;
    r[4] = 16.;

    printf("\n");
    print_matrix1(b, 1, n);
    printf("\n");
    print_matrix1(c, 1, n);
    printf("\n");
    print_matrix1(d, 1, n);
    printf("\n");
    print_matrix1(r, 1, n);

    thomas_algo_verzh(n, b, c, d, r, x);

    printf("\n");
    print_matrix(x, 1, n);

    REQUIRE(x[0] == Approx(-4.));
    REQUIRE(x[1] == Approx(-2.));
    REQUIRE(x[2] == Approx(0.));
    REQUIRE(x[3] == Approx(2.));
    REQUIRE(x[4] == Approx(4.));
    free(r);
    free(b);
    free(c);
    free(d);
    free(x);
}

TEST_CASE("thomas_algo1_test2") {
    const int n = 3;

    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *d = (double *) malloc(n * sizeof(double));
    double *r = (double *) malloc(n * sizeof(double));
    double *x = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        b[i] = b[i] = c[i] = d[i] = x[i] = 0.;

    //под
    b[0] = 0.;
    b[1] = 1.;
    b[2] = 1.;

    //главная
    c[0] = 1.;
    c[1] = 2.;
    c[2] = 3.;

    //над
    d[0] = 1.;
    d[1] = 1.;
    d[2] = 0.;

    r[0] = 2.;
    r[1] = 4.;
    r[2] = 4.;

    printf("\n");
    print_matrix1(b, 1, n);
    printf("\n");
    print_matrix1(c, 1, n);
    printf("\n");
    print_matrix1(d, 1, n);
    printf("\n");
    print_matrix1(r, 1, n);

    thomas_algo_verzh(n, b, c, d, r, x);

    printf("\n");
    print_matrix(x, 1, n);

    REQUIRE(x[0] == Approx(1.));
    REQUIRE(x[1] == Approx(1.));
    REQUIRE(x[2] == Approx(1.));

    free(r);
    free(b);
    free(c);
    free(d);
    free(x);
}

TEST_CASE("thomas_algo1_test3") {
    const int n = 7;

    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *d = (double *) malloc(n * sizeof(double));
    double *r = (double *) malloc(n * sizeof(double));
    double *x = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        b[i] = b[i] = c[i] = d[i] = x[i] = 0.;

    //под
    b[0] = 0.;
    b[1] = 0.;
    b[2] = 2.;
    b[3] = 4.;
    b[4] = 4.;
    b[5] = 2.;
    b[6] = 0.;

    //главная
    c[0] = 0.;
    c[1] = 2.;
    c[2] = 9.;
    c[3] = 17.;
    c[4] = 15.;
    c[5] = 3.;
    c[6] = 0.;

    //над
    d[0] = 0.;
    d[1] = 1.;
    d[2] = 2.;
    d[3] = -4.;
    d[4] = -8.;
    d[5] = 0.;
    d[6] = 0.;

    r[0] = 0;
    r[1] = -10.;
    r[2] = -26.;
    r[3] = -16.;
    r[4] = -2.;
    r[5] = 16.;
    r[6] = 0.;

    printf("\n");
    print_matrix1(b, 1, n);
    printf("\n");
    print_matrix1(c, 1, n);
    printf("\n");
    print_matrix1(d, 1, n);
    printf("\n");
    print_matrix1(r, 1, n);

    thomas_algo_verzh_modified(n, b, c, d, r, x);

    printf("\n");
    print_matrix(x, 1, n);

    REQUIRE(x[0] == Approx(0.));
    REQUIRE(x[1] == Approx(-4.));
    REQUIRE(x[2] == Approx(-2.));
    REQUIRE(x[3] == Approx(0.));
    REQUIRE(x[4] == Approx(2.));
    REQUIRE(x[5] == Approx(4.));
    REQUIRE(x[6] == Approx(0.));

    free(r);
    free(b);
    free(c);
    free(d);
    free(x);
}

TEST_CASE("thomas_algo1_test4") {
    const int n = 5;

    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *d = (double *) malloc(n * sizeof(double));
    double *r = (double *) malloc(n * sizeof(double));
    double *x = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        b[i] = b[i] = c[i] = d[i] = x[i] = 0.;

    //под
    b[0] = 0.;
    b[1] = 0.;
    b[2] = 1.;
    b[3] = 1.;
    b[4] = 0.;

    //главная
    c[0] = 0.;
    c[1] = 1.;
    c[2] = 2.;
    c[3] = 3.;
    c[4] = 0.;

    //над
    d[0] = 0.;
    d[1] = 1.;
    d[2] = 1.;
    d[3] = 0.;
    d[4] = 0.;

    //
    r[0] = 0.;
    r[1] = 2.;
    r[2] = 4.;
    r[3] = 4.;
    r[4] = 0.;

    printf("\n");
    print_matrix1(b, 1, n);
    printf("\n");
    print_matrix1(c, 1, n);
    printf("\n");
    print_matrix1(d, 1, n);
    printf("\n");
    print_matrix1(r, 1, n);

    thomas_algo_verzh_modified(n, b, c, d, r, x);

    printf("\n");
    print_matrix(x, 1, n);

    REQUIRE(x[0] == Approx(0.));
    REQUIRE(x[1] == Approx(1.));
    REQUIRE(x[2] == Approx(1.));
    REQUIRE(x[3] == Approx(1.));
    REQUIRE(x[4] == Approx(0.));

    free(r);
    free(b);
    free(c);
    free(d);
    free(x);
}