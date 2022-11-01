// [[Rcpp::depends(RcppArmadillo)]]
// Ami Falk, 9/8/22

#include "RcppArmadillo.h"
using namespace arma;

void case_1(int k, mat &R, mat &Y, mat &c_hat)
{

    mat A_mat = (R + R(k, k) * eye(size(Y)));

    vec summation(Y.n_rows);

    for (unsigned int j = k + 1; j < Y.n_rows; j++)
    {
        summation += (R(k, j) * Y.col(j));
    }

    vec B_mat = c_hat.col(k) - summation;

    Y.col(k) = solve(A_mat, B_mat);
}

void case_2(int k, mat &R, mat &Y, mat &c_hat)
{

    mat small_R = {{R(k - 1, k - 1), R(k - 1, k)},
                   {R(k, k - 1), R(k, k)}};

    mat A_mat = kron(eye(2, 2), R) + kron(small_R, eye(size(Y)));

    mat summation(Y.n_rows, 2);

    for (unsigned int j = k + 1; j < Y.n_rows; j++)
    {
        summation.col(0) += (R(k - 1, j) * Y.col(j));
        summation.col(1) += (R(k, j) * Y.col(j));
    }

    vec B_vec = vectorise(c_hat.cols(k - 1, k) - summation);

    mat solution = solve(A_mat, B_vec);
    solution.reshape(Y.n_rows, 2);

    Y.cols(k - 1, k) = solution;
}

// [[Rcpp::export]]
mat clyap_internal(mat A, mat C)
{

    mat U; // unitary matrix
    mat R; // upper triangular matrix
    schur(U, R, A);

    mat Y(size(A));
    mat c_hat = U.t() * (C * U);

    for (unsigned int k = Y.n_rows - 1; k >= 0; k--)
    {

        if (k == 0)
        {
            case_1(k, R, Y, c_hat);
            break;
        }

        if (R(k, k - 1) == 0)
        {
            case_1(k, R, Y, c_hat);
            continue;
        }

        if (R(k, k - 1) != 0)
        {
            case_2(k, R, Y, c_hat);
        }
    }

    return U * Y * (U.t());
}