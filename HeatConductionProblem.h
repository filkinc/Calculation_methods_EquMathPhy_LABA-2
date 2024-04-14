#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>

using namespace std;

template<class T>
vector<T> right3diagLinSolve(const vector<T>& a, const vector<T>& b, const vector<T>& c, const vector<T>& d) {
    int n = d.size();
    vector<T> alpha(n + 1);
    vector<T> beta(n + 1);

    if (n == 1) {
        return { d[0] / b[0] };
    }

    alpha[1] = c[0] / b[0];
    beta[1] = d[0] / b[0];

    // прямой ход
    for (int i = 1; i < n; ++i) {
        alpha[i + 1] = c[i] / (b[i] - a[i] * alpha[i]);
        beta[i + 1] = (d[i] + a[i] * beta[i]) / (b[i] - a[i] * alpha[i]);
    }

    vector<T> x(n);
    x[n - 1] = beta[n];

    // обратный ход
    for (int i = n - 2; i >= 0; --i) {
        x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
    }

    return x;
}
