#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>

#include "HeatConductionProblem.h"

using namespace std;

template <class T>
using ThermCondBorderFunc = T(*)(T);//псевдоним типа (type alias)

template <class T>
using ThermCondCoefFunc = T(*)(T, T);

template <class T>
using ThermCondAFunc = T(*)(ThermCondCoefFunc<T>, T, T, T); // a(K(x, u), x[i - 1], x[i], u[i])

enum ThermCondBorderType {
    TEMP,
    HEAT_FLUX
};

template <class Type>
struct ThermCondTaskParams {
    Type L;
    Type T;

    Type c;
    Type rho;

    
    ThermCondBorderFunc<Type> initCond;

    ThermCondBorderType leftBorderType;
    ThermCondBorderFunc<Type> leftCond;

    ThermCondBorderType rightBorderType;
    ThermCondBorderFunc<Type> rightCond;

    ThermCondCoefFunc<Type> coefFunc;
};

template<class Type>
struct ThermCondMethodParams {
    int n; /*число отрезков по х*/
    int m; /*число отрезков по u*/

    ThermCondAFunc<Type> a;
};

template<class T>
void mixedLinThermalCondScheme(const ThermCondTaskParams<double>& taskParams, const ThermCondMethodParams<double>& methodParams, double sigma, char* nameFile) {
    int N = methodParams.n;
    int M = methodParams.m;

    double tau = taskParams.T / M;
    double h = taskParams.L / N;
    double gamma = h * h * taskParams.c * taskParams.rho / tau;

    vector<double> t(M + 1);
    t[0] = 0;
    for (int i = 1; i < M; ++i) {
        t[i] = t[0] + i * tau;
    }

    auto leftCond = taskParams.leftCond;
    auto rightCond = taskParams.rightCond;

    vector<double> a(N + 1);
    vector<double> x(N + 1);
    x[0] = 0.;
    for (int i = 1; i < N + 1; ++i) {
        x[i] = x[0] + i * h;
    }
    for (int i = 1; i < N + 1; ++i) {
        a[i] = methodParams.a(taskParams.coefFunc, x[i - 1], x[i], 0.);
    }

    vector<double> prev_y(N + 1);
    for (int i = 0; i < N + 1; ++i) {
        prev_y[i] = taskParams.initCond(x[i]);
    }

    vector<double> A(N), B(N), C(N), D(N);

    for (int j = 0; j < M; ++j) {

        for (int i = 1; i < N - 1; ++i) {
            A[i] = sigma * a[i];
            B[i] = gamma + sigma * (a[i + 1] + a[i]);
            C[i] = sigma * a[i + 1];
            D[i] = (1 - gamma) * a[i + 1] * prev_y[i + 1] + (gamma - (1 - gamma) * (a[i + 1] + a[i])) * prev_y[i] + (1 - gamma) * a[i] * prev_y[i - 1];
        }

        if (taskParams.leftBorderType == TEMP) {
            A[0] = 0;
            B[0] = -1;
            C[0] = 0;
            D[0] = -leftCond(t[j + 1]);
        }
        else if (taskParams.leftBorderType == HEAT_FLUX) {
            A[0] = 0;
            B[0] = gamma + sigma * 2 * a[1];
            C[0] = 2 * sigma * a[1];
            D[0] = (gamma - 2 * (1 - sigma) * a[1]) * prev_y[0] + 2 * (1 - sigma) * a[1] * prev_y[1] - 2 * sigma * h * leftCond(t[j + 1]) - 2 * (1 - sigma) * h * leftCond(t[j]);
        }

        if (taskParams.rightBorderType == TEMP) {
            A[N] = 0;
            B[N] = -1;
            C[N] = 0;
            D[N] = -rightCond(t[j + 1]);
        }
        else if (taskParams.rightBorderType == HEAT_FLUX) {
            A[N] = 2 * sigma * a[N];
            B[N] = gamma + 2 * sigma * a[N]; 
            C[N] = 0;
            D[N] = 2 * (1 - sigma) * a[N] * prev_y[N - 1] + (gamma - 2 * (1 - sigma) * a[N]) * prev_y[N] + 2 * sigma * h * rightCond(t[j + 1]) + 2 * (1 - sigma) * h * rightCond(t[j]);
        }

        auto cur_y = right3diagLinSolve(A, B, C, D);
        ofstream outFile;
        outFile.open(nameFile);
        for (int k = 0; k < N; ++k) {
            outFile << cur_y[k] << " ";
        }
        outFile << endl;
        outFile.close();
    }
}

template<typename T>
void explicitLinThermalCondScheme(const ThermCondTaskParams<double>& taskParams, const ThermCondMethodParams<double>& methodParams, char* nameFile) {
    int N = methodParams.n;
    int M = methodParams.m;

    double tau = taskParams.T / M;
    double h = taskParams.L / N;
    double gamma = h * h * taskParams.c * taskParams.rho / tau;

    vector<double> t(M + 1);
    t[0] = 0;
    for (int i = 1; i < M; ++i) {
        t[i] = t[0] + i * tau;
    }

    auto leftCond = taskParams.leftCond;
    auto rightCond = taskParams.rightCond;

    vector<double> a(N + 1);
    vector<double> x(N + 1);
    x[0] = 0.;
    for (int i = 1; i < N + 1; ++i) {
        x[i] = x[0] + i * h;
    }
    for (int i = 1; i < N + 1; ++i) {
        a[i] = methodParams.a(taskParams.coefFunc, x[i - 1], x[i], 0.);
    }

    vector<double> prev_y(N + 1);
    for (int i = 0; i < N + 1; ++i) {
        prev_y[i] = taskParams.initCond(x[i]);
    }

    vector<double> cur_y(N + 1);
     
    for (int j = 0; j < M; ++j) {
        for (int i = 1; i < N - 1; ++i) {
            cur_y[i] = a[i] / gamma * prev_y[i - 1] + (1 - a[i + 1] / gamma - a[i] / gamma) * prev_y[i] + a[i + 1] * prev_y[i + 1];
        }

        if(taskParams.leftBorderType == TEMP){
            cur_y[0] = leftCond(t[j + 1]);
        }
        else if (taskParams.leftBorderType == HEAT_FLUX) {
            cur_y[0] = (1 - 2 * a[1] / gamma) * prev_y[0] + 2 * a[1]/gamma * prev_y[1] - 2 * h /gamma * leftCond(t[j]);
        }

        if (taskParams.rightBorderType == TEMP) {
            cur_y[N] = rightCond(t[j + 1]);
        }
        else if (taskParams.rightBorderType == HEAT_FLUX) {
            cur_y[N] = 2 * a[N] / gamma * prev_y[N - 1] + (1 - 2 * a[N] / gamma) * prev_y[N] + 2 * h / gamma * rightCond(t[j]);
        }

        ofstream outFile;         
        outFile.open(nameFile);
            for (int k = 0; k < N; ++k) {
                outFile << cur_y[k] << " ";
            }
            outFile << endl;
        outFile.close();
    }
}

//макросы почитать про них

template<typename T>
T K( T x, T u ) {
    ThermCondTaskParams<T> taskParams;
    
    //Все коэффициенты нужно будет заменить на свой вариант
    taskParams.L = 10;
    T k_1 = 2,
      k_2 = 0.5, 
      x_1 = 1/2, 
      x_2 = 2/3;

    if (x >= 0 && x <= x_1)
        return k_1;
    else if (x > x_1 && x < x_2)
        return k_1 * (x - x_2) / (x_1 - x_2) + k_2 * (x - x_1) / (x_2 - x_1);
    else if (x >= x_2 && x <= taskParams.L)
        return k_2;
    else {
        cout << "х выходит за границы L!";
        return 0;
    }
}