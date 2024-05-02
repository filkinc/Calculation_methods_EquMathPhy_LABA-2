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
using ThermCondAFunc = T(*)(ThermCondCoefFunc<T>, T, T, T, T); // a(K(x, u), x[i - 1], x[i], u[i - 1], u[i])

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
    string newNameFile;
    bool idFuncType;
    int n; /*число отрезков по х*/
    int m; /*число отрезков по u*/
    Type sigma;
    int iterCount;

    ThermCondAFunc<Type> a;
};

template<typename T>
T K(T x, T u) {
    //Данные для нашего варианта
    T k_1 = 2,
        k_2 = 0.5,
        x_1 = 1 / 2,
        x_2 = 2 / 3;

    if (x >= 0 && x <= x_1)
        return k_1;
    else if (x > x_1 && x < x_2)
        return k_1 * (x - x_2) / (x_1 - x_2) + k_2 * (x - x_1) / (x_2 - x_1);
    else if (x >= x_2)
        return k_2;
}

template<typename T>
T constK(T x, T u) {
    return 1;
}

template<typename T>
T quasiK(T x, T u) {
    //Данные для нашего варианта
    return 0.5 + 2 * u * u;
}

template<typename T>
T a1(ThermCondCoefFunc<T> K, T prev_x, T cur_x, T prev_y, T cur_y) { // Здесь prev и cur по пространству
    return 0.5 * (K(cur_x, cur_y) + K(prev_x, prev_y));
}

template<typename T>
T a2(ThermCondCoefFunc<T> K, T prev_x, T cur_x, T prev_y, T cur_y) {
    return K(cur_x - 0.5 * (cur_x - prev_x), cur_y);
}

template<typename T>
T a3(ThermCondCoefFunc<T> K, T prev_x, T cur_x, T prev_y, T cur_y) {
    return sqrt(K(cur_x, cur_y) * K(prev_x, prev_y));
}

//template<typename T>
//T aQuasi(ThermCondCoefFunc<T>& K, T prev_x, T cur_x, T prev_u, T cur_u) {
//    return 0.5 * (K(cur_u) + K(prev_u));
//} 

template<class T>
void mixedLinThermalCondScheme(const ThermCondTaskParams<T>& taskParams, const ThermCondMethodParams<T>& methodParams, const string& nameFile) {
    int N = methodParams.n;
    int M = methodParams.m;

    T sigma = methodParams.sigma;

    T tau = taskParams.T / M;
    T h = taskParams.L / N;
    T gamma = h * h * taskParams.c * taskParams.rho / tau;

    vector<T> t(M + 1);
    t[0] = 0;
    for (int i = 1; i < M + 1; ++i) {
        t[i] = t[0] + i * tau;
    }

    auto leftCond = taskParams.leftCond;
    auto rightCond = taskParams.rightCond;

    vector<T> a(N + 1);
    vector<T> x(N + 1);
    x[0] = 0.;
    for (int i = 1; i < N + 1; ++i) {
        x[i] = x[0] + i * h;
    }
   
    vector<T> prev_y(N + 1);
    for (int i = 0; i < N + 1; ++i) {
        prev_y[i] = taskParams.initCond(x[i]);
    }

    for (int i = 1; i < N + 1; ++i) {
        a[i] = methodParams.a(taskParams.coefFunc, x[i - 1], x[i], prev_y[i - 1], prev_y[i]);
    }

    vector<T> cur_y(N + 1);
    vector<T> A(N + 1), B(N + 1), C(N + 1), D(N + 1);

    ofstream outFile(nameFile + methodParams.newNameFile);
    outFile.precision(8);

    outFile << 0 << " ";
    for (int i = 0; i < N + 1; ++i) {
        outFile << x[i] << " ";
    }
    outFile << endl;

    /*outFile.open;*/
    
    outFile << t[0] << " ";
    for (int k = 0; k < N + 1; ++k) { // законментил чтобы не выводились начальные условия
        outFile << prev_y[k] << " ";
    }
    outFile << endl;

    //for (int i = 0; i < M + 1; ++i) { //Запись в файл времени
    //    outFile << t[i] << " ";
    //}
    //outFile << endl;

    //for (int i = 0; i < N + 1; ++i) { // Запись в файл пространственных координат
    //    outFile << x[i] << " ";
    //}

    for (int j = 0; j < M; ++j) {

        for (int i = 1; i < N; ++i) {
            A[i] = sigma * a[i];
            B[i] = gamma + sigma * (a[i + 1] + a[i]);
            C[i] = sigma * a[i + 1];
            D[i] = (1 - sigma) * a[i + 1] * prev_y[i + 1] + (gamma - (1 - sigma) * (a[i + 1] + a[i])) * prev_y[i] + (1 - sigma) * a[i] * prev_y[i - 1];
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

        cur_y = right3diagLinSolve(A, B, C, D);

        outFile << t[j + 1] << " ";
        
        for (int k = 0; k < N + 1; ++k) {
            outFile << cur_y[k] << " ";
        }
        outFile << endl;

        for (int i = 0; i < N + 1; ++i) {
            prev_y[i] = cur_y[i];
        }
    }
    outFile.close();

    /*ofstream outFileParams(nameFile + methodParams.newNameFile.substr(0, methodParams.newNameFile.length() - 4) + "_task_params.txt");
    outFileParams << taskParams.T << " " << taskParams.L;
    outFileParams.close();*/
}

template<typename T>
void explicitLinThermalCondScheme(const ThermCondTaskParams<T>& taskParams, const ThermCondMethodParams<T>& methodParams, const string& nameFile) {
    int N = methodParams.n;
    int M = methodParams.m;

    T tau = taskParams.T / M;
    T h = taskParams.L / N;
    T gamma = h * h * taskParams.c * taskParams.rho / tau;

    vector<T> t(M + 1);
    t[0] = 0;
    for (int i = 1; i < M; ++i) {
        t[i] = t[0] + i * tau;
    }

    auto leftCond = taskParams.leftCond;
    auto rightCond = taskParams.rightCond;

    vector<T> a(N + 1);
    vector<T> x(N + 1);
    x[0] = 0.;
    for (int i = 1; i < N + 1; ++i) {
        x[i] = x[0] + i * h;
    }
    
    vector<T> prev_y(N + 1);
    for (int i = 0; i < N + 1; ++i) {
        prev_y[i] = taskParams.initCond(x[i]);
    }

    for (int i = 1; i < N + 1; ++i) {
        a[i] = methodParams.a(taskParams.coefFunc, x[i - 1], x[i], prev_y[i - 1], prev_y[i]);
    }

    vector<T> cur_y(N + 1);
    
    ofstream outFile(nameFile + methodParams.newNameFile);

    /*outFile.open;*/

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

        for (int k = 0; k < N + 1; ++k) {
            outFile << cur_y[k] << " ";
        }
        outFile << endl;

        for (int i = 0; i < N + 1; ++i) {
            prev_y[i] = cur_y[i];
        }
    }
    outFile.close();
}


template<typename T>
void quasiLinThermalCondScheme(const ThermCondTaskParams<T>& taskParams, const ThermCondMethodParams<T>& methodParams, const string& nameFile) {
    int N = methodParams.n;
    int M = methodParams.m;
    int iterCount = methodParams.iterCount;

    T tau = taskParams.T / M;
    T h = taskParams.L / N;
    T gamma = h * h * taskParams.c * taskParams.rho / tau;

    vector<T> t(M + 1);
    t[0] = 0;
    for (int i = 1; i < M + 1; ++i) {
        t[i] = t[0] + i * tau;
    }

    auto leftCond = taskParams.leftCond;
    auto rightCond = taskParams.rightCond;

    vector<T> a(N + 1);
    vector<T> x(N + 1);
    x[0] = 0.;
    for (int i = 1; i < N + 1; ++i) {
        x[i] = x[0] + i * h;
    }
    
    vector<T> prev_y(N + 1);
    for (int i = 0; i < N + 1; ++i) {
        prev_y[i] = taskParams.initCond(x[i]);
    }
    for (int i = 1; i < N + 1; ++i) {
        a[i] = methodParams.a(taskParams.coefFunc, x[i - 1], x[i], prev_y[i - 1], prev_y[i]);
    }

    vector<T> cur_y(N + 1);
    vector<T> A(N + 1), B(N + 1), C(N + 1), D(N + 1);

    ofstream outFile(nameFile + methodParams.newNameFile);
    /*outFile.open;*/

    for (int i = 0; i < N + 1; ++i) { // законментил чтобы не выводились начальные условия
        outFile << prev_y[i] << " ";
    }
    outFile << endl;
    
    vector<int> iterVec(M + 1);

    for (int j = 0; j < M; ++j) {
        int iter = 0;

        vector<T> prev_appr_y = prev_y;

            while (true) {
                for (int i = 1; i < N; ++i) {
                    A[i] = a[i];
                    B[i] = gamma + a[i + 1] + a[i];
                    C[i] = a[i + 1];
                    D[i] = gamma * prev_y[i];
                }

                if (taskParams.leftBorderType == TEMP) {
                    A[0] = 0;
                    B[0] = -1;
                    C[0] = 0;
                    D[0] = -leftCond(t[j + 1]);
                }
                else if (taskParams.leftBorderType == HEAT_FLUX) {
                    A[0] = 0;
                    B[0] = gamma + 2 * a[1];
                    C[0] = 2 * a[1];
                    D[0] = gamma * prev_y[0] - 2 * h * leftCond(t[j + 1]);
                }

                if (taskParams.rightBorderType == TEMP) {
                    A[N] = 0;
                    B[N] = -1;
                    C[N] = 0;
                    D[N] = -rightCond(t[j + 1]);
                }
                else if (taskParams.rightBorderType == HEAT_FLUX) {
                    A[N] = 2 * a[N];
                    B[N] = gamma + 2 * a[N];
                    C[N] = 0;
                    D[N] = gamma * prev_y[N] + 2 * h * rightCond(t[j + 1]);
                }

                prev_appr_y = move(cur_y);
                cur_y = right3diagLinSolve(A, B, C, D);

                for (int k = 1; k < N + 1; ++k) {
                    a[k] = methodParams.a(taskParams.coefFunc, x[k - 1], x[k], cur_y[k - 1], cur_y[k]);
                }

                T norm = norm_2(diff(cur_y, prev_appr_y));
                ++iter;

                if (iterCount == 0 && norm_2(diff(cur_y, prev_appr_y)) < 1e-9 && iter <= 100) {
                    break;
                }
                else if (iterCount > 0 && iter >= iterCount) {
                    break;
                }
            }

            iterVec[j] = iter;
        
        for (int k = 0; k < N + 1; ++k) {
            outFile << cur_y[k] << " ";
        }
        outFile << endl;

        for (int i = 0; i < N + 1; ++i) {
            prev_y[i] = cur_y[i];
        }
    }
    outFile.close();

    ofstream outFileIter(nameFile + methodParams.newNameFile.substr(0, methodParams.newNameFile.length() - 4) + "_iter.txt");
    for (int i = 0; i < M; ++i) {
        outFileIter << iterVec[i] << " ";
    }
    outFileIter.close();

    ofstream outFileParams(nameFile + methodParams.newNameFile.substr(0, methodParams.newNameFile.length() - 4) + "_task_params.txt");
    outFileParams << taskParams.T << " " << taskParams.L;
    outFileParams.close();
    
}
