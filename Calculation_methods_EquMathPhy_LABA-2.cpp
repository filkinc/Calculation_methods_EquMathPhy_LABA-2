#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>

#include "ConditionsHeatConductivity.h"
#include "HeatConductionProblem.h"

using namespace std;

const double pi = acos(-1);

// LINER
const string INPUT_TASK_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\file_task_params.txt";
const string INPUT_METHOD_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\file_method_params.txt";
const string OUTPUT_FILE = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\";

// QUAZILINER
//const string INPUT_TASK_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\quasi_file_task_params.txt";
//const string INPUT_METHOD_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\quasi_file_method_params.txt";
//const string OUTPUT_FILE = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\quasi_output_file.txt";


template<typename T>
ThermCondTaskParams<T> inputTaskParams(const string& nameFile) {
    ThermCondTaskParams<T> taskParams; 
    T lenght, time, capacity, density, startCond;
    int  leftBorderType, rightBorderType, initCondType, coefFuncType;

    ifstream inFile(nameFile);
    inFile >> lenght >> time >> capacity >> density >> initCondType >> leftBorderType >> rightBorderType >> coefFuncType;
    inFile.close();

    taskParams.L = lenght;
    taskParams.T = time;
    taskParams.c = capacity;
    taskParams.rho = density;

    //Выбор начальных условий 
    if (initCondType == 1) {
        taskParams.initCond = [](T) {return 0.2; };
    }
    else if (initCondType == 2) {
        taskParams.initCond = [](T x) {return sin(pi * x); };
    }

    //Выбор граничных условий
    //Слева
    if (leftBorderType == 1) {
        taskParams.leftCond = [](T) {return 0.2;};
        taskParams.leftBorderType = TEMP;
    }
    else if (leftBorderType == 2) {
        taskParams.leftCond = [](T) {return 0.; };
        taskParams.leftBorderType = TEMP;
    }
    else if (leftBorderType == 3) {
        taskParams.leftCond = [](T t) { 
            if (t >= 0 && t < 20) return 500.; // задание потока Q, температура которая выключается через t0 сек
            else return 0.;
            };
        taskParams.leftBorderType = HEAT_FLUX;
    }
    else if (leftBorderType == 4) {
        taskParams.leftCond = [](T t) {
            if (t >= 0 && t < 20) return 0.;
            else return 0.;
            };
        taskParams.leftBorderType = HEAT_FLUX;
    }

    //Справа
    if (rightBorderType == 1) {
        taskParams.rightCond = [](T t) {return 0.2; };
        taskParams.rightBorderType = TEMP;
    }
    else if (rightBorderType == 2) {
        taskParams.rightCond = [](T) {return 0.; };
        taskParams.rightBorderType = TEMP;
    }
    else if (rightBorderType == 3) {
        taskParams.rightCond = [](T t) {
            if (t >= 0 && t < 20) return 500.; // задание потока Q, температура которая выключается через t0 сек
            else return 0.;
            };
        taskParams.rightBorderType = HEAT_FLUX;
    }
    else if (rightBorderType == 4) {
        taskParams.rightCond = [](T t) {
            if (t >= 0 && t < 20) return 0.;
            else return 0.;
            };
        taskParams.rightBorderType = HEAT_FLUX;
    }

    //Выбор коэффициента теплопроводности
    if (coefFuncType == 1) {
        taskParams.coefFunc = K;
    }
    else if (coefFuncType == 2) {
        taskParams.coefFunc = constK;
    }
    else if (coefFuncType == 3) {
        taskParams.coefFunc = quasiK;
    }

    return taskParams;
}

template<typename T>
ThermCondMethodParams<T> inputMethodParams(const string& nameFile) {
    ThermCondMethodParams<T> methodParams;
    size_t n, m, ThermCondAFuncType;
    T sigma;
    int iterCount, idFunc;
    string newNameFile;

    ifstream inFile(nameFile);
    inFile >> newNameFile >> idFunc >> n >> m >> sigma >> ThermCondAFuncType >> iterCount;
    inFile.close();

    methodParams.newNameFile = newNameFile;
    methodParams.idFuncType = idFunc;
    methodParams.n = n;
    methodParams.m = m;
    methodParams.sigma = sigma;
    methodParams.iterCount = iterCount;

    if (ThermCondAFuncType == 1) {
        methodParams.a = a1;
    }
    else if (ThermCondAFuncType == 2) {
        methodParams.a = a2;
    }
    else if (ThermCondAFuncType == 3) {
        methodParams.a = a3;
    }

    return methodParams; 
 }

template<typename T>
void solveHeatPromblem(const ThermCondTaskParams<T>& taskParams, const ThermCondMethodParams<T>& methodParams, const string& nameFile) {
    if (methodParams.idFuncType == 1) {
        mixedLinThermalCondScheme<T>(taskParams, methodParams, nameFile);
    }
    else if (methodParams.idFuncType == 2) {
        quasiLinThermalCondScheme<T>(taskParams, methodParams, nameFile);
    }
 }

int main()
{
    auto taskParams = inputTaskParams<double>(INPUT_TASK_PARAMS);
    auto methodParams = inputMethodParams<double>(INPUT_METHOD_PARAMS);

    solveHeatPromblem(taskParams, methodParams, OUTPUT_FILE);

    return 0;
}
