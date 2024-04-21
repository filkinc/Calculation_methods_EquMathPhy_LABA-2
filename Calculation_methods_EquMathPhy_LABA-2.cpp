#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>

#include "ConditionsHeatConductivity.h"
#include "HeatConductionProblem.h"

using namespace std;

const string INPUT_TASK_PARAMS = "file_task_params.txt";
const string INPUT_METHOD_PARAMS = "file_method_params.txt";

template<typename T>
ThermCondTaskParams<T> inputTaskParams(const string& nameFile) {
    ThermCondTaskParams<T> taskParams; 
    T lenght, time, capacity, density, startCond;
    bool leftBorderType, rightBorderType, coefFuncType;

    ifstream inFile;
    inFile.open(nameFile);
    inFile >> lenght >> time >> capacity >> density >> leftBorderType >> rightBorderType >> coefFuncType;
    inFile.close();

    taskParams.L = lenght;
    taskParams.T = time;
    taskParams.c = capacity;
    taskParams.rho = density;
    taskParams.initCond = [](T) {return 0.2; };

    if (leftBorderType == true) { // TEMP (1)
        taskParams.leftCond = [](T) {return 0.2;};
    }
    else if (leftBorderType == false) {
        taskParams.rightCond = [](T t) { // FALSE (0)
            T t0 = 5, Q = 500; //потом сделать ввод этого из файла
            if (t >= 0 && t < t0)
                return Q;
            else
                return 0.;
            };
    }

    if (rightBorderType == TEMP) {
        taskParams.leftCond = [](T) {return 0.2; };
    }
    else if (rightBorderType == HEAT_FLUX) {
        taskParams.rightCond = [](T t) {
            T t0 = 5, Q = 500; //потом сделать ввод этого из файла
            if (t >= 0 && t < t0) return Q;
            else return 0.;
            };
    }

    if (coefFuncType == true) {
        taskParams.coefFunc = K;
    }
    else if (coefFuncType == false) {
        taskParams.coefFunc = quasiK;
    }

    return taskParams;
}

template<typename T>
ThermCondMethodParams<T> inputMethodParams(const string& nameFile) {
    ThermCondMethodParams<T> methodParams;
    size_t n, m, ThermCondAFuncType;

    ifstream inFile;
    inFile.open(nameFile);
    inFile >> n >> m >> ThermCondAFuncType;
    inFile.close();

    methodParams.n = n;
    methodParams.m = m;

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

int main()
{
    /*auto taskParams = inputTaskParams<double>("file_task_params.txt");*/
    auto taskParams = inputTaskParams<double>(INPUT_TASK_PARAMS);
    auto methodParams = inputMethodParams<double>(INPUT_METHOD_PARAMS);

    /*cout << taskParams.leftCond;*/
    /*mixedLinThermalCondScheme<double>(taskParams, methodParams, 1, "output_file_test");*/

    return 0;
}
