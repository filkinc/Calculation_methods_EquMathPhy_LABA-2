#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>

#include "ConditionsHeatConductivity.h"
#include "HeatConductionProblem.h"

using namespace std;

// LINER
const string INPUT_TASK_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\file_task_params.txt";
const string INPUT_METHOD_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\file_method_params.txt";
const string OUTPUT_FILE = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\output_file.txt";

// QUAZILINER
//const string INPUT_TASK_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\quasi_file_task_params.txt";
//const string INPUT_METHOD_PARAMS = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\quasi_file_method_params.txt";
//const string OUTPUT_FILE = "C:\\Users\\filki\\Documents\\Учеба в МГТУ\\6 семестр\\Численные методы МФ\\Лаба №2 ЧМ\\Calculation_methods_EquMathPhy_LABA-2\\quasi_output_file.txt";


template<typename T>
ThermCondTaskParams<T> inputTaskParams(const string& nameFile) {
    ThermCondTaskParams<T> taskParams; 
    T lenght, time, capacity, density, startCond;
    bool leftBorderType, rightBorderType, coefFuncType;

    ifstream inFile(nameFile);
    /*inFile.open;*/
    inFile >> lenght >> time >> capacity >> density >> leftBorderType >> rightBorderType >> coefFuncType;
    inFile.close();

    taskParams.L = lenght;
    taskParams.T = time;
    taskParams.c = capacity;
    taskParams.rho = density;
    taskParams.initCond = [](T) {return 0.2; };

    if (leftBorderType == true) { // TEMP (1)
        taskParams.leftCond = [](T) {return 0.2;};
        taskParams.leftBorderType = TEMP;
    }
    else if (leftBorderType == false) {
        taskParams.leftCond = [](T t) { // FALSE (0)
            if (t >= 0 && t < 20) return 500.;
            else return 0.;
            };
        taskParams.leftBorderType = HEAT_FLUX;
    }

    if (rightBorderType == true) {
        taskParams.rightCond = [](T t) {return 0.2; };
        taskParams.rightBorderType = TEMP;
    }
    else if (rightBorderType == false) {
        taskParams.rightCond = [](T t) {
            if (t >= 0 && t < 20) return 500.; // задание потока Q, температура которая выключается через t0 сек
            else return 0.;
            };
        taskParams.rightBorderType = HEAT_FLUX;
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
    T sigma;
    int iterCount;

    ifstream inFile(nameFile);
    /*inFile.open;*/
    inFile >> n >> m >> sigma >> ThermCondAFuncType >> iterCount;
    inFile.close();

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

int main()
{
    auto taskParams = inputTaskParams<double>(INPUT_TASK_PARAMS);
    auto methodParams = inputMethodParams<double>(INPUT_METHOD_PARAMS);

    /*cout << taskParams.leftCond << endl;
    cout << taskParams.rightCond;*/

    // LINER
    //mixedLinThermalCondScheme<double>(taskParams, methodParams, OUTPUT_FILE);

    //QUAIZILINER
    quasiLinThermalCondScheme<double>(taskParams, methodParams, OUTPUT_FILE);



    return 0;
}
