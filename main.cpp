#define ARMA_NO_DEBUG
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <iostream>

#include "TimeSeries.hpp"

using namespace arma;
using namespace TimeSeries;

int main()
{
    // example 3d
    vec icept = { 0.0, 0.25, 1.0 };
    mat rho = {
        { 0.9, 0.0, 0.0 }, //
        { 0.0, 0.8, 0.0 }, //
        { 0.0, 0.0, 0.0 }
    };
    mat sigma = {
        { 0.03 * 0.03, -0.0005, 0.0 },   //
        { -0.0005, 0.025 * 0.025, 0.0 }, //
        { 0.0, 0.0, 0.01 }
    };

    // example 2d
    /*
    vec icept = { 0.0, 0.25 };
    mat rho = {
        { 0.9, 0.0 }, //
        { 0.0, 0.8 }  //
    };
    mat sigma = {
        { 0.03 * 0.03, -0.0005 }, //
        { -0.0005, 0.02 * 0.02 }  //
    };
    */

    VAR myVarN(icept, rho, sigma);
    if (!myVarN.stationaryQ()) {
        cout << "Provided VAR is not stationary." << endl;
        return 1;
    }

    cout << "Start discretize..." << endl;

    uvec gridSizes = { 15, 15, 21 };
    DiscreteVAR myMCvar(myVarN, gridSizes, true, DiscreteVAR::Method::Cholesky);
    // myMCvar.print();
    myMCvar.save("myMCvar.h5");

    cout << "Done." << endl;
    return 0;
}
