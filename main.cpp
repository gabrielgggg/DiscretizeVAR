#define ARMA_NO_DEBUG
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <iostream>

#include "TimeSeries.hpp"

using namespace arma;
using namespace TimeSeries;

void markovExample();
void discretizeVARexample();
void discretizeStochasticVolVARexample();

int main()
{
  // markovExample();
  discretizeVARexample();
  discretizeStochasticVolVARexample();
  return 0;
}

/*
 *
 * Examples
 *
 */

void markovExample()
{
  cout << "Markov Example." << endl;
  AR vol(0.0, 0.95, 0.005);
  const MarkovChain mc(vol, 5, 3.0, true);
  mc.print();
}

void discretizeVARexample()
{
  cout << "Discretize VAR example. Results saved in myMCvar.h5" << endl;

  // example 3d
  vec icept = { 0.0, 0.25, 1.0 };
  mat rho = { { 0.97, 0.0, 0.0 }, //
    { 0.0, 0.98, 0.0 }, //
    { 0.0, 0.0, 0.3 } };
  mat sigma = { { 0.03 * 0.03, -0.0005, 0.0 }, //
    { -0.0005, 0.025 * 0.025, 0.0 }, //
    { 0.0, 0.0, 0.01 } };

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
    cout << "Provided VAR is not stationary. Stopping." << endl;
    return;
  }

  uvec gridSizes = { 15, 15, 15 };
  DiscreteVAR myMCvar(myVarN, gridSizes, true, OrthoMethod::Cholesky);
  myMCvar.save("myMCvar.h5");
}

void discretizeStochasticVolVARexample()
{
  cout << "Discretize stochatic VAR example. Results saved in stochVol.h5"
       << endl;

  AR vol(0.0, 0.9, 0.05);

  vec icept = { 0.0, 0.0 };
  mat rho = { { 0.9, 0.0 }, //
    { 0.0, 0.9 } };
  mat sigma = { { 0.015 * 0.015, -0.3 * 0.015 * 0.015 }, //
    { -0.3 * 0.015 * 0.015, 0.015 * 0.015 } };
  VAR atMeanVolVAR(icept, rho, sigma);

  uvec gridSizes = { 15, 15 };
  uword volGridSize = 15;
  DiscreteStochVolVAR volVAR(atMeanVolVAR, vol, gridSizes, volGridSize, true);
  volVAR.save("stochVol.h5");
}
