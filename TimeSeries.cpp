// define ARMA_NO_DEBUG
#define ARMA_ALLOW_FAKE_GCC
#include "TimeSeries.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>

using namespace arma;

namespace TimeSeries {

/*
 * Discretize Normal for given grid.
 */
vec discretizeNormal(vec grid, const double mean, const double sigma,
                     bool cutTails)
{
  uword sz = grid.n_elem;
  vec prs(sz);

  if (sz == 1) {
    prs(0) = 1.0;
    return prs;
  } else if (sz == 2) {
    prs(0) = normcdf((grid(0) + grid(1)) / 2.0, mean, sigma);
    prs(1) = 1.0 - prs(0);
    return prs;
  }

#pragma omp parallel for
  for (uword ix = 1; ix < sz - 1; ++ix)
    prs(ix) = normcdf((grid(ix) + grid(ix + 1)) / 2.0, mean, sigma) -
              normcdf((grid(ix - 1) + grid(ix)) / 2.0, mean, sigma);

  if (cutTails) {
    const double stepOut = (grid(1) - grid(0)) / 2.0;
    prs(0) = normcdf(grid(0) + stepOut, mean, sigma) -
             normcdf(grid(0) - stepOut, mean, sigma);
    prs(sz - 1) = normcdf(grid(sz - 1) + stepOut, mean, sigma) -
                  normcdf(grid(sz - 1) - stepOut, mean, sigma);

    prs /= sum(prs);
  } else {
    prs(0) = normcdf((grid(0) + grid(1)) / 2.0, mean, sigma);
    prs(sz - 1) =
        1.0 - normcdf((grid(sz - 2) + grid(sz - 1)) / 2.0, mean, sigma);
  }

  return prs;
}

/*
 * Discretize Normal over uniformly spaced grid
 * spanning +/- pmSigma standard deviations.
 */
DiscreteRV discretizeNormal(const double mean, const double sigma,
                            const uword n_elem, const uword pmSigma,
                            bool cutTails)
{
  DiscreteRV drv;
  drv.support =
      linspace<vec>(mean - pmSigma * sigma, mean + pmSigma * sigma, n_elem);
  drv.probabilities = discretizeNormal(drv.support, mean, sigma, cutTails);
  return drv;
}

/*
 *
 * VAR(1)
 *
 */
const MvNormal VAR::conditional(vec& current) const
{
  MvNormal mvn;
  mvn.mean = m_intercept + m_rho * current;
  mvn.varCovar = m_sigma;
  return mvn;
}

void VAR::findStationary()
{
  const mat eyeRho = eye(m_size, m_size) - m_rho;
  double theDet = det(eyeRho);

  if (fabs(theDet) < TStol) {
    cout << "Determinant " << det(eyeRho) << endl;
    is_stationary = false;
    m_statMean.reset();
    m_statSigma.reset();
    return;
  } else
    is_stationary = true;

  m_statMean = solve(eyeRho, m_intercept);

  mat v0 = m_sigma;
  mat v1 = m_sigma;
  v0.fill(0.0);
  double err = 1.0;
  while (err > TStol) {
    v1 = m_rho * v0 * m_rho.t() + m_sigma;
    err = abs(v1 - v0).max();
    v0 = v1;
  }
  m_statSigma = v1;
}

vec& VAR::stationaryMean()
{
  if (m_statMean.is_empty()) findStationary();
  return m_statMean;
}

mat& VAR::stationarySigma()
{
  if (m_statSigma.is_empty()) findStationary();
  return m_statSigma;
}

bool VAR::stationaryQ()
{
  if (m_statMean.is_empty()) findStationary();
  return is_stationary;
}

void VAR::print()
{
  cout << "Intercept: " << endl;
  intercept().print();
  cout << "Rho: " << endl;
  rho().print();
  cout << "Sigma: " << endl;
  sigma().print();

  cout << "Stationary: " << endl;
  stationaryMean().print();
  cout << "Stationary Sigma: " << endl;
  stationarySigma().print();
}

/*
 *
 * Markov Chain
 *
 */
MarkovChain::MarkovChain(const AR process, const uword sz, const double pmMCsd,
                         bool expSupport)
{
  m_size = sz;
  double pMean = process.stationaryMean();
  double pSd = process.stationarySigma();
  m_support = linspace<rowvec>(pMean - pmMCsd * pSd, pMean + pmMCsd * pSd, sz);

  m_tran.set_size(sz, sz);
#pragma omp parallel for
  for (uword rIx = 0; rIx < sz; ++rIx) {
    double condMean = process.intercept() + process.rho() * m_support(rIx);
    vec condPrs = discretizeNormal(m_support.t(), condMean, process.sigma());
    m_tran.row(rIx) = condPrs.t();
  }

  if (expSupport) m_support = exp(m_support);
}

rowvec stationaryDistribution(const mat& transitionMatrix)
{
  rowvec stat;
  cx_vec eigval;
  cx_mat eigvec;
  // eig_gen(eigval, eigvec, m_tran.t(), "balance");
  eigs_gen(eigval, eigvec, sp_mat(transitionMatrix.t()), 1, 1.0);

  uword unitEigIx = (abs(abs(eigval) - 1.0)).index_min();

  stat = abs(eigvec.col(unitEigIx)).t();
  return stat / sum(stat);
}

const rowvec& MarkovChain::stationary()
{
  if (m_stationary.is_empty()) {
    m_stationary = stationaryDistribution(m_tran);
  }

  return m_stationary;
}

void MarkovChain::save(const std::string fname) const
{
  m_support.save(hdf5_name(fname, "MC/support", hdf5_opts::replace));
  m_stationary.save(hdf5_name(fname, "MC/stationary", hdf5_opts::replace));
  m_tran.save(hdf5_name(fname, "MC/transitions", hdf5_opts::replace));
}

void MarkovChain::print() const
{ /* TODO */
}

void trimMarkovChain(mat& grids, mat& transition)
{
  /*
   * Based on Gordon 2020 wip.
   */
  rowvec probs = stationaryDistribution(transition);
  uword sz = probs.n_elem;

  uvec removePoints;
  uword removeCount = 0;
  for (uword ix = 0; ix < sz; ++ix)
    if (probs(ix) <= minPr) {
      removePoints.resize(removeCount + 1);
      removePoints(removeCount) = ix;
      removeCount++;
    }

  uword newSz = sz - removeCount;
  cout << "Removing " << removeCount << " points out of " << sz << "." << endl;

  grids.shed_rows(removePoints);
  transition.shed_rows(removePoints);
  transition.shed_cols(removePoints);
#pragma omp parallel for
  for (uword fIx = 0; fIx < newSz; ++fIx)
    transition.row(fIx) = transition.row(fIx) / sum(transition.row(fIx));
}

/*
 *
 * Discretized VAR(1)
 *
 */
void DiscreteVAR::impl(bool trimGrids, OrthoMethod method)
{
  // Sizing
  m_size = m_var.size();
  m_flatSize = prod(m_supportSizes);
  m_grids.set_size(m_flatSize, m_size);
  umat m_map(m_flatSize, m_size);

  OrthogonalizedVAR ortho(m_var, method);
  mat LL = ortho.getSupportRotationMatrix();
  VAR orthog = ortho.getVAR();

  const vec Atilde = orthog.intercept();
  const mat Btilde = orthog.rho();
  const vec DD = orthog.sigma().diag();

  const vec uncondE = orthog.stationaryMean();
  const vec uncondV = orthog.stationarySigma().diag();

  // Prepare logical grids
  field<vec> m_orthogGrids(m_size);
#pragma omp parallel for
  for (uword vIx = 0; vIx < m_size; ++vIx)
    m_orthogGrids(vIx) = linspace<vec>(uncondE(vIx) - pmSd * sqrt(uncondV(vIx)),
                                       uncondE(vIx) + pmSd * sqrt(uncondV(vIx)),
                                       m_supportSizes(vIx));

    // Prepare flat grids
#pragma omp parallel for
  for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx) {
    uword tmp = flatIx;
    for (uword vIx = 0; vIx < m_size; ++vIx) {
      const vec logical = m_orthogGrids(vIx);
      m_map(flatIx, vIx) = tmp % m_supportSizes(vIx);
      m_grids(flatIx, vIx) = logical(m_map(flatIx, vIx));
      tmp = tmp / m_supportSizes(vIx);
    }
  }

  // Transition matrix
  m_tran.set_size(m_flatSize, m_flatSize);
  for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx) {
    const rowvec tmpVal = m_grids.row(flatIx);
    const vec condMeans = Atilde + Btilde * tmpVal.t();
    const vec condSd = sqrt(DD);
    field<vec> condPrs(m_size);
    for (uword vIx = 0; vIx < m_size; ++vIx)
      condPrs(vIx) =
          discretizeNormal(m_orthogGrids(vIx), condMeans(vIx), condSd(vIx));

#pragma omp parallel for
    for (uword flatPrIx = 0; flatPrIx < m_flatSize; ++flatPrIx) {
      double goPr = 1.0;
      for (uword vPrIx = 0; vPrIx < m_size; ++vPrIx) {
        vec tmpDist = condPrs(vPrIx);
        goPr *= tmpDist(m_map(flatPrIx, vPrIx));
      }
      m_tran(flatIx, flatPrIx) = goPr;
    }
  }

  // Adjust grids
  m_grids = (LL * m_grids.t()).t();
  cout << "Completed discretization." << endl;

  if (trimGrids) {
    trimMarkovChain(m_grids, m_tran);
    m_flatSize = m_tran.n_rows;
  }

  // Finding midIx
  vec distances(m_flatSize);
  const vec target = m_var.stationaryMean();
#pragma omp parallel for
  for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx) {
    const vec tmp = target - m_grids.row(flatIx).t();
    distances(flatIx) = sum(tmp % tmp);
  }
  m_midIx = distances.index_min();
}

void DiscreteVAR::print() const
{
  cout << "Grid:" << endl;
  m_grids.print();
  cout << endl << "Transition matrix: " << endl;
  m_tran.print();
}

void DiscreteVAR::save(const std::string fname) const
{
  uvec tmp = {m_midIx};
  tmp.save(hdf5_name(fname, "DVAR/midIx", hdf5_opts::replace));
  m_grids.save(hdf5_name(fname, "DVAR/grids", hdf5_opts::replace));
  m_tran.save(hdf5_name(fname, "DVAR/transitions", hdf5_opts::replace));
  m_var.intercept().save(
      hdf5_name(fname, "DVAR/varIntercept", hdf5_opts::replace));
  m_var.rho().save(hdf5_name(fname, "DVAR/varRho", hdf5_opts::replace));
  m_var.sigma().save(hdf5_name(fname, "DVAR/varSigma", hdf5_opts::replace));
}

/*
 *
 * Stochastic volatility VAR
 *
 */
void DiscreteStochVolVAR::impl(bool trimGrids)
{
  m_size = m_var.size();
  uword justDimsSz = prod(m_supportSizes);
  m_flatSize = justDimsSz * m_volGridSize;
  m_grids.set_size(m_flatSize, m_size + 1);
  umat m_map(m_flatSize, m_size + 1);

  MarkovChain volMC(m_vol, m_volGridSize, pmSd, true);
  const rowvec& vols = volMC.support();
  const mat& volTran = volMC.transition();

  OrthogonalizedVAR ortho(m_var);  // Default to Cholesky
  mat rotationLL = ortho.getSupportRotationMatrix();
  VAR orthog = ortho.getVAR();
  const vec uncondE = orthog.stationaryMean();
  const vec uncondV = orthog.stationarySigma().diag();

  const vec Atilde = orthog.intercept();
  const mat Btilde = orthog.rho();
  const vec DD = orthog.sigma().diag();

  // Prepare logical grids
  field<vec> m_orthogGrids(m_size, m_volGridSize);
#pragma omp parallel for collapse(2)
  for (uword volIx = 0; volIx < m_volGridSize; ++volIx)
    for (uword iIx = 0; iIx < m_size; ++iIx)
      m_orthogGrids(iIx, volIx) =
          linspace<vec>(uncondE(iIx) - pmSd * vols(volIx) * sqrt(uncondV(iIx)),
                        uncondE(iIx) + pmSd * vols(volIx) * sqrt(uncondV(iIx)),
                        m_supportSizes(iIx));

      // Prepare flat grids
#pragma omp parallel for collapse(2)
  for (uword volIx = 0; volIx < m_volGridSize; ++volIx) {
    for (uword flatIx = 0; flatIx < justDimsSz; ++flatIx) {
      uword withV = justDimsSz * volIx + flatIx;

      m_map(withV, m_size) = volIx;
      m_grids(withV, m_size) = vols(volIx);

      uword tmp = flatIx;
      for (uword iIx = 0; iIx < m_size; ++iIx) {
        const vec logical = m_orthogGrids(iIx, volIx);
        m_map(withV, iIx) = tmp % m_supportSizes(iIx);
        m_grids(withV, iIx) = logical(m_map(withV, iIx));
        tmp /= m_supportSizes(iIx);
      }
    }
  }

  // Transition matrix
  m_tran.set_size(m_flatSize, m_flatSize);
  m_tran.fill(0.0);
#pragma omp parallel for
  for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx) {
    const uword thisVolIx = m_map(flatIx, m_size);
    const rowvec thisValues = m_grids.row(flatIx).head(m_size);
    const vec condMeans = Atilde + Btilde * thisValues.t();
    const vec condSdReference = sqrt(DD);

    field<vec> condPrs(m_size, m_volGridSize);
    for (uword iIx = 0; iIx < m_size; ++iIx)
      for (uword vPrIx = 0; vPrIx < m_volGridSize; ++vPrIx)
        condPrs(iIx, vPrIx) =
            discretizeNormal(m_orthogGrids(iIx, vPrIx), condMeans(iIx),
                             vols(vPrIx) * condSdReference(iIx));

    for (uword flatPrIx = 0; flatPrIx < m_flatSize; ++flatPrIx) {
      const uword vPrIx = m_map(flatPrIx, m_size);
      double goPr = 1.0;
      for (uword iiIx = 0; iiIx < m_size; ++iiIx) {
        const vec condPrVec = condPrs(iiIx, vPrIx);
        goPr *= condPrVec(m_map(flatPrIx, iiIx));
      }
      m_tran(flatIx, flatPrIx) = goPr * volTran(thisVolIx, vPrIx);
    }
  }

  // Adjust grids
  m_grids.cols(0, m_size - 1) =
      (rotationLL * m_grids.cols(0, m_size - 1).t()).t();
  cout << "Completed discretization." << endl;

  if (trimGrids) {
    trimMarkovChain(m_grids, m_tran);
    m_flatSize = m_tran.n_rows;
  }

  // Finding midIx
  vec distances(m_flatSize);
  vec uncondEextended(m_size + 1);
  uncondEextended.head(m_size) = m_var.stationaryMean();
  uncondEextended(m_size) = vols(m_volGridSize / 2);
#pragma omp parallel for
  for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx) {
    const vec tmp = uncondEextended - m_grids.row(flatIx).t();
    distances(flatIx) = sum(tmp % tmp);
  }
  m_midIx = distances.index_min();
}

void DiscreteStochVolVAR::save(const std::string fname) const
{
  uvec tmp = {m_midIx};
  tmp.save(hdf5_name(fname, "DsvVAR/midIx", hdf5_opts::replace));
  m_grids.save(hdf5_name(fname, "DsvVAR/grids", hdf5_opts::replace));
  m_tran.save(hdf5_name(fname, "DsvVAR/transitions", hdf5_opts::replace));
  m_var.intercept().save(
      hdf5_name(fname, "DsvVAR/varIntercept", hdf5_opts::replace));
  m_var.rho().save(hdf5_name(fname, "DsvVAR/varRho", hdf5_opts::replace));
  m_var.sigma().save(hdf5_name(fname, "DsvVAR/varSigma", hdf5_opts::replace));

  vec volAR = {m_vol.intercept(), m_vol.rho(), m_vol.sigma()};
  volAR.save(
      hdf5_name(fname, "DsvVAR/volatilityAR1params", hdf5_opts::replace));
}

void DiscreteStochVolVAR::print() const
{ /* TODO */
}

}  // namespace TimeSeries
