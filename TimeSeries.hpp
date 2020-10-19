#pragma once
#include <armadillo>

using namespace arma;

namespace TimeSeries {

const double TStol = 1.0E-10;
const double pmSd = 3.0;
const double minPr = 1.0E-4;

struct DiscreteRV {
  vec support;
  vec probabilities;
};

struct MvNormal {
  vec mean;
  mat varCovar;
};

/*
 * Normal helpers
 */
DiscreteRV discretizeNormal(const double mean,
    const double sigma,
    const uword n_elem = 21,
    const double pmSigma = pmSd, bool cutTails = true);
vec discretizeNormal(vec grid, const double mean, const double sigma, bool cutTails = true);

/*
 *
 * AR(1) process
 *
 */
class AR {
  public:
  double intercept() const { return m_intercept; }
  double rho() const { return m_rho; }
  double sigma() const { return m_sigma; }
  double stationaryMean() const { return m_intercept / (1.0 - m_rho); }
  double stationarySigma() const
  {
    return m_sigma / std::sqrt(1.0 - m_rho * m_rho);
  }
  bool stationaryQ() const { return std::abs(m_rho) < 1.0; }

  AR(double intercept, double rho, double sigma)
      : m_intercept(intercept)
      , m_rho(rho)
      , m_sigma(sigma)
  {
  }

  private:
  double m_intercept;
  double m_rho;
  double m_sigma;
};

/*
 *
 * VAR(1) process
 *
 */
class VAR {
  public:
  uword size() const { return m_size; }
  const vec& intercept() const { return m_intercept; }
  const mat& rho() const { return m_rho; }
  const mat& sigma() const { return m_sigma; }

  const MvNormal conditional(vec& current) const;
  vec& stationaryMean();
  mat& stationarySigma();
  bool stationaryQ();
  void print();

  VAR(vec intercept, mat rho, mat sigma)
      : m_intercept(intercept)
      , m_rho(rho)
      , m_sigma(sigma)
  {
    m_size = intercept.n_elem;
  }

  private:
  uword m_size;
  vec m_intercept;
  mat m_rho;
  mat m_sigma;

  vec m_statMean;
  mat m_statSigma;
  bool is_stationary;

  void findStationary();
};

/*
 *
 * Markov Chain
 *
 */
class MarkovChain {
  public:
  uword size() const { return m_size; }
  rowvec& support() { return m_support; }
  const rowvec& support() const { return m_support; }
  mat& transition() { return m_tran; }
  const mat& transition() const { return m_tran; }
  const rowvec& stationary();
  void save(std::string fname) const;

  MarkovChain() { }
  MarkovChain(rowvec support, mat transition)
      : m_support(support)
      , m_tran(transition)
  {
    m_size = support.n_elem;
  }
  MarkovChain(const AR process, const uword sz = 21, const double pmMCsd = pmSd);

  private:
  uword m_size;
  rowvec m_support;
  rowvec m_stationary;
  mat m_tran;
};

/*
 *
 * Discrete approximation of VAR(1)
 *
 */
class DiscreteVAR {
  public:
  uword size() const { return m_size; }
  uword flatSize() const { return m_flatSize; }
  uword midIx() const { return m_midIx; }
  const MarkovChain& markovChain() const { return m_mc; }
  const vec grid(uword varIx) const { return m_grids.col(varIx); }
  const rowvec gridValues(uword ix) const { return m_grids.row(ix); }
  void save(std::string fname) const;
  void print() const;

  DiscreteVAR(const VAR var, const uword supportSize, bool trimGrids = true)
      : m_var(var)
      , m_supportSize(supportSize)
  {
    impl(trimGrids);
  }

  private:
  VAR m_var;
  MarkovChain m_mc;

  uword m_supportSize;
  uword m_size;
  uword m_flatSize;
  uword m_midIx;

  mat m_grids;

  void impl(bool trimGrids);
};

} // namespace TimeSeries
