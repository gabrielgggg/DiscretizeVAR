#pragma once
#include <armadillo>

using namespace arma;

namespace TimeSeries {

const double TStol = 1.0E-12;
const double pmSd = 3.0;
const double minPr = 5.0E-5;

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
DiscreteRV discretizeNormal(const double mean, const double sigma,
                            const uword n_elem = 21,
                            const double pmSigma = pmSd, bool cutTails = true);
vec discretizeNormal(vec grid, const double mean, const double sigma,
                     bool cutTails = true);

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
      : m_intercept(intercept), m_rho(rho), m_sigma(sigma)
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
      : m_intercept(intercept), m_rho(rho), m_sigma(sigma)
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
 * Orthogonalized VAR(1) and "rotation" matirx
 *
 */
enum OrthoMethod {  //
  SVD,
  Cholesky
};

class OrthogonalizedVAR {
 public:
  const VAR& getVAR() const { return m_ortho; }
  const VAR& getOriginalVAR() const { return m_original; }
  const mat& getSupportRotationMatrix() const { return LL; }

  OrthogonalizedVAR(VAR original, OrthoMethod method = OrthoMethod::Cholesky)
      : m_original(original), m_ortho(original)
  {
    const mat mEye = eye(original.size(), original.size());
    mat Atilde = original.intercept();
    mat Btilde = original.rho();
    mat GG = original.sigma();
    vec DD = ones(original.size());
    mat diagDD = mEye;
    LL = mEye;

    switch (method) {
      case OrthoMethod::SVD: {
        mat VV;
        svd(LL, DD, VV, GG);

        Atilde = LL.t() * original.intercept();
        Btilde = LL.t() * original.rho() * LL;
        diagDD = diagmat(DD);
        break;
      }
      case OrthoMethod::Cholesky: {
        LL = chol(GG, "lower");
        mat invLL = inv(LL);

        Atilde = invLL * original.intercept();
        Btilde = invLL * original.rho() * LL;
        break;
      }
      default: {
        cout << "Unknown method." << endl;
        exit(1);
      }
    }

    m_ortho = VAR(Atilde, Btilde, diagDD);
  }

 private:
  VAR m_original;
  VAR m_ortho;
  mat LL;
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

  MarkovChain() {}
  MarkovChain(rowvec support, mat transition)
      : m_support(support), m_tran(transition)
  {
    m_size = support.n_elem;
  }
  MarkovChain(const AR process, const uword sz = 21, const double pmMCsd = pmSd,
              bool expSupport = false);

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

  DiscreteVAR(const VAR var, const uword supportSize, bool trimGrids = true,
              OrthoMethod mm = OrthoMethod::Cholesky)
      : m_var(var)
  {
    m_supportSizes.set_size(m_var.size());
    m_supportSizes.fill(supportSize);
    impl(trimGrids, mm);
  }

  DiscreteVAR(const VAR var, const uvec gridSizes, bool trimGrids = true,
              OrthoMethod mm = OrthoMethod::Cholesky)
      : m_var(var), m_supportSizes(gridSizes)
  {
    impl(trimGrids, mm);
  }

 private:
  VAR m_var;
  MarkovChain m_mc;

  uword m_size;
  uword m_flatSize;
  uword m_midIx;

  mat m_grids;
  uvec m_supportSizes;

  void impl(bool trimGrids, OrthoMethod mm);
};

/*
 *
 * Discrete approximation to VAR(1) with common AR(1) stochastic volatility,
 * Basal-Yaron-like.
 *
 */
class StochasticVolVAR {
 public:
  void save(std::string fname) const;
  void print() const;

  StochasticVolVAR(const VAR var, const AR vol, const uvec gridSizes,
                   const uword volGridSize, bool trimGrids = true)
      : m_var(var)
      , m_vol(vol)
      , m_supportSizes(gridSizes)
      , m_volGridSize(volGridSize)
  {
    impl(trimGrids);
  };

 private:
  VAR m_var;
  AR m_vol;

  mat m_grids;
  uvec m_supportSizes;
  uword m_volGridSize;

  void impl(bool trimGrids);
};

}  // namespace TimeSeries
