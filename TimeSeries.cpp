#include "TimeSeries.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>

using namespace arma;
using namespace std;

namespace TimeSeries {

/*
 * Discretize Normal for given grid.
 */
vec discretizeNormal(vec grid, double mean, double sigma)
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

    prs(0) = normcdf((grid(0) + grid(1)) / 2.0, mean, sigma);
    prs(sz - 1) = 1.0 - normcdf((grid(sz - 2) + grid(sz - 1)) / 2.0, mean, sigma);

    for (uword ix = 1; ix < sz - 1; ++ix) {
        prs(ix) = normcdf((grid(ix) + grid(ix + 1)) / 2.0, mean, sigma) - normcdf((grid(ix - 1) + grid(ix)) / 2.0, mean, sigma);
    }

    return prs;
}

/*
 * Discretize Normal over uniformly spaced grid spanning +/- pmSigma standard
 * deviations.
 */
DiscreteRV discretizeNormal(double mean,
    double sigma,
    uword n_elem,
    uword pmSigma)
{
    DiscreteRV drv;
    drv.support = linspace<vec>(mean - pmSigma * sigma, mean + pmSigma * sigma, n_elem);
    drv.probabilities = discretizeNormal(drv.support, mean, sigma);
    return drv;
}

/*
 * VAR(1) stationary distribution.
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
    m_statMean = solve(eyeRho, m_intercept);

    const mat tmp = inv(eyeRho); // TODO Fail on invert of non-stationary case?
    m_statSigma = tmp * m_sigma * tmp.t();

    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, m_statSigma);
    is_stationary = eigval.min() >= 0.0;
}

vec& VAR::stationaryMean()
{
    if (m_statMean.is_empty())
        findStationary();
    return m_statMean;
}

mat& VAR::stationarySigma()
{
    if (m_statSigma.is_empty())
        findStationary();
    return m_statSigma;
}

bool VAR::stationaryQ()
{
    if (m_statMean.is_empty())
        findStationary();
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
 * Markov Chain
 */
MarkovChain::MarkovChain(const AR process, uword sz, double pmSd)
{
    double pMean = process.stationaryMean();
    double pSd = process.stationarySigma();
    rowvec support = linspace<rowvec>(pMean - pmSd * pSd, pMean + pmSd * pSd, sz);

    mat transition(sz, sz);
    for (uword rIx = 0; rIx < sz; ++rIx) {
        double condMean = process.intercept() + process.rho() * support(rIx);
        vec condPrs = discretizeNormal(support.t(), condMean, process.sigma());
        transition.row(rIx) = condPrs.t();
    }

    MarkovChain(support, transition);
}

const rowvec& MarkovChain::stationary()
{
    if (m_stationary.is_empty()) {
        // TODO
    }

    return m_stationary;
}

void MarkovChain::save(std::string fname) const
{
    // TODO
}

/*
 * Discretized VAR(1) => MarkovChain
 */
void DiscreteVAR::save(std::string fname) const
{
    m_grids.save(hdf5_name(fname, "grids", hdf5_opts::replace));
    m_mc.transition().save(hdf5_name(fname, "transitions", hdf5_opts::replace));
    m_var.intercept().save(hdf5_name(fname, "varIntercept", hdf5_opts::replace));
    m_var.rho().save(hdf5_name(fname, "varRho", hdf5_opts::replace));
    m_var.sigma().save(hdf5_name(fname, "varSigma", hdf5_opts::replace));
}

void DiscreteVAR::impl(bool trimGrids)
{
    // Sizing
    m_size = m_var.size();
    m_flatSize = pow(m_supportSize, m_size);
    m_grids.set_size(m_flatSize, m_size);

    umat m_map(m_flatSize, m_size);
    field<vec> m_orthogGrids(m_size);

    const mat mEye = eye(m_size, m_size);

    const mat AA = m_var.intercept();
    const mat BB = m_var.rho();
    const mat GG = m_var.sigma();

    mat LL, VV;
    vec DD;
    svd(LL, DD, VV, GG);

    // Unconditional expectation and covariance matrix
    const mat Atilde = LL.t() * AA;
    const mat Btilde = LL.t() * BB * LL;
    const mat diagDD = diagmat(DD);

    VAR orthog(Atilde, Btilde, diagDD);
    const vec uncondE = orthog.stationaryMean();
    const vec uncondV = orthog.stationarySigma().diag();

    cout << "Original VAR: " << (m_var.stationaryQ() ? "Stationary" : "Nonstationary") << endl;
    m_var.print();
    cout << endl
         << "Orthog VAR: " << (orthog.stationaryQ() ? "Stationary" : "Nonstationary") << endl;
    orthog.print();

    // Prepare logical grids
    for (uword vIx = 0; vIx < m_size; ++vIx) {
        m_orthogGrids(vIx) = linspace<vec>(uncondE(vIx) - pmSd * sqrt(uncondV(vIx)),
            uncondE(vIx) + pmSd * sqrt(uncondV(vIx)), m_supportSize);
    }

    // Prepare flat grids
    for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx) {
        uword tmp = flatIx;
        for (uword vIx = 0; vIx < m_size; ++vIx) {
            const vec logical = m_orthogGrids(vIx);
            m_map(flatIx, vIx) = tmp % m_supportSize;
            m_grids(flatIx, vIx) = logical(m_map(flatIx, vIx));
            tmp = tmp / m_supportSize;
        }
    }

    // Transition matrix
    mat& tran = m_mc.transition();
    tran.set_size(m_flatSize, m_flatSize);
    for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx) {
        rowvec tmpVal = m_grids.row(flatIx);
        vec condMeans = Atilde + Btilde * tmpVal.t();
        vec condSd = sqrt(DD);
        field<vec> condPrs(m_size);
        for (uword vIx = 0; vIx < m_size; ++vIx)
            condPrs(vIx) = discretizeNormal(m_orthogGrids(vIx), condMeans(vIx), condSd(vIx));

        for (uword flatPrIx = 0; flatPrIx < m_flatSize; ++flatPrIx) {
            double goPr = 1.0;
            for (uword vPrIx = 0; vPrIx < m_size; ++vPrIx) {
                vec tmpDist = condPrs(vPrIx);
                goPr *= tmpDist(m_map(flatPrIx, vPrIx));
            }
            tran(flatIx, flatPrIx) = goPr;
        }
    }

    // Adjust grids
    m_grids = (LL * m_grids.t()).t();
    cout << "Completed discretization." << endl;

    if (trimGrids) {
        cx_vec eigval;
        cx_mat eigvec;
        eig_gen(eigval, eigvec, tran.t());

        uword unitEigIx = (abs(abs(eigval) - 1.0)).index_min();

        vec probs = abs(eigvec.col(unitEigIx));
        probs = probs / sum(probs);

        uvec removePoints;
        uword removeCount = 0;
        for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx)
            if (probs(flatIx) <= minPr) {
                removePoints.resize(removeCount + 1);
                removePoints(removeCount) = flatIx;
                removeCount++;
            }

        uword newFlatSize = m_flatSize - removeCount;
        cout << "Removing " << removeCount << " points out of " << m_flatSize << "." << endl;

        mat newGrids = m_grids;
        newGrids.shed_rows(removePoints);
        mat newTran = tran;
        newTran.shed_rows(removePoints);
        newTran.shed_cols(removePoints);
        for (uword fIx = 0; fIx < newFlatSize; ++fIx)
            newTran.row(fIx) = newTran.row(fIx) / sum(newTran.row(fIx));

        /*
        uword advIx = 0;
        for (uword flatIx = 0; flatIx < m_flatSize; ++flatIx)
            if (removePoints(flatIx) == 0) {
                newGrids.row(advIx) = m_grids.row(flatIx);
                rowvec remainingMass(newFlatSize);
                uword advPrIx = 0;
                for (uword flatPrIx = 0; flatPrIx < m_flatSize; ++flatPrIx)
                    if (removePoints(flatPrIx) == 0) {
                        remainingMass(advPrIx) = tran(flatIx, flatPrIx);
                        advPrIx++;
                    }
                newTran.row(advIx) = remainingMass / sum(remainingMass);
                advIx++;
            }
            */

        m_grids = newGrids;
        tran = newTran;
        m_flatSize = newFlatSize;
    }
}

void DiscreteVAR::print() const
{
    cout << "Grid:" << endl;
    m_grids.print();
    cout << endl
         << "Transition matrix: " << endl;
    m_mc.transition().print();
}

} // namespace TimeSeries
