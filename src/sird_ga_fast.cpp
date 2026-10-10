// Fast SIRD-GA engine: transition quantities and Kalman log-likelihoods.
//
// Over each observation interval [t_{k-1}, t_k] the ODE solution x(t), the
// resolvent Phi(t, t_{k-1}) and M(t) = int_{t_{k-1}}^t Phi(t,u) Sigma(x(u))
// Phi(t,u)' du are integrated jointly by classical RK4. M satisfies the
// Lyapunov equation dM/dt = J M + M J' + Sigma, with M(t_{k-1}) = 0, so that
// A_k = Phi(t_k, t_{k-1}), F_k = x(t_k) - A_k x(t_{k-1}) and Q_k = M(t_k) / N.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

namespace {

struct Rates { double lambda, gamma, alpha; };

inline arma::vec drift(const arma::vec& x, const Rates& r) {
  double s = x(0), i = x(1);
  arma::vec b(3);
  b(0) = -r.lambda * s * i;
  b(1) = r.lambda * s * i - (r.gamma + r.alpha) * i;
  b(2) = r.gamma * i;
  return b;
}

inline arma::mat jacobian(const arma::vec& x, const Rates& r) {
  double s = x(0), i = x(1);
  arma::mat J(3, 3, arma::fill::zeros);
  J(0, 0) = -r.lambda * i;  J(0, 1) = -r.lambda * s;
  J(1, 0) = r.lambda * i;   J(1, 1) = r.lambda * s - r.gamma - r.alpha;
  J(2, 1) = r.gamma;
  return J;
}

inline arma::mat diffusion(const arma::vec& x, const Rates& r) {
  double s = x(0), i = x(1);
  double li = r.lambda * s * i;
  arma::mat S(3, 3, arma::fill::zeros);
  S(0, 0) = li;   S(0, 1) = -li;
  S(1, 0) = -li;  S(1, 1) = li + (r.gamma + r.alpha) * i;  S(1, 2) = -r.gamma * i;
  S(2, 1) = -r.gamma * i;  S(2, 2) = r.gamma * i;
  return S;
}

struct State { arma::vec x; arma::mat P; arma::mat M; };

inline State deriv(const State& z, const Rates& r) {
  State d;
  arma::mat J = jacobian(z.x, r);
  d.x = drift(z.x, r);
  d.P = J * z.P;
  d.M = J * z.M + z.M * J.t() + diffusion(z.x, r);
  return d;
}

inline State axpy(const State& z, const State& d, double h) {
  State o;
  o.x = z.x + h * d.x; o.P = z.P + h * d.P; o.M = z.M + h * d.M;
  return o;
}

struct Transition {
  arma::mat x;          // (n+1) x 3
  arma::cube A, M;      // 3 x 3 x n
  arma::mat F;          // 3 x n
  bool ok;
};

Transition transition(int n, double Delta, const arma::vec& x0, const Rates& r0,
                      int nsub, const NumericVector& lam_vec = NumericVector()) {
  Transition tr;
  tr.x.set_size(n + 1, 3);
  tr.A.set_size(3, 3, n);
  tr.M.set_size(3, 3, n);
  tr.F.set_size(3, n);
  tr.ok = true;
  tr.x.row(0) = x0.t();
  arma::vec xc = x0;
  double h = Delta / nsub;
  for (int k = 0; k < n; ++k) {
    Rates r = r0;
    if (lam_vec.size() == n) r.lambda = lam_vec[k];
    State z;
    z.x = xc; z.P = arma::eye(3, 3); z.M = arma::zeros(3, 3);
    for (int j = 0; j < nsub; ++j) {
      State k1 = deriv(z, r);
      State k2 = deriv(axpy(z, k1, h / 2), r);
      State k3 = deriv(axpy(z, k2, h / 2), r);
      State k4 = deriv(axpy(z, k3, h), r);
      z.x += h / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
      z.P += h / 6 * (k1.P + 2 * k2.P + 2 * k3.P + k4.P);
      z.M += h / 6 * (k1.M + 2 * k2.M + 2 * k3.M + k4.M);
    }
    if (!z.x.is_finite() || !z.P.is_finite() || !z.M.is_finite()) {
      tr.ok = false;
      return tr;
    }
    tr.A.slice(k) = z.P;
    tr.M.slice(k) = 0.5 * (z.M + z.M.t());
    tr.F.col(k) = z.x - z.P * xc;
    tr.x.row(k + 1) = z.x.t();
    xc = z.x;
  }
  return tr;
}

}  // namespace

//' @noRd
// [[Rcpp::export]]
List sird_ga_transition_cpp(int n, double Delta, NumericVector x0,
                            double lambda, double gamma, double alpha,
                            int nsub = 10) {
  Rates r{lambda, gamma, alpha};
  arma::vec x0v(x0.begin(), 3);
  Transition tr = transition(n, Delta, x0v, r, nsub);
  return List::create(_["ok"] = tr.ok, _["x"] = tr.x, _["A"] = tr.A,
                      _["F"] = tr.F, _["M"] = tr.M);
}

// Prevalence observation (Narci et al. 2021 design): y_k / N = p i_k + V_k,
// k = 0, ..., n, with Var(V_k) = (p(1-p) + tau2) i(t_k) / N.
//' @noRd
// [[Rcpp::export]]
List sird_ga_loglik_prev_cpp(NumericVector y, double N, double Delta,
                             NumericVector x0, double lambda, double gamma,
                             double alpha, double p, double tau2,
                             int nsub = 10, double Q0 = 1e-8) {
  int n = y.size() - 1;
  Rates r{lambda, gamma, alpha};
  arma::vec x0v(x0.begin(), 3);
  Transition tr = transition(n, Delta, x0v, r, nsub);
  NumericVector fitted(n + 1, NA_REAL);
  if (!tr.ok) return List::create(_["loglik"] = R_NegInf, _["fitted"] = fitted);
  arma::rowvec B = {0.0, p, 0.0};
  arma::vec X = x0v;
  arma::mat S = Q0 * arma::eye(3, 3);
  double ll = 0.0;
  for (int k = 0; k <= n; ++k) {
    if (k > 0) {
      X = tr.F.col(k - 1) + tr.A.slice(k - 1) * X;
      S = tr.A.slice(k - 1) * S * tr.A.slice(k - 1).t() + tr.M.slice(k - 1) / N;
    }
    double ik = tr.x(k, 1);
    if (!(ik > 0)) return List::create(_["loglik"] = R_NegInf, _["fitted"] = fitted);
    double Rk = (p * (1 - p) + tau2) * ik / N;
    double m = arma::as_scalar(B * X);
    double G = arma::as_scalar(B * S * B.t()) + Rk;
    if (!(G > 0) || !std::isfinite(G))
      return List::create(_["loglik"] = R_NegInf, _["fitted"] = fitted);
    double nu = y[k] / N - m;
    fitted[k] = N * m;
    ll += -0.5 * (std::log(2 * M_PI) + std::log(G) + nu * nu / G);
    arma::vec K = S * B.t() / G;
    X = X + K * nu;
    S = S - K * B * S;
  }
  return List::create(_["loglik"] = ll, _["fitted"] = fitted);
}

// Incidence observation. The augmented state Z_k = (X_k, X_{k-1}) makes the
// daily increments linear in Z_k:
//   new infections  N (s_{k-1} - s_k),
//   new deaths      N (d_k - d_{k-1}),  d = 1 - s - i - r.
// Reported cases  C_k ~ Bin(new infections, p) + N(0, tau2 * new infections),
// reported deaths D_k ~ Bin(new deaths, q)     + N(0, tau2_d * new deaths),
// both Gaussian-approximated with variances evaluated on the ODE increments.
// `use` selects the streams: 1 = cases only, 2 = deaths only, 3 = both.
//' @noRd
// [[Rcpp::export]]
List sird_ga_loglik_inc_cpp(NumericVector cases, NumericVector deaths,
                            double N, double Delta, NumericVector x0,
                            double lambda, double gamma, double alpha,
                            double p, double tau2, double q, double tau2_d,
                            int use = 3, int nsub = 10, double Q0 = 1e-10,
                            double var_floor = 1e-3,
                            NumericVector lambda_t = NumericVector()) {
  int n = std::max(cases.size(), deaths.size());
  Rates r{lambda, gamma, alpha};
  arma::vec x0v(x0.begin(), 3);
  Transition tr = transition(n, Delta, x0v, r, nsub, lambda_t);
  NumericVector fit_c(n, NA_REAL), fit_d(n, NA_REAL), sd_c(n, NA_REAL),
      sd_d(n, NA_REAL);
  List bad = List::create(_["loglik"] = R_NegInf, _["fitted_cases"] = fit_c,
                          _["fitted_deaths"] = fit_d);
  if (!tr.ok) return bad;

  arma::rowvec hc = {-p, 0, 0, p, 0, 0};
  arma::rowvec hd = {-q, -q, -q, q, q, q};

  arma::vec Z(6);
  Z.head(3) = x0v; Z.tail(3) = x0v;
  arma::mat S = Q0 * arma::eye(6, 6);
  double ll = 0.0;
  int nobs = 0;
  for (int k = 1; k <= n; ++k) {
    const arma::mat& A = tr.A.slice(k - 1);
    arma::mat T(6, 6, arma::fill::zeros);
    T.submat(0, 0, 2, 2) = A;
    T.submat(3, 0, 5, 2) = arma::eye(3, 3);
    arma::vec c(6, arma::fill::zeros);
    c.head(3) = tr.F.col(k - 1);
    arma::mat W(6, 6, arma::fill::zeros);
    W.submat(0, 0, 2, 2) = tr.M.slice(k - 1) / N;
    Z = c + T * Z;
    S = T * S * T.t() + W;

    double inc_c = tr.x(k - 1, 0) - tr.x(k, 0);
    double d_prev = 1 - arma::accu(tr.x.row(k - 1));
    double d_now = 1 - arma::accu(tr.x.row(k));
    double inc_d = d_now - d_prev;
    double floor_frac = var_floor / N;

    std::vector<int> idx;
    std::vector<double> yv, rv;
    std::vector<arma::rowvec> hv;
    bool have_c = (use & 1) && k - 1 < cases.size() && !NumericVector::is_na(cases[k - 1]);
    bool have_d = (use & 2) && k - 1 < deaths.size() && !NumericVector::is_na(deaths[k - 1]);
    fit_c[k - 1] = N * arma::as_scalar(hc * Z);
    fit_d[k - 1] = N * arma::as_scalar(hd * Z);
    sd_c[k - 1] = N * std::sqrt(std::max(arma::as_scalar(hc * S * hc.t()) +
                                (p * (1 - p) + tau2) * std::max(inc_c, floor_frac) / N, 0.0));
    sd_d[k - 1] = N * std::sqrt(std::max(arma::as_scalar(hd * S * hd.t()) +
                                (q * (1 - q) + tau2_d) * std::max(inc_d, floor_frac) / N, 0.0));
    if (have_c) {
      hv.push_back(hc); yv.push_back(cases[k - 1] / N);
      rv.push_back((p * (1 - p) + tau2) * std::max(inc_c, floor_frac) / N);
    }
    if (have_d) {
      hv.push_back(hd); yv.push_back(deaths[k - 1] / N);
      rv.push_back((q * (1 - q) + tau2_d) * std::max(inc_d, floor_frac) / N);
    }
    int m = hv.size();
    if (m == 0) continue;
    arma::mat H(m, 6);
    arma::vec yy(m), Rd(m);
    for (int j = 0; j < m; ++j) { H.row(j) = hv[j]; yy(j) = yv[j]; Rd(j) = rv[j]; }
    arma::vec nu = yy - H * Z;
    arma::mat G = H * S * H.t() + arma::diagmat(Rd);
    G = 0.5 * (G + G.t());
    double logdet, sgn;
    arma::log_det(logdet, sgn, G);
    if (!(sgn > 0) || !std::isfinite(logdet)) return bad;
    arma::mat Gi;
    if (!arma::inv_sympd(Gi, G)) return bad;  // ill-conditioned G: penalise, do not throw
    ll += -0.5 * (m * std::log(2 * M_PI) + logdet + arma::as_scalar(nu.t() * Gi * nu));
    arma::mat K = S * H.t() * Gi;
    Z = Z + K * nu;
    arma::mat IKH = arma::eye(6, 6) - K * H;
    S = IKH * S * IKH.t() + K * arma::diagmat(Rd) * K.t();  // Joseph form
    nobs += m;
  }
  return List::create(_["loglik"] = ll, _["nobs"] = nobs,
                      _["fitted_cases"] = fit_c, _["fitted_deaths"] = fit_d,
                      _["sd_cases"] = sd_c, _["sd_deaths"] = sd_d,
                      _["x"] = tr.x);
}
