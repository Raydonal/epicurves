// Generic Gaussian-approximation engine for mass-action compartmental models.
//
// State: d compartment fractions plus C event counters (cumulative number of
// firings of the observed events, divided by N). Event e moves the state by
// column e of the stoichiometry matrix V (counters included) at rate
//   theta[pidx[e]] * x[r1[e]] * x[r2[e]]   (r2[e] < 0: no second factor;
//                                            r1[e] < 0: constant rate).
// Over each observation interval the ODE solution, the resolvent and the
// Lyapunov integral M(t) = int Phi(t,u) Sigma(x(u)) Phi(t,u)' du are
// integrated jointly by classical RK4, exactly as in sird_ga_fast.cpp:
//   dx/dt = b(x),  dPhi/dt = J Phi,  dM/dt = J M + M J' + Sigma.
// Observations at t_1..t_n are linear in Z_k = (X_k, X_{k-1}):
//   incidence of event e: p (c_k - c_{k-1}),   prevalence of x_j: p x_{j,k},
// with variance (p(1-p) + tau2) m_k / N, m_k the ODE value of the latent
// count fraction (floored at var_floor / N).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

namespace {

struct Spec {
  arma::mat V;               // D x E stoichiometry (compartments + counters)
  std::vector<int> pidx, r1, r2;
  int D, E;
};

inline arma::vec rates(const arma::vec& x, const arma::rowvec& th, const Spec& s) {
  arma::vec r(s.E);
  for (int e = 0; e < s.E; ++e) {
    double v = th(s.pidx[e]);
    if (s.r1[e] >= 0) v *= x(s.r1[e]);
    if (s.r2[e] >= 0) v *= x(s.r2[e]);
    r(e) = v;
  }
  return r;
}

// E x D matrix of partial derivatives of the rates
inline arma::mat drates(const arma::vec& x, const arma::rowvec& th, const Spec& s) {
  arma::mat G(s.E, s.D, arma::fill::zeros);
  for (int e = 0; e < s.E; ++e) {
    double k = th(s.pidx[e]);
    int a = s.r1[e], b = s.r2[e];
    if (a < 0) continue;
    if (b < 0) { G(e, a) += k; continue; }
    G(e, a) += k * x(b);
    G(e, b) += k * x(a);
  }
  return G;
}

// Workspace for one RK4 step; avoids allocating temporaries in the inner loop.
struct Work {
  arma::vec r, x1, kx[4];
  arma::mat G, J, P1, M1, kP[4], kM[4], VJ;
  Work(int D, int E) : r(E), x1(D), G(E, D), J(D, D), P1(D, D), M1(D, D), VJ(D, E) {
    for (int i = 0; i < 4; ++i) { kx[i].set_size(D); kP[i].set_size(D, D); kM[i].set_size(D, D); }
  }
};

// derivative of (x, P, M) at (x, P, M), written into slot i of the workspace
inline void deriv(const arma::vec& x, const arma::mat& P, const arma::mat& M,
                  const arma::rowvec& th, const Spec& s, Work& w, int i) {
  w.r = rates(x, th, s);
  w.G = drates(x, th, s);
  w.J = s.V * w.G;
  w.kx[i] = s.V * w.r;
  w.kP[i] = w.J * P;
  w.VJ = s.V.each_row() % w.r.t();
  w.kM[i] = w.J * M;
  w.kM[i] += w.kM[i].t();
  w.kM[i] += w.VJ * s.V.t();
}

struct Transition {
  arma::mat x;      // (n+1) x D
  arma::cube A, M;  // D x D x n
  arma::mat F;      // D x n
  bool ok;
};

Transition transition(int n, double Delta, const arma::vec& x0, const arma::mat& theta,
                      const Spec& s, int nsub) {
  Transition tr;
  int D = s.D;
  tr.x.set_size(n + 1, D); tr.A.set_size(D, D, n); tr.M.set_size(D, D, n); tr.F.set_size(D, n);
  tr.ok = true;
  tr.x.row(0) = x0.t();
  arma::vec xc = x0, x(D);
  arma::mat P(D, D), M(D, D);
  Work w(D, s.E);
  double h = Delta / nsub;
  for (int k = 0; k < n; ++k) {
    arma::rowvec th = theta.row(theta.n_rows == 1 ? 0 : k);
    x = xc; P.eye(); M.zeros();
    for (int j = 0; j < nsub; ++j) {
      deriv(x, P, M, th, s, w, 0);
      deriv(x + (h / 2) * w.kx[0], P + (h / 2) * w.kP[0], M + (h / 2) * w.kM[0], th, s, w, 1);
      deriv(x + (h / 2) * w.kx[1], P + (h / 2) * w.kP[1], M + (h / 2) * w.kM[1], th, s, w, 2);
      deriv(x + h * w.kx[2], P + h * w.kP[2], M + h * w.kM[2], th, s, w, 3);
      x += (h / 6) * (w.kx[0] + 2 * w.kx[1] + 2 * w.kx[2] + w.kx[3]);
      P += (h / 6) * (w.kP[0] + 2 * w.kP[1] + 2 * w.kP[2] + w.kP[3]);
      M += (h / 6) * (w.kM[0] + 2 * w.kM[1] + 2 * w.kM[2] + w.kM[3]);
    }
    if (!x.is_finite() || !P.is_finite() || !M.is_finite()) { tr.ok = false; return tr; }
    tr.A.slice(k) = P;
    tr.M.slice(k) = 0.5 * (M + M.t());
    tr.F.col(k) = x - P * xc;
    tr.x.row(k + 1) = x.t();
    xc = x;
  }
  return tr;
}

Spec make_spec(const NumericMatrix& V, const IntegerVector& pidx, const IntegerVector& r1,
               const IntegerVector& r2) {
  Spec s;
  s.V = arma::mat(const_cast<double*>(V.begin()), V.nrow(), V.ncol(), true);
  s.D = V.nrow(); s.E = V.ncol();
  s.pidx.assign(pidx.begin(), pidx.end());
  s.r1.assign(r1.begin(), r1.end());
  s.r2.assign(r2.begin(), r2.end());
  return s;
}

}  // namespace

//' @noRd
// [[Rcpp::export]]
List ga_transition_cpp(int n, double Delta, NumericVector x0, NumericMatrix theta,
                       NumericMatrix V, IntegerVector pidx, IntegerVector r1, IntegerVector r2,
                       int nsub = 10) {
  Spec s = make_spec(V, pidx, r1, r2);
  arma::vec x0v(x0.begin(), x0.size());
  arma::mat th(theta.begin(), theta.nrow(), theta.ncol());
  Transition tr = transition(n, Delta, x0v, th, s, nsub);
  return List::create(_["ok"] = tr.ok, _["x"] = tr.x, _["A"] = tr.A, _["F"] = tr.F, _["M"] = tr.M);
}

// Y: n x m matrix of observed counts (NA allowed). For stream j:
//   obs_type[j] = 1 (incidence of counter state obs_index[j]) or 2 (prevalence of
//   compartment obs_index[j]); obs_p[j], obs_tau2[j] its reporting probability and
//   extra-binomial variance scale.
//' @noRd
// [[Rcpp::export]]
List ga_loglik_cpp(NumericMatrix Y, double N, double Delta, NumericVector x0,
                   NumericMatrix theta, NumericMatrix V, IntegerVector pidx,
                   IntegerVector r1, IntegerVector r2, IntegerVector obs_type,
                   IntegerVector obs_index, NumericVector obs_p, NumericVector obs_tau2,
                   int nsub = 10, double Q0 = 1e-10, double var_floor = 1e-3) {
  Spec s = make_spec(V, pidx, r1, r2);
  int n = Y.nrow(), m = Y.ncol(), D = s.D;
  arma::vec x0v(x0.begin(), x0.size());
  arma::mat th(theta.begin(), theta.nrow(), theta.ncol());
  NumericMatrix fit(n, m), sdv(n, m);
  std::fill(fit.begin(), fit.end(), NA_REAL);
  std::fill(sdv.begin(), sdv.end(), NA_REAL);
  List bad = List::create(_["loglik"] = R_NegInf, _["fitted"] = fit, _["sd"] = sdv);
  Transition tr = transition(n, Delta, x0v, th, s, nsub);
  if (!tr.ok) return bad;

  std::vector<arma::rowvec> h(m);
  for (int j = 0; j < m; ++j) {
    arma::rowvec hj(2 * D, arma::fill::zeros);
    if (obs_type[j] == 1) { hj(obs_index[j]) = obs_p[j]; hj(D + obs_index[j]) = -obs_p[j]; }
    else hj(obs_index[j]) = obs_p[j];
    h[j] = hj;
  }
  arma::vec Z(2 * D);
  Z.head(D) = x0v; Z.tail(D) = x0v;
  arma::mat S = Q0 * arma::eye(2 * D, 2 * D);
  double ll = 0.0, floor_frac = var_floor / N;
  int nobs = 0;
  for (int k = 1; k <= n; ++k) {
    arma::mat T(2 * D, 2 * D, arma::fill::zeros);
    T.submat(0, 0, D - 1, D - 1) = tr.A.slice(k - 1);
    T.submat(D, 0, 2 * D - 1, D - 1) = arma::eye(D, D);
    arma::vec c(2 * D, arma::fill::zeros);
    c.head(D) = tr.F.col(k - 1);
    arma::mat W(2 * D, 2 * D, arma::fill::zeros);
    W.submat(0, 0, D - 1, D - 1) = tr.M.slice(k - 1) / N;
    Z = c + T * Z;
    S = T * S * T.t() + W;

    std::vector<int> use;
    std::vector<double> rv;
    for (int j = 0; j < m; ++j) {
      int i = obs_index[j];
      double lat = obs_type[j] == 1 ? tr.x(k, i) - tr.x(k - 1, i) : tr.x(k, i);
      double r = (obs_p[j] * (1 - obs_p[j]) + obs_tau2[j]) * std::max(lat, floor_frac) / N;
      fit(k - 1, j) = N * arma::as_scalar(h[j] * Z);
      sdv(k - 1, j) = N * std::sqrt(std::max(arma::as_scalar(h[j] * S * h[j].t()) + r, 0.0));
      if (!NumericMatrix::is_na(Y(k - 1, j))) { use.push_back(j); rv.push_back(r); }
    }
    int q = use.size();
    if (q == 0) continue;
    arma::mat H(q, 2 * D);
    arma::vec yy(q), Rd(q);
    for (int a = 0; a < q; ++a) { H.row(a) = h[use[a]]; yy(a) = Y(k - 1, use[a]) / N; Rd(a) = rv[a]; }
    arma::vec nu = yy - H * Z;
    arma::mat G = H * S * H.t() + arma::diagmat(Rd);
    G = 0.5 * (G + G.t());
    double logdet, sgn;
    arma::log_det(logdet, sgn, G);
    if (!(sgn > 0) || !std::isfinite(logdet)) return bad;
    arma::mat Gi;
    if (!arma::inv_sympd(Gi, G)) return bad;
    ll += -0.5 * (q * std::log(2 * M_PI) + logdet + arma::as_scalar(nu.t() * Gi * nu));
    arma::mat K = S * H.t() * Gi;
    Z = Z + K * nu;
    arma::mat IKH = arma::eye(2 * D, 2 * D) - K * H;
    S = IKH * S * IKH.t() + K * arma::diagmat(Rd) * K.t();  // Joseph form
    nobs += q;
  }
  return List::create(_["loglik"] = ll, _["nobs"] = nobs, _["fitted"] = fit, _["sd"] = sdv,
                      _["x"] = tr.x);
}
