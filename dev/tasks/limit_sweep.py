"""
Limit-experiment calibration of the field interval, K = 2, Sigma = I (unit SE scale).

True field mu = (delta, 0). Observed field W = mu + zeta, zeta ~ N(0, I).
Winner G = argmax W. De-biased error (Theorem 2(ii) limit): Lambda = zeta_G - m(W),
with the closed form m(w) = exp(-(w1-w2)^2/4)/sqrt(pi) for K = 2, Sigma = I.
Re-selection probability of the winner: p_hat(W) = Phi(|W1-W2|/sqrt(2)).

Field: w = W with the winner's entry shrunk, w_G = W_G - s * m(W), s from the rule.
Lambda* = zeta*_{G(w+zeta*)} - m(w+zeta*), zeta* ~ N(0, I).
Second-order error: Lambda - mean(Lambda*).
Two-sided 95% covers iff q025(Lambda*) <= Lambda <= q975(Lambda*); one-sided iff Lambda <= q95(Lambda*).

Rules: A s = 1 (current); B s = 1 - p_hat; C s = (1 - p_hat)^2; D s = 0.5 (constant);
       E "tie-set": if |W1 - W2| <= c (c = 1.5), set both entries to the shrunk winner level (an
         exact tie at the de-biased level), else s = 0 (plug-in) -- the Guo-He-type construction.
Outputs per rule and delta: two-sided and one-sided coverage, mean second-order error, SE ratio.
"""
import numpy as np, pandas as pd, sys
from scipy.stats import norm
rng = np.random.default_rng(7)
R_REP, R_OUT = 4000, 1000
def m2(w):                       # m(w) for K=2, Sigma=I, vectorised on last axis
    d = w[..., 0] - w[..., 1]
    return np.exp(-d * d / 4) / np.sqrt(np.pi)

def field_quantiles(W, G, z, s, tie_c=None, rng=None):
    n=W.shape[0]; w=W.copy(); w[np.arange(n), G] -= s*m2(W)
    if tie_c is not None:
        d=np.abs(W[:,0]-W[:,1]); loser=1-G; tie=d<=tie_c
        w[np.arange(n)[tie], loser[tie]] = w[np.arange(n)[tie], G[tie]]
    zs=rng.standard_normal((n, R_OUT, 2)); V=w[:,None,:]+zs; Gs=V.argmax(2)
    zsG=np.take_along_axis(zs, Gs[:,:,None], axis=2)[:,:,0]; Ls=zsG-m2(V)
    q=np.quantile(Ls,[0.025,0.95,0.975],axis=1); return q[0],q[1],q[2],Ls.mean(1),Ls.std(1)

def sweep(delta, rule, c_tie=1.5):
    mu = np.array([delta, 0.0])
    z = rng.standard_normal((R_REP, 2)); W = mu + z
    G = W.argmax(1); zG = z[np.arange(R_REP), G]
    Lam = zG - m2(W)
    d = np.abs(W[:, 0] - W[:, 1]); phat = norm.cdf(d / np.sqrt(2))
    if rule == 'F':   # envelope over {plug-in, shrunk, tie-set}
        qs=[field_quantiles(W,G,z,sv,tc,rng) for sv,tc in [(0.0,None),(1.0,None),(1.0,1.5)]]
        q025=np.min([q[0] for q in qs],axis=0); q95=np.max([q[1] for q in qs],axis=0); q975=np.max([q[2] for q in qs],axis=0)
        lam_bar=qs[1][3]; err2=Lam-lam_bar
        return dict(delta=delta, rule=rule, cov2=np.mean((q025<=Lam)&(Lam<=q975)), cov1=np.mean(Lam<=q95),
                    bias_est2=err2.mean(), sd_est2=err2.std(), se_ratio=((q975-q025)/(2*1.96)).mean()/err2.std(), s_mean=np.nan, phat_mean=phat.mean())
    if rule == 'A': s = np.ones(R_REP)
    elif rule == 'B': s = 1 - phat
    elif rule == 'C': s = (1 - phat) ** 2
    elif rule == 'D': s = np.full(R_REP, 0.5)
    elif rule == 'E': s = np.where(d <= c_tie, 1.0, 0.0)
    w = W.copy(); w[np.arange(R_REP), G] -= s * m2(W)
    if rule == 'E':                           # tie-set: pull the loser up to the shrunk winner level
        loser = 1 - G
        tie = d <= c_tie
        w[np.arange(R_REP)[tie], loser[tie]] = w[np.arange(R_REP)[tie], G[tie]]
    zs = rng.standard_normal((R_REP, R_OUT, 2))
    V = w[:, None, :] + zs
    Gs = V.argmax(2)
    zsG = np.take_along_axis(zs, Gs[:, :, None], axis=2)[:, :, 0]
    Ls = zsG - m2(V)                          # (R_REP, R_OUT)
    q025, q95, q975 = np.quantile(Ls, [0.025, 0.95, 0.975], axis=1)
    lam_bar = Ls.mean(1); lam_sd = Ls.std(1)
    cov2 = np.mean((q025 <= Lam) & (Lam <= q975)); cov1 = np.mean(Lam <= q95)
    err2 = Lam - lam_bar
    return dict(delta=delta, rule=rule, cov2=cov2, cov1=cov1, bias_est2=err2.mean(),
                sd_est2=err2.std(), se_ratio=lam_sd.mean() / err2.std(), s_mean=s.mean(), phat_mean=phat.mean())

if __name__ == '__main__':
    deltas = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4]
    rules = sys.argv[1].split(',') if len(sys.argv) > 1 else list('ABCDE')
    rows = [sweep(dl, r) for r in rules for dl in deltas]
    res = pd.DataFrame(rows)
    res.to_csv('/home/claude/limit_sweep_K2.csv', index=False)
    pd.set_option('display.width', 200)
    piv = res.pivot(index='delta', columns='rule', values='cov2'); print("two-sided coverage\n", piv.round(3))
    piv1 = res.pivot(index='delta', columns='rule', values='cov1'); print("\none-sided coverage\n", piv1.round(3))
    pb = res.pivot(index='delta', columns='rule', values='bias_est2'); print("\nbias of est2 (SE units)\n", pb.round(3))
    pr = res.pivot(index='delta', columns='rule', values='se_ratio'); print("\nlambda-sd / SD(est2)\n", pr.round(2))
    print("\nworst-case over delta: two-sided", piv.min().round(3).to_dict(), " one-sided", piv1.min().round(3).to_dict())
