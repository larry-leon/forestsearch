import numpy as np, pandas as pd, sys, time
rng=np.random.default_rng(3)
R_REP, R_OUT, R_IN = 800, 400, 200
def m_hat(V, L, rng):           # V: (..., K) fields; m by inner MC with shared inner draws
    K=V.shape[-1]; zi=rng.standard_normal((R_IN,K))@L.T
    U=V[...,None,:]+zi; Gi=U.argmax(-1)
    return np.take_along_axis(np.broadcast_to(zi,U.shape), Gi[...,None], axis=-1)[...,0].mean(-1)
def run(K, rho, delta, kappas):
    Sig=np.full((K,K),rho); np.fill_diagonal(Sig,1.0); L=np.linalg.cholesky(Sig)
    mu=np.zeros(K); mu[0]=delta
    z=rng.standard_normal((R_REP,K))@L.T; W=mu+z; G=W.argmax(1); zG=z[np.arange(R_REP),G]
    mW=m_hat(W,L,rng); Lam=zG-mW
    w=W.copy(); w[np.arange(R_REP),G]-=mW                    # rule A: shrink winner by m(W)
    zs=rng.standard_normal((R_REP,R_OUT,K))@L.T; V=w[:,None,:]+zs; Gs=V.argmax(2)
    zsG=np.take_along_axis(zs,Gs[:,:,None],axis=2)[:,:,0]
    Ls=np.empty((R_REP,R_OUT))
    for i in range(0,R_REP,100): Ls[i:i+100]=zsG[i:i+100]-m_hat(V[i:i+100],L,rng)
    q025,q95,q975=np.quantile(Ls,[0.025,0.95,0.975],axis=1); lb=Ls.mean(1)
    out=dict(K=K,rho=rho,delta=delta,cov1=np.mean(Lam<=q95),bias_bt=Lam.mean(),bias_est2=(Lam-lb).mean(),phat=np.mean(Gs==G[:,None]))
    for k in kappas:
        lo=lb+k*(q025-lb); hi=lb+k*(q975-lb); out[f'cov2_k{k}']=np.mean((lo<=Lam)&(Lam<=hi))
    return out
if __name__=='__main__':
    K=int(sys.argv[1]); rho=float(sys.argv[2]); kappas=[1.0,1.2,1.3]
    rows=[]
    for d in [0,1,2,3]:
        t=time.time(); rows.append(run(K,rho,d,kappas)); print(d,'%.0fs'%(time.time()-t),flush=True)
    df=pd.DataFrame(rows); pd.set_option('display.width',200); print(df.round(3).to_string())
    df.to_csv(f'/home/claude/limit_sweep_K{K}_rho{rho}.csv',index=False)
