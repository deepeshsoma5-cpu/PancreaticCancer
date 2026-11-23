# small deepsurv style net: genes -> 16 hidden -> 1 risk score, trained with cox loss
# same idea as the other models: gene filter happens inside each training fold
import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

torch.manual_seed(42)
np.random.seed(42)
nrep = int(os.environ.get("N_REPEATS", "3"))

X = pd.read_csv("data/nn_X.csv", index_col=0)
meta = pd.read_csv("data/nn_meta.csv", index_col=0).loc[X.index]
folds = pd.read_csv("data/nn_folds.csv")
sel = pd.read_csv("data/nn_selected.csv")
t_all = meta["time_months"].values
e_all = meta["status"].values
print("patients:", len(X), " genes:", X.shape[1])


def cox_loss(h, t, e):
    # negative cox partial likelihood (ignoring tied times)
    order = torch.argsort(t, descending=True)
    h = h[order]
    e = e[order]
    lse = torch.logcumsumexp(h, dim=0)
    return -((h - lse) * e).sum() / e.sum().clamp(min=1.0)


def cindex(risk, t, e):
    num = 0.0
    den = 0.0
    for i in range(len(t)):
        if e[i] == 1:
            for j in range(len(t)):
                if t[j] > t[i]:
                    den += 1
                    if risk[i] > risk[j]:
                        num += 1
                    elif risk[i] == risk[j]:
                        num += 0.5
    return num / den


def fit_predict(xtr, ttr, etr, xte):
    # hold out 20% of the training fold for early stopping
    n = len(ttr)
    idx = np.random.permutation(n)
    nval = int(0.2 * n)
    va, tr = idx[:nval], idx[nval:]
    f = lambda a: torch.tensor(a, dtype=torch.float32)
    x_tr, t_tr, e_tr = f(xtr[tr]), f(ttr[tr]), f(etr[tr])
    x_va, t_va, e_va = f(xtr[va]), f(ttr[va]), f(etr[va])

    model = nn.Sequential(nn.Linear(xtr.shape[1], 16), nn.ReLU(),
                          nn.Dropout(0.3), nn.Linear(16, 1))
    opt = torch.optim.Adam(model.parameters(), lr=0.005, weight_decay=0.01)

    best = 1e9
    best_state = None
    wait = 0
    for ep in range(300):
        model.train()
        opt.zero_grad()
        loss = cox_loss(model(x_tr).squeeze(1), t_tr, e_tr)
        loss.backward()
        opt.step()

        model.eval()
        with torch.no_grad():
            vl = cox_loss(model(x_va).squeeze(1), t_va, e_va).item()
        if vl < best - 1e-4:
            best = vl
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            wait = 0
        else:
            wait += 1
            if wait >= 25:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        return model(f(xte)).squeeze(1).numpy()


scores = []
for r in range(1, nrep + 1):
    fr = folds[folds.rep == r].set_index("patient")["fold"].loc[X.index].values
    for k in range(1, 6):
        genes = sel[(sel.rep == r) & (sel.fold == k)]["gene"].tolist()
        tr = np.where(fr != k)[0]
        te = np.where(fr == k)[0]
        xtr = X.iloc[tr][genes].values
        xte = X.iloc[te][genes].values
        risk = fit_predict(xtr, t_all[tr], e_all[tr], xte)
        scores.append(cindex(risk, t_all[te], e_all[te]))
    print("done with repeat", r, "of", nrep)

res = pd.read_csv("results/07_cv_results.csv")
print("\n===== CROSS-VALIDATED C-INDEX (TCGA only) =====")
for m, g in res.groupby("model"):
    print(f"{m:<22} mean {g.cindex.mean():.3f} | sd {g.cindex.std():.3f}")
print(f"{'NeuralNetwork':<22} mean {np.mean(scores):.3f} | sd {np.std(scores, ddof=1):.3f}")
pd.DataFrame({"model": "NeuralNetwork", "cindex": scores}).to_csv("results/09_nn_cv_results.csv", index=False)
