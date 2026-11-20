# Author: Giacomo Bignardi
# Partially adapted from the BrainSpace tutorial (https://brainspace.readthedocs.io)
# Description: decompose individual functional connectomes (Schaefer 400) into 10 gradients by
#              diffusion map embedding and align them to the group-level template by Procrustes rotation.
# Revision:
#   - template computed from the cortical 400 x 400 block of the group-average connectome,
#     as the individual matrices (previously 400 cortical rows x all 419 columns);
#   - standard BrainSpace instead of a locally modified copy;
#   - template sign fixed by convention (G1: default mode > somatomotor; G2: visual > somatomotor);
#   - additionally writes sign-matched, unrotated G1 (sensitivity analysis) and per-subject
#     alignment diagnostics.
# To run reticulate::source_python("Python/03_FC_G.py")
import os
import glob
import numpy as np
import pandas as pd
from brainspace.gradient import GradientMaps
from brainspace.gradient.alignment import procrustes

# directories (data and outputs live outside the repository)
wd_up    = os.path.dirname(os.getcwd())
wd_out   = os.path.join(wd_up, "11_noa_output")
wd_noa   = os.path.join(wd_up, "03_rawData", "noa")
wd_fc    = os.path.join(wd_up, "03_rawData", "private", "rsfMRIavg")
ids_file = os.path.join(wd_noa, "subj_inclusion.csv")

GM_KW = dict(n_components=10, kernel="normalized_angle", approach="dm", random_state=42)
N = 400  # cortical parcels

# GROUP LEVEL ####
fc_mean = pd.read_csv(os.path.join(wd_out, "02_FC_m.csv")).values[:N, :N]
gm_t = GradientMaps(**GM_KW)
gm_t.fit(fc_mean)
T = gm_t.gradients_.copy()

# fix the arbitrary sign of the template by convention (Yeo-Krienen 7 labels)
ann = pd.read_csv(os.path.join(wd_noa, "cortical_types.csv"))
def parcels(net):
    return ann.loc[ann["label_Yeo_7"].str.contains(net), "parcel_Yeo"].values.astype(int) - 1
dmn, som, vis = parcels("_Default_"), parcels("_SomMot_"), parcels("_Vis_")
if T[dmn, 0].mean() < T[som, 0].mean():
    T[:, 0] *= -1                                             # G1: default mode positive
if T[vis, 1].mean() < T[som, 1].mean():
    T[:, 1] *= -1                                             # G2: visual positive
print("variance explained, G1:", round(gm_t.lambdas_[0] / gm_t.lambdas_.sum(), 4),
      "| G1-G3:", round(gm_t.lambdas_[:3].sum() / gm_t.lambdas_.sum(), 4))
pd.DataFrame(T).to_csv(os.path.join(wd_out, "03_FC_m_G.csv"), index=False)

# INDIVIDUAL LEVEL ####
ids = set(pd.read_csv(ids_file, header=None)[0].astype(int))

def individual(fc):
    """Gradients of one connectome: aligned to the template, sign-matched G1, diagnostics."""
    gm = GradientMaps(**GM_KW)
    gm.fit(fc)
    U, lam = gm.gradients_, gm.lambdas_
    A = procrustes(U, T)                                         # rotation onto the fixed template
    r_raw = np.array([np.corrcoef(U[:, k], T[:, 0])[0, 1] for k in range(U.shape[1])])
    S = U[:, :1] * np.sign(r_raw[0])                             # first component, sign-matched only
    diag = dict(expl_g1=lam[0] / lam.sum(),                      # explanatory ratio of G1
                r_raw_g1=abs(r_raw[0]),                          # raw comp 1 vs template G1
                best_raw_comp=int(np.argmax(np.abs(r_raw))) + 1, # 1 = raw comp 1 is the S-A gradient
                r_aligned_g1=np.corrcoef(A[:, 0], T[:, 0])[0, 1],
                resid=np.linalg.norm(A - T) / np.linalg.norm(T)) # normalised Procrustes residual
    return A, S, diag

def run(day_dir, session, suffix):
    files = glob.glob(os.path.join(wd_fc, day_dir, "*.csv"))
    subs = sorted(s for s in (int(os.path.basename(f).rsplit(".", 1)[0]) for f in files) if s in ids)
    out_a, out_s, out_d = [], [], []
    for s in subs:
        fc = pd.read_csv(os.path.join(wd_fc, day_dir, f"{s}.csv"), header=None).values[:N, :N]
        if np.isnan(fc).any():
            A, S, d = np.full((N, 3), np.nan), np.full((N, 1), np.nan), {}
        else:
            A, S, d = individual(fc)
        meta = pd.DataFrame({"Sub": s, "Session": session, "Parcel": np.arange(1, N + 1)})
        out_a.append(pd.concat([meta, pd.DataFrame(A[:, :3])], axis=1))
        out_s.append(pd.concat([meta, pd.DataFrame(S)], axis=1))
        out_d.append({"Sub": s, "Session": session, **d})
    pd.concat(out_a).to_csv(os.path.join(wd_out, f"03_GFC_i{suffix}.csv"), index=False)
    pd.concat(out_s).to_csv(os.path.join(wd_out, f"03_GFC_i{suffix}_signmatched.csv"), index=False)
    pd.DataFrame(out_d).to_csv(os.path.join(wd_out, f"03_GFC_i{suffix}_diag.csv"), index=False)
    print(f"{day_dir}: {len(subs)} subjects written")

run("day1", 1.5, "_d1")      # sessions averaged within day 1
run("day2", 2.5, "_d2")      # sessions averaged within day 2
run("day_avg", 0, "")        # averaged across days
