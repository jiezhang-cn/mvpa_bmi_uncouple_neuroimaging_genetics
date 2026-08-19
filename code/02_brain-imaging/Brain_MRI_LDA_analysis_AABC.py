"""
AABC replication of UKB Brain_MRI_LDA_analysis.ipynb (sMRI track only).

Inputs
------
- Cortical_Areal_Volumes.csv (HCP MMP 1.0, 360 cortical ROIs)
- asegstats.csv             (FreeSurfer subcortical volumes; we keep 17 true
                             subcortical ROIs per hemisphere where applicable,
                             matching UKB's Harvard-Oxford subcortical scope:
                             Hippocampus, Amygdala, Caudate, Putamen, Pallidum,
                             Thalamus, Accumbens, VentralDC, Cerebellum-Cortex,
                             Cerebellum-WM, Brain-Stem)
- aabc_pa_obesity_data.csv  (already has high/low_mvpa_bmi_uncouple at
                             BUFFER=0.25, age 40-69, sex_num, race_bin,
                             education)

Pipeline (mirrors the UKB notebook)
-----------------------------------
1. Z-score features, then deconfound with eTIV + site dummies (UKB's
   head-motion / head-position fields are not collected in AABC, so we
   substitute the available QC metadata).
2. Bootstrap LDA, n_BS = 500. Each iteration:
     - resample minority class with replacement to balance the two classes,
     - within the bootstrap fold, regress out age + sex + race_bin +
       education from the brain features (matches the UKB notebook's
       in-fold covariate cleaning),
     - fit LinearDiscriminantAnalysis, store coef_.
3. Mean coef + 5/95-percentile CIs across the 500 bootstraps. Features
   whose CI does not span 0 are flagged as significant.
4. Save coefficient table per phenotype.
"""
import os
import warnings
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.utils import resample
from nilearn.signal import clean

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---- paths ----
AABC_IDP_DIR = r"C:/Users/张杰/Desktop/Cohort Data/AABC/AABC_Release2_StructuralIDPs"
AABC_NIMG    = r"C:/Users/张杰/Desktop/Cohort Data/AABC/AABC_Release2_Non-imaging_Data-XL.csv"
PHENO_CSV    = r"C:/Users/张杰/Desktop/PA_obesity_biological/data/aabc_pa_obesity_data.csv"
OUT_DIR      = r"C:/Users/张杰/Desktop/PA_obesity_biological/results/AABC_replication"
os.makedirs(OUT_DIR, exist_ok=True)

# ---- 1. Load cortical volumes (HCP MMP 1.0, 360 ROIs) ----
cort = pd.read_csv(os.path.join(AABC_IDP_DIR, "Cortical_Areal_Volumes.csv"))
cort["id"]    = cort["x___"].str.split("_").str[0]
cort["event"] = cort["x___"].str.split("_").str[1]
cort = cort[cort["event"] == "V1"].drop(columns=["x___", "event"]).set_index("id")
cort.columns = [f"cort_{c}" for c in cort.columns]
cort = cort.dropna()
print(f"Cortical volumes (V1): {cort.shape}")

# ---- 2. Load FreeSurfer aseg subcortical volumes ----
aseg = pd.read_csv(os.path.join(AABC_IDP_DIR, "asegstats.csv"))
aseg["id"]    = aseg["Session"].str.split("_").str[0]
aseg["event"] = aseg["Session"].str.split("_").str[1]
aseg = aseg[aseg["event"] == "V1"].drop(columns=["Session", "event"]).set_index("id")

# Keep only true subcortical/cerebellum ROI volumes; drop totals, eTIV-related
# meta, surface holes, hypointensities, ventricles, choroid plexus, vessels.
SUBCORT_KEEP = [
    "FS_L_Cerebellum-White-Matter", "FS_R_Cerebellum-White-Matter",
    "FS_L_Cerebellum-Cortex",       "FS_R_Cerebellum-Cortex",
    "FS_L_Thalamus-Proper",         "FS_R_Thalamus-Proper",
    "FS_L_Caudate",                 "FS_R_Caudate",
    "FS_L_Putamen",                 "FS_R_Putamen",
    "FS_L_Pallidum",                "FS_R_Pallidum",
    "FS_L_Hippocampus",             "FS_R_Hippocampus",
    "FS_L_Amygdala",                "FS_R_Amygdala",
    "FS_L_Accumbens-area",          "FS_R_Accumbens-area",
    "FS_L_VentralDC",               "FS_R_VentralDC",
    "Brain-Stem",
]
META_FOR_DECONF = ["EstimatedTotalIntraCranialVol"]
keep_cols = SUBCORT_KEEP + META_FOR_DECONF
aseg = aseg[keep_cols].dropna()
aseg.columns = [f"sub_{c}" for c in aseg.columns]
print(f"Subcortical (aseg) at V1: {aseg.shape}")

# Pull eTIV out as a deconfound covariate, drop it from features
etiv = aseg["sub_EstimatedTotalIntraCranialVol"].copy()
aseg = aseg.drop(columns=["sub_EstimatedTotalIntraCranialVol"])

# ---- 3. Site (recruitment site) for deconfound ----
ni = pd.read_csv(AABC_NIMG, low_memory=False)
ni.columns = [c.split(" - ")[0] for c in ni.columns]
ni = ni[ni["id"] != "id"]
ni_v1 = ni[ni["event"] == "V1"][["id", "site"]].drop_duplicates(subset="id").set_index("id")

# ---- 4. Merge sMRI matrix ----
sMRI = cort.join(aseg, how="inner")
sMRI_name = sMRI.columns.tolist()
print(f"sMRI feature matrix (V1): {sMRI.shape}")

# ---- 5. Z-score + deconfound (eTIV + site dummies) ----
common = sMRI.index.intersection(etiv.index).intersection(ni_v1.index)
sMRI = sMRI.loc[common]
etiv = etiv.loc[common]
sites = pd.get_dummies(ni_v1.loc[common, "site"], drop_first=False).astype(float).values

sc = StandardScaler()
X_z = sc.fit_transform(sMRI.values)
conf_mat = np.hstack([np.atleast_2d(etiv.values).T, sites])
X_dec = clean(X_z, confounds=conf_mat, detrend=False, standardize=False)
df_input_sMRI = pd.DataFrame(X_dec, index=sMRI.index, columns=sMRI_name)
df_input_sMRI.to_csv(os.path.join(OUT_DIR, "aabc_sMRI_deconfounded.csv"))
print(f"Deconfounded sMRI saved. Shape: {df_input_sMRI.shape}")

# ---- 6. Phenotype + covariates ----
pheno = pd.read_csv(PHENO_CSV)
pheno = pheno.rename(columns={"id": "sub_id"}).set_index("sub_id")
edu_cat = pd.Categorical(pheno["education"]).codes  # numeric encoding for confounding
pheno["edu_num"] = edu_cat
race_dum = pd.get_dummies(pheno["race_bin"], drop_first=True).astype(float)
race_dum.columns = [f"raceD_{c}" for c in race_dum.columns]
pheno = pd.concat([pheno, race_dum], axis=1)

COV_COLS = ["age_open", "sex_num", "edu_num"] + race_dum.columns.tolist()
print("Covariate columns for in-fold cleaning:", COV_COLS)

# ---- 7. Bootstrap LDA helpers ----
def perform_bootstrap(Y, n_samples):
    df_y_test  = Y.sample(frac=0.1, random_state=0)
    df_y_train = Y[~Y.index.isin(df_y_test.index)]
    minority   = df_y_train.value_counts().min()
    samples = []
    for s in range(n_samples):
        np.random.seed(s)
        c0 = resample(df_y_train[df_y_train == 0].dropna(), n_samples=minority, replace=True)
        c1 = resample(df_y_train[df_y_train == 1].dropna(), n_samples=minority, replace=True)
        samples.append(pd.concat([c0, c1], axis=0).index.to_list())
    return samples, df_y_test.index


def calculate_feature_coef(X, Y, n_BS, covariates=None):
    samples, _ = perform_bootstrap(Y, n_BS)
    # Regularised LDA: lsqr solver with Ledoit-Wolf automatic shrinkage —
    # required because N (~200) << p (381) makes the within-class covariance
    # matrix near-singular under the default solver.
    lda = LDA(solver="lsqr", shrinkage="auto")
    out = np.zeros((n_BS, X.shape[1]))
    for i, idxs in enumerate(samples):
        X_train = X.loc[idxs].values
        y_train = np.ravel(Y.loc[idxs].values)
        if covariates is not None:
            cov_train = covariates.loc[idxs].values
            X_train = clean(X_train, confounds=cov_train, detrend=False, standardize=False)
        out[i, :] = lda.fit(X_train, y_train).coef_
    return out


def calculate_significant_coef(coef_mat, feat_names, alpha):
    coef_mean = coef_mat.mean(axis=0)
    lo = np.percentile(coef_mat, alpha / 2,        axis=0)
    hi = np.percentile(coef_mat, 100 - alpha / 2,  axis=0)
    df_all = pd.DataFrame({"coef": coef_mean, "lower": lo, "upper": hi},
                          index=feat_names)
    sig_mask = (lo > 0) | (hi < 0)
    df_sig = df_all[sig_mask].copy()
    return df_sig, df_all


# ---- 8. Run for HD ----
def run_lda_for(target_col, tag):
    Y = pheno[target_col].dropna()
    common_idx = Y.index.intersection(df_input_sMRI.index)
    Y = Y.loc[common_idx].astype(int)
    X = df_input_sMRI.loc[common_idx]
    cov = pheno.loc[common_idx, COV_COLS].astype(float).fillna(pheno[COV_COLS].astype(float).median())
    print(f"\n=== {tag} ===")
    print(f"  N (overlap with sMRI): {len(common_idx)}")
    print(f"  class counts:\n{Y.value_counts()}")
    n_BS = 500
    alpha = 10
    coef = calculate_feature_coef(X, Y, n_BS, covariates=cov)
    df_sig, df_all = calculate_significant_coef(coef, sMRI_name, alpha)
    out_csv = os.path.join(OUT_DIR, f"LDA_sMRI_{tag}_coef_5-95CI_all.csv")
    df_all.to_csv(out_csv, index_label="feature")
    print(f"  Significant features (CI excludes 0): {len(df_sig)}")
    print(f"  Saved: {out_csv}")
    if len(df_sig):
        top = df_sig.reindex(df_sig["coef"].abs().sort_values(ascending=False).index).head(15)
        print(top.to_string())


run_lda_for("high_mvpa_bmi_uncouple", "high_mvpa_bmi")
run_lda_for("low_mvpa_bmi_uncouple",  "low_mvpa_bmi")

print("\nDone.")






""
AABC rfMRI LDA — replication of UKB rfMRI ICA LDA, on Glasser-360 FC matrices
collapsed to a network-level matrix.

Two network choices supported, selectable via CLI:

    --net cole12 (default)
        Cole-Anticevic Brain-wide Network Partition (Ji et al. 2019 NeuroImage).
        Each Glasser parcel is uniquely assigned to one of 12 networks.
        Direct lookup, 78-D upper triangle.

    --net yeo17
        Yeo et al. 2011 17-network parcellation, mapped to Glasser parcels by
        majority-vote of the Schaefer-400 17-network labels in conte69 vertex
        space (uses ENIGMA toolbox Schaefer/Glasser CSVs as the geometry
        reference). 153-D upper triangle.

Both run the same regularised LDA pipeline (lsqr + Ledoit-Wolf shrinkage,
n_BS=500, in-fold deconfound for age + sex + site dummies + race_bin +
edu_num, 5/95 percentile bootstrap CIs).
"""
import os
import re
import io
import ssl
import time
import argparse
import urllib.request
import warnings
import numpy as np
import pandas as pd
import scipy.io as sio
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.utils import resample
from nilearn.signal import clean

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---- paths ----
RFMRI_DIR = (r"C:/Users/张杰/Desktop/Cohort Data/AABC/"
             r"AABC_Release2_rfMRI_REST_FullCovarianceConnectivity")
AABC_NIMG = (r"C:/Users/张杰/Desktop/Cohort Data/AABC/"
             r"AABC_Release2_Non-imaging_Data-XL.csv")
PHENO_CSV = r"C:/Users/张杰/Desktop/PA_obesity_biological/data/aabc_pa_obesity_data.csv"
OUT_DIR   = r"C:/Users/张杰/Desktop/PA_obesity_biological/results/AABC_replication"
CACHE_DIR = r"C:/Users/张杰/Desktop/PA_obesity_biological/data"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

# ---- canonical Glasser order (LH 1..180 then RH 181..360) ----
# Same convention as Cole-Anticevic netassignments.  AABC's column order is
# RH (R_V1..R_p24) then LH (L_V1..L_p24).  Build a re-ordering vector.
def aabc_to_glasser_order(aabc_cols):
    """Return integer indices that, applied to a length-380 AABC vector
    (after dropping the x___ subj-id column), yield Glasser canonical order
    LH 1..180 then RH 181..360, then the 19 subcortical IDs preserved at the end.
    Returns: idx (length 360 cortical), subcort_names (list of 19 strs)."""
    cort_names = [c for c in aabc_cols if c.startswith(("L_", "R_"))]
    subc_names = [c for c in aabc_cols if not c.startswith(("L_", "R_")) and c != "x___"]
    # Build LH then RH cortical order:
    lh = [c for c in cort_names if c.startswith("L_")]
    rh = [c for c in cort_names if c.startswith("R_")]
    glasser_order = lh + rh                                    # length 360
    idx = [aabc_cols.index(n) for n in glasser_order]          # AABC col idx
    return idx, glasser_order, subc_names


# ---- Cole-Anticevic 12-network mapping ----
def fetch_netassignments():
    cache = os.path.join(CACHE_DIR, "ColeAnticevic_netassignments.npy")
    if os.path.exists(cache):
        return np.load(cache)
    url = ("https://raw.githubusercontent.com/ColeLab/ColeAnticevicNetPartition/"
           "master/cortex_parcel_network_assignments.mat")
    ctx = ssl.create_default_context()
    blob = urllib.request.urlopen(url, context=ctx, timeout=20).read()
    mat = sio.loadmat(io.BytesIO(blob))
    nets = np.asarray(mat["netassignments"]).squeeze().astype(int)
    np.save(cache, nets)
    return nets


NET_NAMES_COLE = {
    1: "Visual1", 2: "Visual2", 3: "Somatomotor",
    4: "Cingulo-Opercular", 5: "Dorsal-attention",
    6: "Language", 7: "Frontoparietal", 8: "Auditory",
    9: "Default", 10: "Posterior-Multimodal",
    11: "Ventral-Multimodal", 12: "Orbito-Affective",
}

# ---- Yeo 17 networks (canonical ordering, Schaefer 2018 Cereb Cortex Table S4)
NET_NAMES_YEO17 = {
    1:  "VisCent",   2:  "VisPeri",   3:  "SomMotA",   4:  "SomMotB",
    5:  "DorsAttnA", 6:  "DorsAttnB", 7:  "SalVentAttnA", 8: "SalVentAttnB",
    9:  "LimbicA",   10: "LimbicB",   11: "ContA",     12: "ContB",
    13: "ContC",     14: "DefaultA",  15: "DefaultB",  16: "DefaultC",
    17: "TempPar",
}


def fetch_yeo17_for_glasser():
    """Build Glasser-360 -> Yeo 17 mapping by majority vote of the Schaefer-400
    Yeo-17 labels in conte69 vertex space.

    Geometry reference: ENIGMA toolbox bundles
      - glasser_360_conte69.csv   (length 64984, label per vertex, 1..360)
      - schaefer_400_conte69.csv  (length 64984, label per vertex, 1..400)
    Schaefer-400 parcel IDs are ordered to match the 17Networks_order_info.txt
    on the CBIG GitHub repo, which we fetch once and parse to get parcel->net.
    """
    cache = os.path.join(CACHE_DIR, "Glasser360_to_Yeo17.npy")
    if os.path.exists(cache):
        return np.load(cache)

    parc_dir = (r"C:/Users/张杰/AppData/Local/Packages/PythonSoftwareFoundation."
                r"Python.3.11_qbz5n2kfra8p0/LocalCache/local-packages/Python311/"
                r"site-packages/enigmatoolbox/datasets/parcellations")
    glasser = np.genfromtxt(os.path.join(parc_dir, "glasser_360_conte69.csv"),
                            delimiter=",").astype(int)
    schaefer = np.genfromtxt(os.path.join(parc_dir, "schaefer_400_conte69.csv"),
                             delimiter=",").astype(int)

    # Schaefer-400 parcel-id -> Yeo 17 net (1..17), via official order file.
    url = ("https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/"
           "stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/"
           "Parcellations/HCP/fslr32k/cifti/"
           "Schaefer2018_400Parcels_17Networks_order_info.txt")
    txt = urllib.request.urlopen(url, context=ssl.create_default_context(),
                                 timeout=15).read().decode("utf-8", errors="replace")
    # Records appear in pairs: line1 = name, line2 = "<id> r g b a"
    yeo_label_to_int = {name: i + 1 for i, name in enumerate(NET_NAMES_YEO17.values())}
    schaefer_to_yeo = {}
    lines = [l.strip() for l in txt.splitlines() if l.strip()]
    for i in range(0, len(lines), 2):
        name = lines[i]
        # Format: 17Networks_LH_VisCent_ExStr_3
        m = re.match(r"17Networks_(LH|RH)_([A-Za-z]+)_", name)
        if not m:
            continue
        yeo_short = m.group(2)  # e.g. VisCent
        yeo_id = yeo_label_to_int.get(yeo_short)
        if yeo_id is None:
            continue
        parc_id = int(lines[i + 1].split()[0])
        schaefer_to_yeo[parc_id] = yeo_id

    # Vertex-wise: where Schaefer is 0, drop. Where Glasser is 0, drop.
    nets_per_vertex = np.zeros_like(schaefer, dtype=int)
    for sid, yid in schaefer_to_yeo.items():
        nets_per_vertex[schaefer == sid] = yid

    nets = np.zeros(360, dtype=int)
    for pid in range(1, 361):
        m = (glasser == pid)
        labels_in_parcel = nets_per_vertex[m]
        labels_in_parcel = labels_in_parcel[labels_in_parcel > 0]
        if len(labels_in_parcel) == 0:
            nets[pid - 1] = 0  # unmapped — will be excluded
            continue
        vals, counts = np.unique(labels_in_parcel, return_counts=True)
        nets[pid - 1] = vals[np.argmax(counts)]

    n_unmapped = int((nets == 0).sum())
    if n_unmapped:
        print(f"  warning: {n_unmapped} Glasser parcels have no Yeo-17 vertex "
              f"overlap (will be excluded from the network matrix).")
    np.save(cache, nets)
    return nets


def get_network_assignment(net_choice):
    if net_choice == "cole12":
        nets = fetch_netassignments()
        names = NET_NAMES_COLE
    elif net_choice == "yeo17":
        nets = fetch_yeo17_for_glasser()
        names = NET_NAMES_YEO17
    else:
        raise ValueError(f"unknown net_choice: {net_choice}")
    return nets, names


def build_network_matrix(fc360, nets, n_nets):
    """fc360: (360,360) symmetric matrix in Glasser canonical order.
    nets:  (360,) integer network labels 1..n_nets (0 = unassigned).
    Returns: (n_nets, n_nets) matrix.  Within-network: mean of upper-triangle
    pairs (i!=j) inside net.  Between-network: mean of all i-in-A x j-in-B."""
    out = np.zeros((n_nets, n_nets), dtype=np.float64)
    masks = [nets == k for k in range(1, n_nets + 1)]
    for a in range(n_nets):
        ma = masks[a]
        if not ma.any():
            continue
        for b in range(a, n_nets):
            mb = masks[b]
            if not mb.any():
                continue
            sub = fc360[np.ix_(ma, mb)]
            if a == b:
                idx = np.triu_indices_from(sub, k=1)
                vals = sub[idx]
            else:
                vals = sub.ravel()
            m = float(np.nanmean(vals)) if vals.size else 0.0
            out[a, b] = m
            out[b, a] = m
    return out


def upper_tri_with_diag(M, n):
    """Return upper triangle (incl. diag) of (n,n) matrix as (n*(n+1)/2,)."""
    iu = np.triu_indices(n, k=0)
    return M[iu]


# ---- Filename helper: column names like 'L_p9_46v_ROI' map to file names with
#      'L_p9-46v_ROI' (dash instead of underscore).  Build the mapping once.
def col_to_filename_seed(col_name, available_filenames):
    """Try a few simple substitutions to find the seed file matching a column."""
    if col_name in available_filenames:
        return col_name
    # Replace one '_' with '-' in patterns like 'p9_46v', 'a9_46v', '9_46d',
    # 'i6_8', 's6_8', 'OP2_3' (Glasser irregular IDs).
    for n_swap in (1,):
        # Try replacing each '_' (other than leading L_/R_ or trailing _ROI) with '-'
        # Simpler: regex sub for digit_digit where the pair is part of an irregular Glasser name
        candidates = []
        candidates.append(re.sub(r"(\d)_(\d)", r"\1-\2", col_name))
        candidates.append(re.sub(r"OP2_3", "OP2-3", col_name))
        candidates.append(re.sub(r"i6_8",  "i6-8",  col_name))
        candidates.append(re.sub(r"s6_8",  "s6-8",  col_name))
        candidates.append(re.sub(r"9_46",  "9-46",  col_name))
        for c in candidates:
            if c in available_filenames:
                return c
    raise KeyError(f"Cannot find seed file for col {col_name}")


# ---- Stage 2: deconfound + LDA ----
def perform_bootstrap(Y, n_samples):
    df_y_test  = Y.sample(frac=0.1, random_state=0)
    df_y_train = Y[~Y.index.isin(df_y_test.index)]
    minority   = df_y_train.value_counts().min()
    samples = []
    for s in range(n_samples):
        np.random.seed(s)
        c0 = resample(df_y_train[df_y_train == 0].dropna(),
                      n_samples=minority, replace=True)
        c1 = resample(df_y_train[df_y_train == 1].dropna(),
                      n_samples=minority, replace=True)
        samples.append(pd.concat([c0, c1]).index.to_list())
    return samples, df_y_test.index


def calculate_feature_coef(X, Y, n_BS, covariates=None):
    samples, _ = perform_bootstrap(Y, n_BS)
    lda = LDA(solver="lsqr", shrinkage="auto")
    out = np.zeros((n_BS, X.shape[1]))
    for i, idxs in enumerate(samples):
        X_train = X.loc[idxs].values
        y_train = np.ravel(Y.loc[idxs].values)
        if covariates is not None:
            cov_train = covariates.loc[idxs].values
            X_train = clean(X_train, confounds=cov_train,
                            detrend=False, standardize=False)
        out[i, :] = lda.fit(X_train, y_train).coef_
    return out


def calculate_significant_coef(coef_mat, feat_names, alpha):
    coef_mean = coef_mat.mean(axis=0)
    lo = np.percentile(coef_mat, alpha / 2,        axis=0)
    hi = np.percentile(coef_mat, 100 - alpha / 2,  axis=0)
    df_all = pd.DataFrame({"coef": coef_mean, "lower": lo, "upper": hi},
                          index=feat_names)
    sig = df_all[(lo > 0) | (hi < 0)].copy()
    return sig, df_all


def run_lda_for(target_col, tag, X_df, pheno, covars, net_choice,
                n_BS=500, alpha=10):
    Y = pheno[target_col].dropna()
    common_idx = Y.index.intersection(X_df.index)
    Y = Y.loc[common_idx].astype(int)
    X = X_df.loc[common_idx]
    cov = pheno.loc[common_idx, covars].astype(float).fillna(
        pheno[covars].astype(float).median())
    print(f"\n=== {tag} [{net_choice}] ===")
    print(f"  N (overlap with rfMRI): {len(common_idx)}")
    print(f"  class counts:\n{Y.value_counts()}")
    coef = calculate_feature_coef(X, Y, n_BS, covariates=cov)
    sig, df_all = calculate_significant_coef(coef, X.columns, alpha)
    out_csv = os.path.join(OUT_DIR,
        f"LDA_rfMRI_{net_choice}_{tag}_coef_5-95CI_all.csv")
    df_all.to_csv(out_csv, index_label="feature")
    print(f"  Significant features (CI excludes 0): {len(sig)}")
    print(f"  Saved: {out_csv}")
    if len(sig):
        top = sig.reindex(sig["coef"].abs().sort_values(ascending=False).index).head(15)
        print(top.to_string())


# ---- Stage 3: deconfound the n_feat-D feature matrix at the brain-feature stage
def deconfound_brain(X_df, ni_v1):
    """X_df: (N,n_features) DataFrame with subject IDs as index.
    Standard-scale, then regress out site dummies."""
    sites = pd.get_dummies(ni_v1.loc[X_df.index, "site"],
                           drop_first=False).astype(float).values
    sc = StandardScaler()
    Xz = sc.fit_transform(X_df.values)
    X_dec = clean(Xz, confounds=sites, detrend=False, standardize=False)
    return pd.DataFrame(X_dec, index=X_df.index, columns=X_df.columns)


# ---- Stage 1: per-subject 360x360 FC -> nxn network matrix ----
def assemble_network_features(aabc_dir, net_choice):
    nets, names = get_network_assignment(net_choice)
    n_nets = len(names)
    n_feat = n_nets * (n_nets + 1) // 2

    out_cache = os.path.join(CACHE_DIR, f"aabc_rfmri_v1_{net_choice}_fc{n_feat}.npz")
    if os.path.exists(out_cache):
        data = np.load(out_cache, allow_pickle=True)
        return data["X"], data["ids"], data["labels"], names

    sample_seed = "L_V1_ROI"
    sample_path = os.path.join(aabc_dir,
        f"rfMRI_REST_FullCovarianceConnectivity_{sample_seed}.csv")
    sample = pd.read_csv(sample_path)
    aabc_cols = sample.columns.tolist()
    cort_idx, glasser_order, subc_names = aabc_to_glasser_order(aabc_cols)

    sample["subj"] = sample["x___"].str.split("_").str[0]
    sample["event"] = sample["x___"].str.split("_").str[1]
    v1_mask = sample["event"] == "V1"
    ids = sample.loc[v1_mask, "subj"].values
    n_subj = len(ids)
    print(f"V1 baseline rows: {n_subj}")
    print(f"Network choice: {net_choice}, {n_nets} networks, {n_feat} features.")

    chunk_size = 200
    n_chunks = int(np.ceil(n_subj / chunk_size))

    pair_labels = []
    for a in range(1, n_nets + 1):
        for b in range(a, n_nets + 1):
            pair_labels.append(f"{names[a]}__{names[b]}")
    pair_labels = np.array(pair_labels)

    X = np.zeros((n_subj, n_feat), dtype=np.float64)

    print(f"Streaming {n_subj} subjects in {n_chunks} chunk(s)...")
    avail_seeds = {f.replace("rfMRI_REST_FullCovarianceConnectivity_", "").replace(".csv", "")
                   for f in os.listdir(aabc_dir)
                   if f.startswith("rfMRI_REST_FullCovarianceConnectivity_")}
    col_to_file = {c: col_to_filename_seed(c, avail_seeds) for c in glasser_order}
    for ci in range(n_chunks):
        s0 = ci * chunk_size
        s1 = min((ci + 1) * chunk_size, n_subj)
        chunk_ids = ids[s0:s1]
        chunk_n = s1 - s0
        fc_chunk = np.zeros((chunk_n, 360, 360), dtype=np.float32)
        t0 = time.time()
        for seed_pos, seed_name in enumerate(glasser_order):
            seed_csv = os.path.join(aabc_dir,
                f"rfMRI_REST_FullCovarianceConnectivity_{col_to_file[seed_name]}.csv")
            df = pd.read_csv(seed_csv)
            df["subj"]  = df["x___"].str.split("_").str[0]
            df["event"] = df["x___"].str.split("_").str[1]
            df = df[df["event"] == "V1"]
            df = df.set_index("subj")
            sub = df.loc[df.index.intersection(chunk_ids), glasser_order]
            row_map = {s: i for i, s in enumerate(chunk_ids)}
            for s in sub.index:
                fc_chunk[row_map[s], seed_pos, :] = sub.loc[s].values.astype(np.float32)
        fc_chunk = 0.5 * (fc_chunk + fc_chunk.transpose(0, 2, 1))
        for k in range(chunk_n):
            M = build_network_matrix(fc_chunk[k], nets, n_nets)
            X[s0 + k] = upper_tri_with_diag(M, n_nets)
        print(f"  chunk {ci+1}/{n_chunks}  ({s1-s0} subj)  elapsed {time.time()-t0:.0f}s")

    np.savez_compressed(out_cache, X=X, ids=ids, labels=pair_labels)
    print(f"Cache saved: {out_cache}")
    return X, ids, pair_labels, names


# ============================================================
# Main
# ============================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", choices=["cole12", "yeo17"], default="cole12",
                    help="Network parcellation to collapse Glasser-360 onto.")
    args = ap.parse_args()
    net_choice = args.net

    print(f">> Stage 1: build per-subject network matrix [{net_choice}]")
    X_arr, ids, labels, _names = assemble_network_features(RFMRI_DIR, net_choice)
    print(f"X shape: {X_arr.shape}, N subjects: {len(ids)}, features: {len(labels)}")

    X_df = pd.DataFrame(X_arr, index=pd.Index(ids, name="sub_id"), columns=labels)
    n_before = len(X_df)
    X_df = X_df.dropna(how="any")
    print(f"After dropping NaN rows: {len(X_df)} (was {n_before})")

    print("\n>> Stage 2: site-deconfound brain features")
    ni = pd.read_csv(AABC_NIMG, low_memory=False)
    ni.columns = [c.split(" - ")[0] for c in ni.columns]
    ni = ni[ni["id"] != "id"]
    ni_v1 = ni[ni["event"] == "V1"][["id", "site"]].drop_duplicates(subset="id").set_index("id")
    common = X_df.index.intersection(ni_v1.index)
    X_df = X_df.loc[common]
    print(f"After site lookup: {len(X_df)}")
    X_df_dec = deconfound_brain(X_df, ni_v1)
    out_dec = os.path.join(OUT_DIR, f"aabc_rfMRI_{net_choice}_deconfounded.csv")
    X_df_dec.to_csv(out_dec)
    print(f"Saved deconfounded matrix: {out_dec}  shape={X_df_dec.shape}")

    print("\n>> Stage 3: load phenotypes + covariates")
    pheno = pd.read_csv(PHENO_CSV).rename(columns={"id": "sub_id"}).set_index("sub_id")
    pheno["edu_num"] = pd.Categorical(pheno["education"]).codes
    race_dum = pd.get_dummies(pheno["race_bin"], drop_first=True).astype(float)
    race_dum.columns = [f"raceD_{c}" for c in race_dum.columns]
    pheno = pd.concat([pheno, race_dum], axis=1)
    COVS = ["age_open", "sex_num", "edu_num"] + race_dum.columns.tolist()

    print("\n>> Stage 4: run LDA for HD and LD")
    run_lda_for("high_mvpa_bmi_uncouple", "high_mvpa_bmi", X_df_dec, pheno, COVS, net_choice)
    run_lda_for("low_mvpa_bmi_uncouple",  "low_mvpa_bmi",  X_df_dec, pheno, COVS, net_choice)
    print("\nDone.")


if __name__ == "__main__":
    main()
