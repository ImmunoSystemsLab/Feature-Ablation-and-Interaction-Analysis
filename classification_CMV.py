import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import roc_auc_score
from scipy.sparse import hstack

# قراءة بيانات Human baseline + CMV
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_cmv   = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv", sep="\t")

df_human["Label"] = 0
df_cmv["Label"]   = 1
df = pd.concat([df_human, df_cmv]).reset_index(drop=True)

# عينة للتجربة
df_small = df.sample(n=50000, random_state=42).reset_index(drop=True)

# Features رقمية
df_small["len_CDR3a"] = df_small["CDR3a"].astype(str).apply(len)
df_small["len_CDR3b"] = df_small["CDR3b"].astype(str).apply(len)

# استخراج motif (مثال: آخر 4 residues من CDR3β)
df_small["motif"] = df_small["CDR3b"].astype(str).str[-4:]

# تعريف مجموعات الـ features
feature_sets = {
    "TRBV only": ["TRBV"],
    "Length only": ["len_CDR3b"],
    "Motif only": ["motif"],
    "TRBV + Length": ["TRBV","len_CDR3b"],
    "TRBV + Motif": ["TRBV","motif"],
    "TRBV + TRAV": ["TRBV","TRAV"],
    "TRBV + Motif + Length": ["TRBV","motif","len_CDR3b"],
    "Full model": ["TRAV","TRBV","motif","len_CDR3a","len_CDR3b"]
}

def crossval_eval(features, seeds=range(3)):
    aucs = []
    for seed in seeds:
        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
        for train_idx, test_idx in skf.split(df_small[features], df_small["Label"]):
            X_train, X_test = df_small.iloc[train_idx][features], df_small.iloc[test_idx][features]
            y_train, y_test = df_small.iloc[train_idx]["Label"], df_small.iloc[test_idx]["Label"]

            cat_cols = [c for c in X_train.columns if X_train[c].dtype == "object"]
            num_cols = [c for c in X_train.columns if c not in cat_cols]

            enc = OneHotEncoder(handle_unknown="ignore", min_frequency=5)
            X_train_cat = enc.fit_transform(X_train[cat_cols].astype(str))
            X_test_cat  = enc.transform(X_test[cat_cols].astype(str))

            X_train_final = hstack([X_train_cat, X_train[num_cols].values])
            X_test_final  = hstack([X_test_cat, X_test[num_cols].values])

            clf = SGDClassifier(loss="log_loss", max_iter=1000)
            clf.fit(X_train_final, y_train)
            y_prob = clf.predict_proba(X_test_final)[:,1]

            aucs.append(roc_auc_score(y_test, y_prob))

    return np.mean(aucs)

# تشغيل كل الموديلات
results = {}
for name, feats in feature_sets.items():
    results[name] = crossval_eval(feats)

# عرض النتائج في جدول
df_results = pd.DataFrame(results, index=["AUROC"]).T
print("\n=== نتائج CMV (Feature Ablation/Interaction) ===")
print(df_results)
