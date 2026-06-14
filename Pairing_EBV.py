import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_ebv = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

### 4) Interaction: TRAV–TRBV pairing
df_human["TRAV_TRBV"] = df_human["TRAV"].astype(str) + "_" + df_human["TRBV"].astype(str)
df_ebv["TRAV_TRBV"] = df_ebv["TRAV"].astype(str) + "_" + df_ebv["TRBV"].astype(str)

human_freq = df_human["TRAV_TRBV"].value_counts(normalize=True)
ebv_freq = df_ebv["TRAV_TRBV"].value_counts(normalize=True)

df_pair = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": ebv_freq}).fillna(0)
df_pair["fold_change"] = (df_pair["EBV_freq"]+1e-9) / (df_pair["Human_freq"]+1e-9)
df_pair["log2_enrichment"] = np.log2(df_pair["fold_change"])
df_pair = df_pair[df_pair["Human_freq"] > 0.01].sort_values("log2_enrichment", ascending=False)

print("=== TRAV–TRBV Pairing Enrichment (EBV vs Human) ===")
print(df_pair[["Human_freq","EBV_freq","fold_change","log2_enrichment"]].head(20))

plt.figure(figsize=(12,8))
sns.heatmap(df_pair.pivot_table(values="log2_enrichment",
                                index=df_pair.index.str.split("_").str[0],
                                columns=df_pair.index.str.split("_").str[1],
                                fill_value=0),
            cmap="coolwarm", center=0)
plt.title("TRAV–TRBV Pairing Enrichment (EBV vs Human)")
plt.tight_layout()
plt.show()









### 5) Interaction: TRBV + CDR3β length
df_human["TRBV_CDR3bLen"] = df_human["TRBV"].astype(str) + "_len" + df_human["CDR3b"].str.len().astype(str)
df_ebv["TRBV_CDR3bLen"] = df_ebv["TRBV"].astype(str) + "_len" + df_ebv["CDR3b"].str.len().astype(str)

human_freq = df_human["TRBV_CDR3bLen"].value_counts(normalize=True)
ebv_freq = df_ebv["TRBV_CDR3bLen"].value_counts(normalize=True)

df_pair2 = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": ebv_freq}).fillna(0)
df_pair2["fold_change"] = (df_pair2["EBV_freq"]+1e-9) / (df_pair2["Human_freq"]+1e-9)
df_pair2["log2_enrichment"] = np.log2(df_pair2["fold_change"])
df_pair2 = df_pair2[df_pair2["Human_freq"] > 0.01].sort_values("log2_enrichment", ascending=False)

print("\n=== TRBV + CDR3β Length Interaction Enrichment (EBV vs Human) ===")
print(df_pair2[["Human_freq","EBV_freq","fold_change","log2_enrichment"]].head(20))

plt.figure(figsize=(12,6))
sns.barplot(x=df_pair2.index, y=df_pair2["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("TRBV + CDR3β length")
plt.title("TRBV–CDR3β Length Interaction Enrichment (EBV vs Human)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()











### V–J pairing (TRBV–TRBJ)
df_human["TRBV_TRBJ"] = df_human["TRBV"].astype(str) + "_" + df_human["TRBJ"].astype(str)
df_ebv["TRBV_TRBJ"] = df_ebv["TRBV"].astype(str) + "_" + df_ebv["TRBJ"].astype(str)

human_freq = df_human["TRBV_TRBJ"].value_counts(normalize=True)
ebv_freq = df_ebv["TRBV_TRBJ"].value_counts(normalize=True)

df_vj = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": ebv_freq}).fillna(0)
df_vj["fold_change"] = (df_vj["EBV_freq"]+1e-9) / (df_vj["Human_freq"]+1e-9)
df_vj["log2_enrichment"] = np.log2(df_vj["fold_change"])
df_vj = df_vj[df_vj["Human_freq"] > 0.01].sort_values("log2_enrichment", ascending=False)

print("=== TRBV–TRBJ Pairing Enrichment (EBV vs Human) ===")
print(df_vj[["Human_freq","EBV_freq","fold_change","log2_enrichment"]].head(20))

plt.figure(figsize=(12,8))
sns.heatmap(df_vj.pivot_table(values="log2_enrichment",
                              index=df_vj.index.str.split("_").str[0],
                              columns=df_vj.index.str.split("_").str[1],
                              fill_value=0),
            cmap="coolwarm", center=0)
plt.title("TRBV–TRBJ Pairing Enrichment (EBV vs Human)")
plt.tight_layout()
plt.show()








from collections import Counter

# دالة استخراج موتيف (مثلاً آخر 4 حروف من CDR3β)
def extract_motif(seq, k=4):
    if isinstance(seq, str) and len(seq) >= k:
        return seq[-k:]   # ناخد آخر k residues
    return "NA"

# إضافة عمود للموتيف في Human و EBV
df_human["Motif"] = df_human["CDR3b"].apply(lambda x: extract_motif(x, k=4))
df_ebv["Motif"]  = df_ebv["CDR3b"].apply(lambda x: extract_motif(x, k=4))

# بناء الـ triple مع الموتيف
df_human["TRAV_TRBV_CDR3bLen_Motif"] = (
    df_human["TRAV"].astype(str) + "_" +
    df_human["TRBV"].astype(str) + "_len" +
    df_human["CDR3b"].str.len().astype(str) + "_" +
    df_human["Motif"].astype(str)
)

df_ebv["TRAV_TRBV_CDR3bLen_Motif"] = (
    df_ebv["TRAV"].astype(str) + "_" +
    df_ebv["TRBV"].astype(str) + "_len" +
    df_ebv["CDR3b"].str.len().astype(str) + "_" +
    df_ebv["Motif"].astype(str)
)

# حساب frequencies
human_freq = df_human["TRAV_TRBV_CDR3bLen_Motif"].value_counts(normalize=True)
ebv_freq   = df_ebv["TRAV_TRBV_CDR3bLen_Motif"].value_counts(normalize=True)

# دمج
df_triple_motif = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": ebv_freq}).fillna(0)
df_triple_motif["fold_change"] = (df_triple_motif["EBV_freq"]+1e-9) / (df_triple_motif["Human_freq"]+1e-9)
df_triple_motif["log2_enrichment"] = np.log2(df_triple_motif["fold_change"])

# فلترة: نحتفظ بالـ combos اللي Human_freq > 0.001
df_triple_motif = df_triple_motif[df_triple_motif["Human_freq"] > 0.001].sort_values("log2_enrichment", ascending=False)

print("=== TRAV × TRBV × CDR3β length × Motif Enrichment (EBV vs Human) ===")
print(df_triple_motif[["Human_freq","EBV_freq","fold_change","log2_enrichment"]].head(20))

# رسم barplot
plt.figure(figsize=(14,6))
sns.barplot(x=df_triple_motif.index[:20], y=df_triple_motif["log2_enrichment"].head(20), palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("TRAV × TRBV × CDR3β length × Motif")
plt.title("Triple Interaction + Motif Enrichment (EBV vs Human)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
