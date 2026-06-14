import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_EBV = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل gene
human_freq = df_human["TRBV"].value_counts(normalize=True)
EBV_freq = df_EBV["TRBV"].value_counts(normalize=True)

print("Human Freq TRBV Gene Enrichment", human_freq)
print("EBV Freq TRBV Gene Enrichment", EBV_freq)

# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": EBV_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["EBV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالـ genes اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

print("TRBV Gene Enrichment in EBV relative to Human baseline",df)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("TRBV Gene")
plt.title("TRBV Gene Enrichment in EBV relative to Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()








# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_EBV = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل طول
human_freq = df_human["CDR3b"].str.len().value_counts(normalize=True)
EBV_freq = df_EBV["CDR3b"].str.len().value_counts(normalize=True)

print("Human Freq CDR3b Length Enrichment", human_freq)
print("EBV Freq CDR3b Length Enrichment", EBV_freq)

# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": EBV_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["EBV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالأطوال اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

print("CDR3 Length Enrichment in EBV relative to Human baseline",df)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("CDR3 Length (aa)")
plt.title("CDR3 Length Enrichment in EBV relative to Human baseline")
plt.tight_layout()
plt.show()










import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_ebv = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# ============================
# 1) TRAV Gene Enrichment (α chain)
# ============================
human_freq = df_human["TRAV"].value_counts(normalize=True)
ebv_freq = df_ebv["TRAV"].value_counts(normalize=True)

print("Human Freq TRAV Gene Enrichment", human_freq)
print("EBV Freq TRAV Gene Enrichment", ebv_freq)

df_trav = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": ebv_freq}).fillna(0)
df_trav["log2_enrichment"] = np.log2((df_trav["EBV_freq"]+1e-9) / (df_trav["Human_freq"]+1e-9))
df_trav = df_trav[df_trav["Human_freq"] > 0.01].sort_values("log2_enrichment", ascending=False)

print("TRAV Gene Enrichment in EBV relative to Human baseline", df_trav)

plt.figure(figsize=(12,6))
sns.barplot(x=df_trav.index, y=df_trav["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("TRAV Gene")
plt.title("TRAV Gene Enrichment in EBV relative to Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()





# ============================
# 2) CDR3α Length Enrichment
# ============================
human_freq = df_human["CDR3a"].str.len().value_counts(normalize=True)
ebv_freq = df_ebv["CDR3a"].str.len().value_counts(normalize=True)

print("Human Freq CDR3a Length Enrichment", human_freq)
print("EBV Freq CDR3a Length Enrichment", ebv_freq)

df_len = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": ebv_freq}).fillna(0)
df_len["log2_enrichment"] = np.log2((df_len["EBV_freq"]+1e-9) / (df_len["Human_freq"]+1e-9))
df_len = df_len[df_len["Human_freq"] > 0.01].sort_values("log2_enrichment", ascending=False)

print("CDR3α Length Enrichment in EBV relative to Human baseline", df_len)

plt.figure(figsize=(12,6))
sns.barplot(x=df_len.index, y=df_len["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("CDR3α Length (aa)")
plt.title("CDR3α Length Enrichment in EBV relative to Human baseline")
plt.tight_layout()
plt.show()







import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_ebv = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# ============================
# دالة استخراج موتيفات بطول k
# ============================
def extract_motifs(seq_list, k=4):
    motifs = []
    for seq in seq_list:
        if isinstance(seq, str):  # نتأكد إن القيمة نصية
            for i in range(len(seq)-k+1):
                motifs.append(seq[i:i+k])
    return motifs

# ============================
# 1) Motifs من CDR3β (beta chain)
# ============================
cdr3b_human = df_human["CDR3b"].dropna().tolist()
cdr3b_ebv = df_ebv["CDR3b"].dropna().tolist()

motifs_b_human = extract_motifs(cdr3b_human, k=4)
motifs_b_ebv = extract_motifs(cdr3b_ebv, k=4)

human_counts_b = Counter(motifs_b_human)
ebv_counts_b = Counter(motifs_b_ebv)

motif_df_b = pd.DataFrame({
    "Human_freq": {m: human_counts_b[m]/len(motifs_b_human) for m in human_counts_b},
    "EBV_freq": {m: ebv_counts_b[m]/len(motifs_b_ebv) for m in ebv_counts_b}
}).fillna(0)

motif_df_b["log2_enrichment"] = np.log2((motif_df_b["EBV_freq"]+1e-9) / (motif_df_b["Human_freq"]+1e-9))
motif_df_b = motif_df_b[motif_df_b["Human_freq"] > 0.001].sort_values("log2_enrichment", ascending=False).head(20)

print("\nTop β-chain motifs enriched in EBV vs Human:\n", motif_df_b)

plt.figure(figsize=(12,6))
sns.barplot(x=motif_df_b.index, y=motif_df_b["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("Motif (CDR3β)")
plt.title("Top β-chain motifs enriched in EBV vs Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# ============================
# 2) Motifs من CDR3α (alpha chain)
# ============================
cdr3a_human = df_human["CDR3a"].dropna().tolist()
cdr3a_ebv = df_ebv["CDR3a"].dropna().tolist()

motifs_a_human = extract_motifs(cdr3a_human, k=4)
motifs_a_ebv = extract_motifs(cdr3a_ebv, k=4)

human_counts_a = Counter(motifs_a_human)
ebv_counts_a = Counter(motifs_a_ebv)

motif_df_a = pd.DataFrame({
    "Human_freq": {m: human_counts_a[m]/len(motifs_a_human) for m in human_counts_a},
    "EBV_freq": {m: ebv_counts_a[m]/len(motifs_a_ebv) for m in ebv_counts_a}
}).fillna(0)

motif_df_a["log2_enrichment"] = np.log2((motif_df_a["EBV_freq"]+1e-9) / (motif_df_a["Human_freq"]+1e-9))
motif_df_a = motif_df_a[motif_df_a["Human_freq"] > 0.001].sort_values("log2_enrichment", ascending=False).head(20)

print("\nTop α-chain motifs enriched in EBV vs Human:\n", motif_df_a)

plt.figure(figsize=(12,6))
sns.barplot(x=motif_df_a.index, y=motif_df_a["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("Motif (CDR3α)")
plt.title("Top α-chain motifs enriched in EBV vs Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
