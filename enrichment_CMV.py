import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_cmv = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل gene
human_freq = df_human["TRBV"].value_counts(normalize=True)
cmv_freq = df_cmv["TRBV"].value_counts(normalize=True)

print("human Freq TRBV Gene Enrichment",human_freq)
print("CMV Freq TRBV Gene Enrichment",cmv_freq)


# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "CMV_freq": cmv_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["CMV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالـ genes اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

print("TRBV Gene Enrichment in CMV relative to Human baseline",df)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (CMV vs Human)")
plt.xlabel("TRBV Gene")
plt.title("TRBV Gene Enrichment in CMV relative to Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()






# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_cmv = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل طول
human_freq = df_human["CDR3b"].str.len().value_counts(normalize=True)
cmv_freq = df_cmv["CDR3b"].str.len().value_counts(normalize=True)

print("human Freq CDR3b Gene Enrichment",human_freq)
print("CMV Freq CDR3b Gene Enrichment",cmv_freq)

# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "CMV_freq": cmv_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["CMV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالأطوال اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

print("CDR3 Length Enrichment in CMV relative to Human baseline",df)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (CMV vs Human)")
plt.xlabel("CDR3 Length (aa)")
plt.title("CDR3 Length Enrichment in CMV relative to Human baseline")
plt.tight_layout()
plt.show()















import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_cmv = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل gene
human_freq = df_human["TRAV"].value_counts(normalize=True)
cmv_freq = df_cmv["TRAV"].value_counts(normalize=True)

print("human Freq TRAV Gene Enrichment",human_freq)
print("CMV Freq TRAV Gene Enrichment",cmv_freq)


# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "CMV_freq": cmv_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["CMV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالـ genes اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

print("TRBV Gene Enrichment in CMV relative to Human baseline",df)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (CMV vs Human)")
plt.xlabel("TRBV Gene")
plt.title("TRBV Gene Enrichment in CMV relative to Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()








# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_cmv = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل طول
human_freq = df_human["CDR3a"].str.len().value_counts(normalize=True)
cmv_freq = df_cmv["CDR3a"].str.len().value_counts(normalize=True)

print("human Freq CDR3a Gene Enrichment",human_freq)
print("CMV Freq CDR3a Gene Enrichment",cmv_freq)

# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "CMV_freq": cmv_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["CMV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالأطوال اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

print("CDR3 Length Enrichment in CMV relative to Human baseline",df)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (CMV vs Human)")
plt.xlabel("CDR3 Length (aa)")
plt.title("CDR3 Length Enrichment in CMV relative to Human baseline")
plt.tight_layout()
plt.show()








import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_cmv   = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv", sep="\t")

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
cdr3b_cmv   = df_cmv["CDR3b"].dropna().tolist()

motifs_b_human = extract_motifs(cdr3b_human, k=4)
motifs_b_cmv   = extract_motifs(cdr3b_cmv, k=4)

human_counts_b = Counter(motifs_b_human)
cmv_counts_b   = Counter(motifs_b_cmv)

motif_df_b = pd.DataFrame({
    "Human_freq": {m: human_counts_b[m]/len(motifs_b_human) for m in human_counts_b},
    "CMV_freq":   {m: cmv_counts_b[m]/len(motifs_b_cmv) for m in cmv_counts_b}
}).fillna(0)

motif_df_b["log2_enrichment"] = np.log2((motif_df_b["CMV_freq"]+1e-9) / (motif_df_b["Human_freq"]+1e-9))
motif_df_b = motif_df_b[motif_df_b["Human_freq"] > 0.001].sort_values("log2_enrichment", ascending=False).head(20)

print("\nTop β-chain motifs enriched in CMV vs Human:\n", motif_df_b)

plt.figure(figsize=(12,6))
sns.barplot(x=motif_df_b.index, y=motif_df_b["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (CMV vs Human)")
plt.xlabel("Motif (CDR3β)")
plt.title("Top β-chain motifs enriched in CMV vs Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# ============================
# 2) Motifs من CDR3α (alpha chain)
# ============================
cdr3a_human = df_human["CDR3a"].dropna().tolist()
cdr3a_cmv   = df_cmv["CDR3a"].dropna().tolist()

motifs_a_human = extract_motifs(cdr3a_human, k=4)
motifs_a_cmv   = extract_motifs(cdr3a_cmv, k=4)

human_counts_a = Counter(motifs_a_human)
cmv_counts_a   = Counter(motifs_a_cmv)

motif_df_a = pd.DataFrame({
    "Human_freq": {m: human_counts_a[m]/len(motifs_a_human) for m in human_counts_a},
    "CMV_freq":   {m: cmv_counts_a[m]/len(motifs_a_cmv) for m in cmv_counts_a}
}).fillna(0)

motif_df_a["log2_enrichment"] = np.log2((motif_df_a["CMV_freq"]+1e-9) / (motif_df_a["Human_freq"]+1e-9))
motif_df_a = motif_df_a[motif_df_a["Human_freq"] > 0.001].sort_values("log2_enrichment", ascending=False).head(20)

print("\nTop α-chain motifs enriched in CMV vs Human:\n", motif_df_a)

plt.figure(figsize=(12,6))
sns.barplot(x=motif_df_a.index, y=motif_df_a["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (CMV vs Human)")
plt.xlabel("Motif (CDR3α)")
plt.title("Top α-chain motifs enriched in CMV vs Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
