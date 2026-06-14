import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer

print(">>> Starting training script...")

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_ebv   = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")
df_cmv   = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv", sep="\t")
print(">>> Files loaded successfully")

# إضافة label
df_human["label"] = 0
df_ebv["label"]   = 1
df_cmv["label"]   = 2

# دمج
df = pd.concat([df_human, df_ebv, df_cmv], ignore_index=True)
print(">>> Data merged, shape:", df.shape)

# ✅ أخذ عينة صغيرة للتجربة (50 ألف صف)
df = df.sample(n=50000, random_state=42)
print(">>> Sampled subset, shape:", df.shape)

# استخراج موتيف
def extract_motif(seq, k=4):
    if isinstance(seq, str) and len(seq) >= k:
        return seq[-k:]
    return "Unknown"

df["Motif"] = df["CDR3b"].apply(lambda x: extract_motif(x, k=4))
df["CDR3a_len"] = df["CDR3a"].str.len()
df["CDR3b_len"] = df["CDR3b"].str.len()
print(">>> Features engineered")

# اختيار الأعمدة
X = df[["TRAV","TRBV","CDR3a_len","CDR3b_len","Motif"]].copy()
y = df["label"]

# معالجة NaN
X = X.fillna({
    "TRAV": "Unknown",
    "TRBV": "Unknown",
    "Motif": "Unknown",
    "CDR3a_len": 0,
    "CDR3b_len": 0
})
print(">>> NaN handled, shape:", X.shape)

# One-Hot Encoding
categorical = ["TRAV","TRBV","Motif"]
numeric = ["CDR3a_len","CDR3b_len"]

preprocessor = ColumnTransformer(
    transformers=[
        ("cat", OneHotEncoder(handle_unknown="ignore"), categorical),
        ("num", SimpleImputer(strategy="median"), numeric)
    ]
)

# Random Forest pipeline (20 شجرة فقط)
model = Pipeline(steps=[
    ("preprocessor", preprocessor),
    ("classifier", RandomForestClassifier(
        n_estimators=20,
        random_state=42,
        class_weight="balanced"
    ))
])

# Train/Test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)
print(">>> Train/Test split done:", X_train.shape, X_test.shape)

# تدريب النموذج
print(">>> Training Random Forest...")
model.fit(X_train, y_train)
print(">>> Training complete")

# تقييم النموذج
y_pred_prob = model.predict_proba(X_test)
print("\nROC-AUC (macro):", roc_auc_score(y_test, y_pred_prob, multi_class="ovr"))
print("\nClassification Report:\n", classification_report(y_test, model.predict(X_test)))

# مثال prediction
sample = pd.DataFrame([{
    "TRAV": "TRAV13-1",
    "TRBV": "TRBV7-9",
    "CDR3a_len": 12,
    "CDR3b_len": 14,
    "Motif": "TGEL"
}])

prob = model.predict_proba(sample)[0]
print("\nPrediction Probabilities (Human, EBV, CMV):", prob)
print("P(EBV specificity) =", prob[1])
print("P(CMV specificity) =", prob[2])

# استخراج feature importance
feature_names = model.named_steps["preprocessor"].get_feature_names_out()
importances = model.named_steps["classifier"].feature_importances_

feat_imp = pd.DataFrame({"feature": feature_names, "importance": importances})
feat_imp = feat_imp.sort_values("importance", ascending=False).head(20)

print("\nTop 20 important features:\n", feat_imp)

# Visualization
plt.figure(figsize=(10,6))
plt.barh(feat_imp["feature"], feat_imp["importance"], color="skyblue")
plt.xlabel("Importance")
plt.title("Top 20 Feature Importances (Random Forest)")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()





import shap
import numpy as np

rf_model = model.named_steps["classifier"]
preprocessor = model.named_steps["preprocessor"]

# تحويل البيانات إلى dense float array
X_transformed = preprocessor.transform(X_test)
X_transformed = X_transformed.toarray().astype(np.float64)

# TreeExplainer مع تعطيل additivity check
explainer = shap.TreeExplainer(rf_model, feature_perturbation="interventional")

# أول 3 عينات
samples = X_transformed[:3]
shap_values = explainer.shap_values(samples, check_additivity=False)

feature_names = preprocessor.get_feature_names_out()

for i in range(len(samples)):
    print(f"\n>>> SHAP explanation for test sample {i+1}:")
    for cls, values in enumerate(shap_values):
        print(f"\nClass {cls} contributions:")
        for f, val in zip(feature_names, values[i]):
            if abs(val) > 0.05:
                print(f"{f} contributed {val:+.2f}")

