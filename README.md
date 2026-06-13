# 🧬 Explainable AI for TCR Antigen Specificity Prediction

## Overview

This project implements an interpretable machine learning pipeline for predicting T cell receptor (TCR) antigen specificity (EBV, CMV, and Human baseline) using explainable AI approaches.

The goal is not only to classify TCRs but also to understand **why** a model makes a prediction by analyzing V gene usage, CDR3 features, motifs, and sequence patterns.

---

# Workflow

```
TCR sequences
      |
      v
Feature Engineering
      |
      +----------------+
      |                |
      v                v
V/J gene usage     CDR3 properties
(TRAV/TRBV)        (length + motifs)
      |
      v
Enrichment Analysis
      |
      v
Machine Learning Model
(Random Forest)
      |
      v
Prediction Probability
(Human / EBV / CMV)
      |
      v
Explainability
(SHAP feature contributions)
```

---

# Dataset

The analysis compares:

* Human baseline TCR repertoire
* EBV-specific TCRs
* CMV-specific TCRs

The baseline repertoire is used to identify enriched or depleted features.

---

# Feature Engineering

## V Gene Usage

Analyzed:

* TRAV genes
* TRBV genes

Examples:

```
TRBV7-9
TRBV27
TRAV13-1
```

---

## CDR3 Length

Analyzed:

* CDR3 alpha length
* CDR3 beta length

---

## Motif Extraction

Short amino acid patterns are extracted from CDR3 regions.

Examples:

```
TQYF
EQYF
KLFF
TGEL
```

These motifs are analyzed for enrichment.

---

# Enrichment Analysis

Features are compared between antigen-specific TCRs and human baseline.

The enrichment score:

```
log2(enrichment)
```

Positive values:

```
Feature enriched in antigen-specific TCRs
```

Negative values:

```
Feature depleted
```

---

# Machine Learning Model

## Model

Random Forest classifier

Classes:

```
0 = Human
1 = EBV
2 = CMV
```

The model outputs prediction probabilities:

Example:

```
Human : 0.05
EBV   : 0.00
CMV   : 0.95
```

---

# Model Performance

Example:

```
ROC-AUC (macro): 0.788
Accuracy: 0.88
```

---

# Explainable AI (SHAP)

SHAP analysis explains individual predictions by showing how each feature contributes to the final decision.

Example:

```
Feature:
TRAV10

Contribution:
Positive

Effect:
Pushes prediction toward a class
```

This converts the model from a black box into an interpretable biological model.

---

# Key Biological Insights

The analysis suggests that TCR specificity is influenced by:

* V gene usage bias
* CDR3 length preference
* Conserved sequence motifs
* TRAV/TRBV pairing patterns

Specificity is not determined by CDR3 sequence alone.

---

# Project Structure

```
CMV_TCR_TRBV_CDR3LengthStudy/

├── enrichment_CMV.py
├── enrichment_EBV.py
├── Pairing_CMV.py
├── Pairing_EBV.py
├── training.py
└── README.md
```

---

# Future Improvements

* Increase training dataset size
* Improve EBV class balance
* Add SHAP visualization plots
* Test additional ML models
* Integrate structural TCR-pMHC information

---

# Conclusion

This project demonstrates an explainable AI framework for TCR antigen specificity prediction by combining:

* Immunological feature discovery
* Machine learning prediction
* Model interpretability

The goal is to understand not only **which TCRs recognize an antigen**, but also **why**.
