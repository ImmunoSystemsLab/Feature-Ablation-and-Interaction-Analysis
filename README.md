# 🧬 Explainable AI Framework for TCR Antigen Specificity

## Overview

This project develops an explainable machine learning framework to predict T cell receptor (TCR) antigen specificity and identify the biological features responsible for recognition.

The main idea is to move beyond simple classification and answer:

**"Why does a TCR recognize a specific antigen?"**

The pipeline integrates:

* Germline V gene usage
* CDR3 sequence features
* Motif enrichment
* αβ chain pairing patterns
* Machine learning prediction
* SHAP-based explanation

---

# Biological Hypothesis

Antigen specificity is not determined by a single TCR feature.

Instead, antigen recognition emerges from combinations of:

* V gene usage bias
* V-J pairing
* CDR3 length constraints
* Conserved amino acid motifs

These combined features create an antigen-associated TCR repertoire fingerprint.

---

# Analysis Workflow

```
TCR repertoire
       |
       v
Feature extraction
       |
       +----------------------+
       |                      |
       v                      v
V gene usage            CDR3 properties
(TRAV/TRBV)             length + motifs
       |
       v
Pairing analysis
(TRAV-TRBV interactions)
       |
       v
Machine learning model
(Random Forest)
       |
       v
Prediction probability
       |
       v
SHAP explanation
```

---

# Feature Engineering

## 1. Germline V Gene Usage

Analyzed:

* TRAV genes
* TRBV genes

Example:

```
TRBV7-9
TRBV27
TRAV13-1
```

The goal is to identify whether certain germline regions are preferentially used during antigen-specific responses.

---

## 2. CDR3 Length Analysis

CDR3 regions were analyzed because they directly contribute to peptide recognition.

Features:

* CDR3 alpha length
* CDR3 beta length

Specific length preferences may indicate structural constraints for antigen binding.

---

## 3. Motif Discovery

Short amino acid patterns were extracted from CDR3 sequences.

Examples:

```
TQYF
EQYF
KLFF
TGEL
```

Enriched motifs may represent conserved recognition patterns.

---

# Pairing Analysis and Signal Enhancement

A key analysis was testing whether combining features improves antigen-associated signals.

Instead of asking:

"Is this V gene enriched?"

the analysis asks:

"Does this V gene become more informative when combined with other TCR features?"

---

## Example: TRBV7-9 Feature Combination

When analyzing TRBV7-9 alone:

```
TRBV7-9 only

log2 enrichment ≈ 0.5
```

The signal was relatively weak.

After adding biological context:

```
TRBV7-9 + Motif

log2 enrichment ≈ 1.5
```

The association became stronger.

When combined with additional TCR context:

```
TRBV7-9 + TRAV partner

log2 enrichment ≈ 2.4
```

and:

```
TRBV7-9 + CDR3 length

log2 enrichment ≈ 2.4
```

---

## Interpretation

This suggests that:

A single gene usage pattern is not enough to define specificity.

The antigen recognition signal becomes clearer when considering the complete TCR architecture:

* germline encoded regions
* junctional diversity
* sequence motifs
* chain pairing

---

# Machine Learning Model

## Model

Random Forest classifier

Classes:

```
0 = Human baseline
1 = EBV-specific
2 = CMV-specific
```

The model predicts:

```
P(TCR specificity | features)
```

Example:

```
Human : 0.05
EBV   : 0.00
CMV   : 0.95
```

---

# Explainable AI (SHAP)

SHAP analysis identifies which features influenced each prediction.

Example:

```
Feature:
TRAV10

Effect:
Increases probability of Human-like repertoire
```

or

```
Feature:
CDR3 length + motif combination

Effect:
Supports antigen-specific prediction
```

This converts the model from a black box into a biologically interpretable system.

---

# Key Biological Insights

The analysis suggests:

1. Antigen-specific repertoires have characteristic feature combinations.

2. Germline V genes contribute to specificity but are not sufficient alone.

3. Pairing and sequence context can amplify weak signals.

4. Motifs may represent convergent recognition patterns.

5. TCR specificity is an emergent property of multiple interacting features.

---

# Project Structure

```
CMV_TCR_TRBV_CDR3LengthStudy/

├── enrichment_CMV.py
├── enrichment_EBV.py
├── Pairing_CMV.py
├── Pairing_EBV.py
├── training.py
├── SHAP_explanation.py
└── README.md
```

---

# Future Directions

* Increase dataset size
* Improve class balance
* Add SHAP visualization plots
* Test additional interpretable models
* Integrate structural TCR-pMHC information

---

# Conclusion

This project demonstrates an explainable AI approach for understanding TCR antigen specificity.

The goal is not only to predict:

**"Which TCR recognizes an antigen?"**

but also:

**"Which biological features create the recognition fingerprint?"**
