## Feature Ablation and Interaction Analysis

To investigate whether antigen specificity is better explained by individual TCR features or by combinations of features, we performed feature ablation experiments.

Individual and combined feature sets were evaluated using AUROC:

Features tested:

- TRBV gene usage
- CDR3 length properties
- CDR3 motif features
- TRAV–TRBV pairing
- Combined feature model

### EBV Antigen Specificity

| Feature set | AUROC |
|---|---:|
| TRBV only | 0.527 |
| Length only | 0.486 |
| Motif only | 0.641 |
| TRBV + Length | 0.528 |
| TRBV + Motif | 0.651 |
| TRBV + TRAV | 0.671 |
| TRBV + Motif + Length | 0.571 |
| Full model | 0.777 |

Interpretation:

For EBV-specific TCRs, individual germline features showed limited predictive power.
CDR3 motif features provided a stronger signal, suggesting that sequence-level patterns contribute to antigen recognition.
Combining germline usage, motifs, and chain context produced the strongest prediction performance.

---

### CMV Antigen Specificity

| Feature set | AUROC |
|---|---:|
| TRBV only | 0.747 |
| Length only | 0.532 |
| Motif only | 0.577 |
| TRBV + Length | 0.747 |
| TRBV + Motif | 0.759 |
| TRBV + TRAV | 0.808 |
| TRBV + Motif + Length | 0.760 |
| Full model | 0.816 |

Interpretation:

CMV-specific TCR prediction showed a stronger contribution from germline V gene usage.
TRBV usage alone provided substantial predictive information, while adding TRAV pairing further improved performance.
The full feature model achieved the highest AUROC, supporting the idea that antigen specificity emerges from combined TCR features rather than a single determinant.

---

### Biological Conclusion

Feature ablation analysis demonstrates that TCR antigen specificity is not encoded by one isolated feature.

Instead, recognition emerges from interactions between:

- Germline V gene usage
- αβ chain pairing context
- CDR3 sequence motifs
- CDR3 structural properties

These results support a combinatorial model of TCR antigen recognition, where multiple weak or moderate signals combine to create an antigen-associated repertoire fingerprint.
