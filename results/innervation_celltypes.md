# Innervation compartment — cell types

Reconstructed `innervation_selection`: 870 cells (nine `IN_excl_*` clusters removed from `innervation.rds`). Pseudotime from a 30-node ElPiGraph curve on PC1–2, rooted at the SCP tip, `n_map=100`, `seed=42`. Markers ranked by Wilcoxon test on `.raw` (all 31,136 genes), top six shown.

| cell type | n | pseudotime (SCP-tip root) | top markers |
|---|---:|---:|---|
| Bridge_state | 34 | 39.6 | *ASCL1*, *TBX2*, *LINC00682*, *TLX2*, *ZNF704*, *DPYD* |
| SCP_1 | 168 | 47.5 | *PLS3*, *MEF2C*, *POSTN*, *PTPRJ*, *ADAMTSL1*, *IL1RAPL2* |
| SCP_3 | 38 | 49.8 | *HSP90AB1*, *HSP90AA1*, *HSPA1B*, *JUND*, *CREB5*, *PEAK1* |
| SCP_5 | 195 | 49.8 | *SPARC*, *COL1A2*, *RELN*, *COL1A1*, *GFRA3*, *PLEKHA4* |
| My_SC | 28 | 50.4 | *COL3A1*, *PMP22*, *MBP*, *COL1A1*, *COL1A2*, *S100B* |
| SCP_4 | 72 | 53.9 | *TRPM3*, *NRXN3*, *MIR99AHG*, *OLFML2A*, *APP*, *PLP1* |
| SCP_2 | 80 | 54.0 | *LRRC7*, *PLEKHG1*, *MEF2C*, *SOX5*, *ANXA1*, *PTPRZ1* |
| Chrom_C | 41 | 78.1 | *DLK1*, *CHGB*, *FAM162B*, *NDUFA4L2*, *RAMP1*, *RGS5* |
| Aut_Neu_1 | 145 | 90.9 | *EML5*, *SYT1*, *ELAVL4*, *ARHGEF28*, *CADPS*, *SLIT3* |
| Aut_Neu_2 | 69 | 99.7 | *TUBB2B*, *TUBB3*, *TUBA1A*, *STMN2*, *TMSB10*, *TUBB* |

Rows ordered by mean pseudotime. `SCP_3`'s signature is entirely heat-shock (*HSP90AB1*, *HSP90AA1*, *HSPA1B*, *JUND*), i.e. dissociation stress rather than a cell identity.
