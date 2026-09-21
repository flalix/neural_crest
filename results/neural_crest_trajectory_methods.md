# Neural crest — scFates trajectory analysis: methods and deviations

**Subject.** The neural crest–derived compartment of the developing human heart: the
Schwann cell precursor pool and its glial, chromaffin and autonomic-neuron
derivatives (SCP_1–5, My_SC, Bridge_state, Chrom_C, Aut_Neu_1–2).

**Data and published analysis being reproduced.**
Lázár, E. et al. *Spatiotemporal gene expression and cellular dynamics of the
developing human heart.* Nature Genetics 57, 2756-2771 (2025).
doi:10.1038/s41588-025-02352-6 — 36 hearts, PCW 5.5-14. The innervation trajectory
is that paper's Fig. 4c,d. Object analysed: `innervation.rds` (Mendeley
`w65jtfsvpr`); trajectory recipe: `06_3_HDCA_heart_scFates_innervation.py`
(Zenodo 15912657, github.com/rmauron/HDCA_heart_dev).

Species/stage confirmed on the object itself: human symbols (1,591 all-caps vs 27
title-case in the first 2,000 genes), `age` 5.5-14 PCW, 21 `TenX*` samples — matching
the published cohort.

**Conceptual reference for the lineage.**
Soldatov, R. et al. *Spatiotemporal structure of cell fate decisions in murine neural
crest.* Science 364, eaas9536 (2019). doi:10.1126/science.aas9536 — mouse E8.5-E10.5,
Wnt1/Sox10 lineage tracing; source of the SCP -> bridge -> chromaffin model and of
the *Htr3a* bridge marker tested in §10. Data: GSE129114,
pklab.med.harvard.edu/ruslan/neural.crest.html. This is a different dataset and
species; nothing here is fitted to it.

**Other compartments of the same atlas.** The cardiomyocyte, endothelial and
fibroblast compartments were fitted with the same pipeline (§6). They belong to the
same published atlas but are not neural crest derived, so they appear here only as
evidence for which methodological choices generalise.

## 1. Data provenance and the reconstructed selection

`innervation.rds` (53.7 MB, 1,106 cells x 31,136 genes) was the only input required.
`innervation_selection.rds` — the object `06_3` loads — is **not deposited**; the
deposit README states derived objects can be recreated from the shared ones.

It was reconstructed exactly rather than approximated: the `cell.types` column
contains nine clusters named `IN_excl_1`...`IN_excl_9` (236 cells). Dropping them
leaves **870 cells across ten types** — SCP_1-5, My_SC, Bridge_state, Chrom_C,
Aut_Neu_1-2 — which is `innervation_selection`. Cell counts per type and the
`clusters_subset` id each maps to are in [innervation_celltypes.md]({{artifact:6a76a9c3-f26a-4924-9970-de096c8701dd}}).

Context objects (not neural crest): `cardiomyocytes02.rds`, `endothelial.rds`,
`fibroblasts.rds`, each filtered on the same `*_excl_*` convention.

## 2. Seurat -> h5ad conversion

`06_2` used SeuratDisk (`SaveH5Seurat` + `Convert`). SeuratDisk is unmaintained and
its `Convert()` fails on Seurat v5 assay layouts, so conversion goes through a
Matrix Market bundle instead: [convert_seurat_to_h5ad.R]({{artifact:df55edf0-5f12-4d37-b837-385816293ca9}}) writes
`counts.mtx` + `barcodes.csv` + `features.csv` + `metadata.csv` + one CSV per
reduction, and [assemble_h5ad.py]({{artifact:58a7a918-10a9-41c8-be2d-229574b7119d}}) reassembles it, asserting row-order
agreement rather than trusting it. Counts land in both `X` and `layers['counts']`;
Seurat reductions become `obsm['X_<name>']`. All seven reductions carried through
(`pca`, `umap_df`, `harmony_oi`, `umap_oi`, `harmony_subset`, `umap_subset`,
`umap_subset_3d`), integer counts preserved.

## 3. Preprocessing

As `06_3`: `normalize_total`, `log1p(base=10)`, 5,000 `cell_ranger` HVGs, `.raw` set
before subsetting, scale, PCA, 25-neighbour graph. Implemented in
[run_scfates.py]({{artifact:525b7fef-84e5-4156-82c3-7c04a9cc9d25}}) (`preprocess`).

**Deviation:** in `06_3`, `X` is the *spliced* matrix (`adata = spliced_data.copy()`).
Here `X` is total counts, since the spliced matrices are not public. The embedding
therefore differs slightly from the published Fig. 4D and results are not numerically
identical.

## 4. Root choice — the main methodological finding

`06_3` roots by marker gene: `scf.tl.root(adata, "PENK")`, which selects the node of
highest mean PENK. **PENK is a terminal marker here, not a progenitor marker:**

| type | n | PENK mean | % expressing |
|---|---|---|---|
| Chrom_C | 41 | 0.88 | 65.9 |
| Aut_Neu_1 | 145 | 0.13 | 15.9 |
| SCP_1 | 168 | 0.02 | 2.4 |
| Bridge_state | 34 | 0.00 | 0.0 |

So the root lands on **node 12, an internal node**, splitting the curve into two
branches radiating from the chromaffin population:

```
seg 1:  12 -> 0   d = 96.9   all 553 SCPs, plus neurons between root and SCP tip
seg 2:  12 -> 2   d = 34.5   the neuronal tip
```

Pseudotime then measures distance from a midpoint in *both* directions — it runs
progenitor-ward on segment 1 — and the neuronal types are split across both
branches (Aut_Neu_1 49/96, Chrom_C 23/18). This is a property of the published
recipe, not of the reconstruction: all four candidate cell selections tested put the
PENK root on an internal node.

**Resolution:** refit rooted at the tip whose 30 nearest cells are richest in SCPs
(27% SCP at tip 0 vs 0% at tip 2). That gives a single segment spanning the full
graph (96.9 + 34.5 = 131.5 units; observed cell pseudotime range 129.8) and the
ordering SCP -> Chrom_C -> Aut_Neu. Both fits are reported in
[innervation_pseudotime.csv]({{artifact:42e8cc6b-d2a0-4fd3-ab32-146308f7918e}}); the PENK version is retained for comparability
with the published figure.

## 5. Topology: curve vs tree for the neural crest lineage

Lázár et al. report (Fig. 4c,d) "a fork-like transition from early SCPs towards two
parallel trajectories of neuronal-chromaffin and glial cell lineages" — i.e. the fork
is **neuronal-chromaffin vs glial**, and chromaffin is grouped *with* the neurons.

On `X_harmony_subset` (10 dims, `ppt`, lambda=100) that fork is recovered. At both 10
and 15 nodes the tree has 3 tips and 1 fork, with one terminal branch taking the
entire neuronal-chromaffin lineage and the other arm staying SCP:

| branch (15 nodes) | contents |
|---|---|
| seg 2 (terminal) | Bridge_state 34/34, Chrom_C 41/41, Aut_Neu_1 145/145, Aut_Neu_2 69/69, SCP_3 38, My_SC 12 |
| seg 1 (root) | SCP_5 191, SCP_4 56, SCP_2 36, My_SC 16, SCP_1 13 |
| seg 3 (terminal) | SCP_1 99, SCP_2 22 |

So the published topology holds: 100% of the bridge, chromaffin and autonomic-neuron
cells occupy one branch, separated from the SCP arm at a single fork.

Two things the data does **not** resolve, neither of them claimed by the paper:

- **Chromaffin does not separate from neuronal.** No configuration splits Chrom_C
  from Aut_Neu_1/2 into distinct branches — consistent with the paper treating them
  as one "neuronal-chromaffin" trajectory, but it means the sympathoadrenal step
  cannot be read off this topology (see §10).
- **The glial terminal is not resolved.** My_SC (28 cells) splits 16/12 between the
  SCP arm and the neuronal branch rather than forming its own tip, so the glial arm
  of the published fork is represented by the SCP branch generally, not by a
  My_SC-specific branch.

`X_pca` with `ndims_rep=2` — the recipe's setting — behaves much worse: 4 tips,
2 forks, and every fate fragmented across segments (Chrom_C 10/11/0/18, Aut_Neu_1
9/13/2/103). `ndims_rep=10` on `X_pca` fails outright with `ZeroDivisionError` from
empty nodes. Root purity on the harmony embedding is 1.00 vs 0.77-0.97 on PCA. Full
sweep: [innervation_topology_sweep.csv]({{artifact:8ae76e88-5985-4c21-9abe-0196d51ae583}}).

The single tip-rooted curve (§4) is retained for the pseudotime and gene-association
work because it gives one monotone axis over all 870 cells; the tree above is the
topology result.

For the three larger compartments (§6) trees were likewise used.

## 6. Other compartments (non–neural crest, context only)

| compartment | cells | hvgs | embedding | ndims | method | nodes | ppt_lambda | root_node | root_is_tip | tips | forks | branches |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cardiomyocytes | 20366 | 5000 | X_harmony_subset | 10 | ppt (SimplePPT) | 20 | 100 | 19 | 1 | 8 | 5 | 12 |
| Endothelium | 14473 | 5000 | X_harmony_subset | 10 | ppt (SimplePPT) | 20 | 100 | 7 | 1 | 8 | 4 | 11 |
| Fibroblasts | 19843 | 5000 | X_harmony_subset | 10 | ppt (SimplePPT) | 20 | 100 | 19 | 0 | 6 | 4 | 9 |

Roots were chosen from progenitor annotations, not markers. Cardiomyocytes and
endothelium root at progenitor-rich tips (progenitor fraction 1.00).
**Fibroblasts cannot be tip-rooted:** all six tips are differentiated states (VIC,
Adv_FB_1, Peric_MC, Valve_MC_2, PDE4C-high) and the EPDC progenitors sit in the tree
interior, so the tip picker returned 0.00 and would have silently ordered the lineage
backwards from valve interstitial cells. Rerooting at the node carrying the most
EPDC/Prol_FB assignment weight (node 19, internal) gives EPDC_2 (1.9) -> EPDC_1 (4.0)
-> adventitial/inflammatory (8-10) -> interstitial (20-25) -> valve mesenchyme
(39-41) -> VIC (42.7). An internal root is appropriate on a *tree*, where `seg`
identifies which path a cell is on, and inappropriate on a *curve*, where it merely
conflates two directions.

## 7. Pseudotime uncertainty

| compartment | n_map | t_range | t_sd_median | t_sd_p90 | t_sd_max | pct_cells_t_sd_over_5pct_range |
|---|---|---|---|---|---|---|
| Cardiomyocytes | 100 | 46.6 | 0 | 0.58 | 19.14 | 2.3 |
| Endothelium | 100 | 68 | 0 | 0.88 | 24.22 | 2.4 |
| Fibroblasts | 100 | 45.8 | 0 | 1.08 | 16.01 | 4.7 |

Innervation (single curve): median `t_sd` 0.95 on a 129.8-unit axis (0.7% of range),
95th percentile 2.32, max 2.90; no cell exceeds 5% of curve length. `n_map=100` is
therefore not buying much for a curve — the ordering is tightly determined *given*
the curve, and what is uncertain is the topology, which replication cannot address.
On the trees it does earn its cost: medians are 0 (14-24% of cells map
deterministically) but maxima reach 16-24 units on ranges of 46-68, with 2.3-4.7% of
cells above 5% of range — cells near forks whose branch assignment varies between
mappings. Per-type values: [innervation_tsd_by_celltype.csv]({{artifact:7db6e2db-1c85-48f2-a3e2-7f67cef10a00}}).

## 8. Genes associated with pseudotime

`scf.tl.test_association(fdr_cut=1e-4, A_cut=0.3)` then `scf.tl.fit`, both fitting
GAMs in R via rpy2 (R 4.5.3, mgcv 1.9.4). Significant genes were clustered into six
modules by KMeans on z-scored fitted profiles ([trends_figure.py]({{artifact:89cc06df-1563-4d34-9014-527c030fb75f}})).

| compartment | genes tested | associated | module sizes |
|---|---|---|---|
| Innervation | 5,000 | 1,725 | 248, 279, 300, 75, 632, 191 |
| Cardiomyocytes | 5,000 | 204 | 12, 8, 64, 83, 13, 24 |
| Endothelium | 5,000 | 1,310 | 332, 119, 179, 151, 356, 173 |
| Fibroblasts | 5,000 | 769 | 119, 114, 127, 111, 120, 178 |

Binning for the module heatmaps is **per segment**, not pooled over pseudotime: on a
branching tree, cells on different branches share `t` values, so pooling averages
distinct fates together and produces spurious spikes in the module profiles.

Cardiomyocytes yield only 204 associated genes despite being the
largest compartment, and those that pass are the contractile set (*MYH6/7*, *MYL2*,
*NPPA*, *TTN*, *RYR2*) — consistent with the tree's axis capturing chamber identity
rather than a maturation gradient.

## 9. Technical signal to exclude before interpreting

Recurring non-biological modules: heat-shock (*HSPA1B*, *HSP90AA1*, *HMGB2*,
*UBE2C*) and immediate-early genes (*FOS*, *FOSB*, *JUNB*) form their own modules in
innervation (M1, 248 genes), fibroblasts and cardiomyocytes. These are dissociation
artifacts that track pseudotime.

**Ambient cardiomyocyte mRNA (innervation M4, 75 genes).** *METTL7A*, *TLX2*, *CKM*,
*TNNC1*, *PLN*, *CSRP3*, *MYL7*, *TNNT2* — sarcomere and muscle genes in a neural
crest compartment. They are droplet background, not biology: cells in the top M4
decile have a median 2,281 UMIs vs 12,385 overall, and M4 score correlates -0.51 with
log10 nCount. The load is uneven and falls hardest on the chromaffin population —
61% of Chrom_C cells sit in the top M4 decile (median 4,009 UMIs), vs 1-4% of SCP
cells. Exclude the 75 genes rather than the cells; proper correction (SoupX,
CellBender) needs the unfiltered droplet matrices from Mendeley part 1.
Details: [innervation_M4_ambient_check.csv](/home/flavio/.claude-science/orgs/abb023bb-bd59-418e-bf97-0ebfba525efa/artifacts/proj_c74c86fe51af/9e531553-0a3e-485a-8621-16bb90c38fe5/vb44bf866_innervation_M4_ambient_check.csv).

Technical total for the innervation compartment: 323 of 1,725 associated genes
(M1 stress 248 + M4 ambient 75). The neural crest biology is M2/M3 (glial and matrix:
*S100B*, *COL4A1*, *TGFBR2*) and M5/M6 (neuronal: *DLGAP2*, *GNAO1*, *NTRK3*).

At cell-type level the same signal identifies one cluster as technical
([innervation_qc_by_celltype.csv]({{artifact:2aceacef-bb9d-4a01-9d0e-4824757b7abc}})):

| cell type | n | percent_hsp | percent_mito | nCount |
|---|---|---|---|---|
| My_SC | 28 | 4.31 | 6.46 | 12228 |
| SCP_3 | 38 | 2.17 | 8.41 | 10643 |
| SCP_1 | 168 | 0.92 | 5.41 | 13443 |
| SCP_5 | 195 | 0.81 | 7 | 16486 |

`SCP_3` (38 cells) has 2.4x the heat-shock fraction of the other SCPs, the highest
mitochondrial fraction of all ten types, and a top-ranked marker set that is entirely
stress response — it should be excluded from pseudotime interpretation. `My_SC` has
the highest %hsp in the compartment but a genuine myelin identity (*MBP*, *PMP22*,
*S100B*) on 28 cells.

## 10. Sympathoadrenal transition — what the data supports

Marker gradients along the reconstruction match the SCP -> bridge -> chromaffin model:
*ASCL1* peaks in Bridge_state (0.95, ~10x the SCPs), glial markers decline stepwise
(*PLP1* 1.23 -> 0.75 -> 0.17; *SOX10* 0.51 -> 0.33 -> 0.09), and the catecholamine
program appears only in chromaffin cells (*TH* 0.66, *CHGA* 0.81, *CHGB* 1.33).

Two caveats. *HTR3A* (Soldatov et al.'s murine bridge marker) is effectively absent
(0.01) and does not transfer to human fetal heart here. And the tree places
Bridge_state, Chrom_C and both Aut_Neu populations on a single branch (§5) while the
tip-rooted curve puts Bridge_state at the low extreme rather than as an intermediate
— so the *ordering* within the neuronal-chromaffin lineage is not resolved, even
though the lineage itself separates cleanly from the SCP arm. With 41 chromaffin and
34 bridge cells that is expected, and the published analysis makes no claim beyond
the joint neuronal-chromaffin trajectory. Independent support for the transition here
is the marker gradient above, not the topology.

## 11. RNA velocity: not possible on public data

scVelo requires intron-level quantification. `06_1` reads `shoji/Spliced` and
`shoji/Unspliced` from per-sample `.h5` files produced by a collaborator's pipeline;
those files are not deposited. Cell Ranger and Space Ranger outputs in the deposits
carry a single count matrix (`matrix/` group only, no `shoji` group), and intron
inclusion is folded into it rather than split — it cannot be unmixed after the fact.

Routes, in order of cost: request `data/h5_files/*.h5` (or the spliced/unspliced RDS
pair) from the authors; or apply to EGA (`EGAS50000001029`) and re-quantify BAMs with
velocyto or `kallisto|bustools nac`. [build_velocity_adata.py]({{artifact:5205a7df-404f-43c6-b539-0f7cec323e0d}}) provides
`attach_spliced_unspliced()` (name-matched, with a barcode-overlap guard) and
`build_from_shoji()` for whenever the counts arrive; [check_velocity_inputs.py]({{artifact:e863b43d-f878-43eb-aec6-63daa9b46720}})
pre-flights an object before the dynamical model.

## 12. Version deviations from the published environments

**scVelo** (`scVelo.yml` pins scvelo 0.3.0, pandas 2.0.3). On scvelo 0.3.4:
`scv.pp.filter_genes_dispersion` and `scv.pp.log1p` no longer exist, and
`filter_and_normalize` no longer accepts `n_top_genes`. Replacement built on scanpy
equivalents in `preprocess_for_dynamical()`. On pandas 3, `recover_dynamics` raises
`TypeError` from `pandas.unique()` losing list support (GH#52986);
[scvelo_pandas_patch.py]({{artifact:db5709f3-371c-4cc8-adc1-0f93da81e22d}}) carries the pipeline as far as `velocity_graph` but not
`latent_time`. The supported fix is `pandas<3`: verified clean on
**scvelo 0.3.4 + pandas 2.2.3 + anndata 0.11.4 + zarr<3**, Python 3.12. Note
anndata 0.13.x assumes pandas 3 and fails on pandas 2.2 with
`StringDtype.__init__() got an unexpected keyword argument 'na_value'`.

**scFates** (`scFates.yml` pins 1.0.6). On 1.2.5: `test_association` takes `fdr_cut`,
not `fdr`. `scf.tl.pseudotime` raises `IndexError` in its colour-assignment step when
a milestone receives no cells (30-node trees on these compartments); every numeric
output is assigned before that line, so it is safe to wrap — 20 nodes avoids it
entirely. `scf.pl.trends(annot="seg")` raises `UnboundLocalError` on a single-segment
fit, because the colour map is only built when more than one segment exists; use
`annot="milestones"`.

**rpy2 / R.** scFates fits its GAMs in R, and needs **two** independent things:
`R_HOME` for rpy2 (it dlopens `libR.so` at import), and the `R` **binary on `PATH`**
for scFates itself (`importeR` decides availability with `shutil.which("R")`).
Setting only `R_HOME` gives working rpy2 and a scFates that reports
"R installation is necessary". Both are cached at import and cannot be fixed
afterwards: `openrlib.R_HOME` stores `None` silently and raises later in `initr()` as
`Unable to determine R_HOME`, while scFates freezes its availability flags at module
import. Re-importing rpy2 over an already-initialised embedded R produces
`R was initialized outside of rpy2` and inconsistent behaviour between modules — a
kernel restart with both variables already set is the only reliable fix.
`run_scfates.py` resolves and repairs both before importing scFates.

## 13. Artifacts

Figures: [innervation_trajectory.png]({{artifact:cde07fcc-f8f2-4906-bc28-25988d34eced}}), [pseudotime_uncertainty.png]({{artifact:c587f3e0-6de6-4261-b510-8bc222c1be95}}),
[trends_innervation.png]({{artifact:91026804-804d-4a23-b307-cb6b6cd9699d}}), [trajectory_cardiomyocytes02.png]({{artifact:a0ddc43c-a086-477a-8723-5713fc0bce9d}}),
[trajectory_endothelial.png]({{artifact:3e5de5e2-d197-4e85-b06b-da765791d2c6}}), [trajectory_fibroblasts.png]({{artifact:56bf1f90-3489-4183-9909-d6a43b626166}}),
[trends_cardiomyocytes02.png]({{artifact:f76f8d43-acf8-48c3-9cb3-194a9226aa9a}}), [trends_endothelial.png]({{artifact:0780329c-69fc-403b-9983-210e6b463a60}}),
[trends_fibroblasts.png]({{artifact:ee9c7cdf-4142-45ca-8ab9-a18584d487f4}})

Tables: [roots_and_parameters.csv]({{artifact:7f4e9219-a79c-438e-8569-7fb53aeb1a3d}}), [innervation_pseudotime.csv]({{artifact:42e8cc6b-d2a0-4fd3-ab32-146308f7918e}}),
[cardiomyocytes02_pseudotime.csv]({{artifact:6a4f46e3-8867-4692-bcca-34401b1e3a66}}), [endothelial_pseudotime.csv]({{artifact:7c1bd525-afa3-4d2f-9d07-5f04aad5fca0}}),
[fibroblasts_pseudotime.csv]({{artifact:b058e3f9-f082-4972-aa37-a1dff1661561}}), [innervation_gene_modules.csv]({{artifact:03461c2e-4257-4b96-99b6-2d6c7a1a5a9c}}),
[cardiomyocytes02_gene_modules.csv]({{artifact:3961115b-47a1-4b03-8fbd-547c53d12869}}), [endothelial_gene_modules.csv]({{artifact:10bb6610-55a8-41ce-968d-bbee2567644c}}),
[fibroblasts_gene_modules.csv]({{artifact:2d058aa7-0ada-4f49-ba67-f30b4f08bca6}}), [innervation_association.csv]({{artifact:6255fe2e-072f-4fe5-820d-d1460af1bf75}}),
[cardiomyocytes02_association.csv]({{artifact:d3c375f7-adc0-434e-81ce-3feebcc17167}}), [endothelial_association.csv]({{artifact:36231347-bca9-488b-a428-68d760d6ce92}}),
[fibroblasts_association.csv]({{artifact:edffe0ae-55c9-44b2-add0-a6bc98d66b83}}), [innervation_celltype_markers.csv]({{artifact:2ad9e9ae-c378-46f8-8d97-2179bfbde6b1}}),
[innervation_qc_by_celltype.csv]({{artifact:2aceacef-bb9d-4a01-9d0e-4824757b7abc}}), [innervation_tsd_by_celltype.csv]({{artifact:7db6e2db-1c85-48f2-a3e2-7f67cef10a00}}),
[innervation_topology_sweep.csv]({{artifact:8ae76e88-5985-4c21-9abe-0196d51ae583}}), [scfates_pseudotime_timing.csv]({{artifact:ff0e823a-3ac0-42f9-8bb7-fe45ab21463c}})

Code: [convert_seurat_to_h5ad.R]({{artifact:df55edf0-5f12-4d37-b837-385816293ca9}}), [assemble_h5ad.py]({{artifact:58a7a918-10a9-41c8-be2d-229574b7119d}}), [run_scfates.py]({{artifact:525b7fef-84e5-4156-82c3-7c04a9cc9d25}}),
[trends_figure.py]({{artifact:89cc06df-1563-4d34-9014-527c030fb75f}}), [build_velocity_adata.py]({{artifact:5205a7df-404f-43c6-b539-0f7cec323e0d}}),
[check_velocity_inputs.py]({{artifact:e863b43d-f878-43eb-aec6-63daa9b46720}}), [scvelo_pandas_patch.py]({{artifact:db5709f3-371c-4cc8-adc1-0f93da81e22d}})

Fitted object: [innervation_selection_scfates.h5ad]({{artifact:fd40245c-49ab-481b-a389-e70e2e6a9405}}) (870 x 5,000, both roots in
`obs`; `uns['epg']` dropped — elpigraph's ragged `Edges` list cannot be written to
h5ad)
