# Botello-2026-ribosome-aging
A collection of scripts, notebooks and image-analysis pipelines for "Ribosome Molecular Aging Shapes Translation Dynamics".
Additional information available upon request.


## Contents
Each script, notebook and pipeline in the `scripts` folder is explained below.
Software versions are given where they were recorded; Python packages are
otherwise listed per section and were not version-pinned.

### Quantification of nucleolar core/rim enrichment of ribosome signal
Images were analyzed with the `nucleolar_core_rim.cppipe` using CellProfiler
(v4.2.6) and cellpose (v2.3.2).  The distribution of ribosomes were measured
with the MeasureObjectIntensityDistribution module and saved as csv files
for analysis with `nucleolar_core_rim.ipynb`.

The jupyter notebook additionally requires matplotlib, pandas and seaborn.  The
radial distribution is normalized by the area of each section to generate an
enrichment measurement representing the increase or decrease in intensity relative
to a uniform distribution.  While the code can accommodate any number of bins,
3 bins were used to define the core and rim region of each nucleolus.

### Quantification of cytoplasmic ribosome signal
Images were analyzed with the `cytoplasm_measurement.cppipe` using CellProfiler
(v4.2.6). Nuclei were identified from a blurred, Hoechst stain and cell boundaries
estimated by propagating the nuclei with the cytoplasmic stain.  The final
cytoplasmic region was determined by removing the nucleus from the cell boundary
and used to measure the intensity of stained ribosomes.  Resulting csv files were
analyzed in Prism.

### Ribosome pausing and occupancy analysis
Detection of age-dependent pause sites or polybasic codon runs and subsequent
occupancy measurements were performed in the jupyter notebook `pause_occupancy.ipynb`
using ribopy databases as inputs.  Each sample is extracted from the ribo files
using `ribopy dump coverage`.  The coverage, location and gene is read from the
dumped files and used to calculate pause scores or detect age-dependent pause
sites.  The notebook requires ribopy (for `ribopy dump`, run in its own
environment) and biopython, matplotlib, numpy, pandas, scipy, seaborn and
statsmodels.

### Occupancy dynamics and identification of dropoff codons
Display of relative ribosome occupancy and detection of the dropoff codon is
contained in the `occupancy_dynamics.ipynb` notebook.  pysam is required
to read the genome sequence and alignment files.  Biopython Seq objects are
used for performing codon translation and scipy and sklearn are used for fitting
the occupancy profiles to logistic functions.

### Determination of age-dependent RNA modifications
pseudoU-BIDseq or 2'-O-methylation scores were further analyzed in R to assess
which sites were significantly changed upon ribosome aging, as presented in
`Modification_Analysis.Rmd`.  For each site, an ANOVA across the time points
is run, p-values are corrected by Benjamini-Hochberg, and Tukey post-hoc tests
are used for pairwise comparisons.  The notebook requires R and the tidyverse.

### Permutation test of solvent accessibility at modified sites
Whether aging-related 8-oxo-guanine modified rRNA sites are more solvent-exposed
than other sites of the same base was tested with `perm_test_sasa_4col.py`.  The
input is a csv file with four columns of SASA values: 28S modified sites, 28S all
other sites of that base, 18S modified sites and 18S all other sites of that base.
Because each background contains only the modifiable base, random draws are
automatically base-matched.  For each rRNA, the mean or median SASA of the modified
sites is compared to a null distribution built from random draws (without
replacement) of the same number of sites from the background.  Mann-Whitney U
p-values and AUC are reported as a cross-check.  The script requires numpy, pandas
and matplotlib, and optionally scipy.  Results are written as csv, png and pdf
files.

### Statistical determination of molecular age-dependent protein interactions
Proteomic abundances from the RPL10A and RPS17 pulldowns (three biological
replicates each at T0, T24 and T72) were analyzed with the marimo notebook
`anova_trend.py`.  For each protein, differences among the three time points
are assessed by one-way ANOVA, requiring at least two non-missing values per
time point; proteins with zero variance across all samples are excluded.
P-values are adjusted for multiple testing within each pulldown using the
Benjamini-Hochberg method (q < 0.05).  The direction of change is assigned
from the group means: proteins are classified as increasing (+) if mean
abundance rose at each successive time point (T0 < T24 < T72), decreasing (-)
if it fell at each successive time point (T0 > T24 > T72), and otherwise as
having no monotonic trend (0).  The notebook also plots volcano plots and
abundance profiles of the top-ranked increasing and decreasing proteins.

The inputs are two excel files, `L10Aproteomics_08182026.xlsx` and
`S17proteomics_08182026.xlsx`, expected in the same directory as the notebook,
with the gene name (`Gene Name`) in the first column and columns `T0-1`...`T72-3` for each
replicate.  Results are written to `results.csv`.  Dependencies (pandas, scipy,
statsmodels, plotly, openpyxl and marimo) are declared inline, so
`uv run anova_trend.py` will set up the environment, or use
`marimo edit anova_trend.py` to explore interactively.

### Delta-delta G analysis of rRNA secondary structure stability
To estimate the impact of 8-oxo-guanine lesions on the secondary-structure
stability of human 28S and 18S rRNAs, G-to-U substitutions were introduced in
silico to approximate the altered base-pairing properties of 8-oxo-G.  This was
done with `scripts/rna_mutation_screen/g2u_sweep.py`.  For each guanosine, the
change in folding free energy (ddG) relative to a fixed reference structure is
calculated with the ViennaRNA Python bindings (v2.7.2, nearest-neighbor Turner
energy model; equivalent to `RNAeval`).  The wild-type free energy
dG(WT, S) is calculated once, every guanosine is then individually replaced with
uracil, and the mutant is re-evaluated against the same structure S:
ddG = dG(mutant, S) - dG(WT, S).  Positive values indicate destabilization of the
reference fold and unpaired positions are expected to be ~0.  The structure is
held fixed rather than re-folded because, in the assembled ribosome, rRNA folding
is constrained by protein contacts and tertiary interactions that a thermodynamic
re-fold would not capture.  Only canonical, nested secondary-structure energetics
are evaluated; tertiary contacts, protein-binding contributions and pseudoknots
are not accounted for.

```bash
cd scripts/rna_mutation_screen
python g2u_sweep.py --name 18S --fasta data/18S.fasta --structure data/18S.dbn --out output/18S_G-to-U_ddG.csv
python g2u_sweep.py --name 28S --fasta data/28S.fasta --structure data/28S.dbn --out output/28S_G-to-U_ddG.csv
```

The output csv has one row per mutation (position, whether it is paired, its
partner position and base, and ddG), with the input file hashes and ViennaRNA
version in the header.  The only dependency is ViennaRNA >= 2.6
(`pip install -r requirements.txt`).  The script stops with an error if a sequence
is not a single A/C/G/U record, if the structure contains anything but `.()` or
is unbalanced, or if the sequence and structure lengths differ.

The sequences and structures are not included and must be downloaded to a
`data/` folder in `scripts/rna_mutation_screen` as `18S.fasta`, `18S.dbn`,
`28S.fasta` and `28S.dbn`.  They were exported from RNAcentral for human 18S
([URS00025D0EC0/9606](https://rnacentral.org/rna/URS00025D0EC0/9606)) and 28S
([URS000075EC78/9606](https://rnacentral.org/rna/URS000075EC78/9606)) rRNA, using
"Copy as RNA" for the sequence and the dot-bracket copy button for the structure.
The structures are R2DT predictions built on the RiboVision SSU and LSU templates
and are pseudoknot-free.  The sha256 hashes of the files used were:

| File | sha256 |
|---|---|
| `18S.fasta` | `7c722d546899bb5764e7bfe26ebbda3249e453f33996f02e1bbaaa01bdf13ae2` |
| `18S.dbn` | `1bb7975e0f2772e39e9adb4f14ea1db11c91699efd5c429378ea615d06697e6b` |
| `28S.fasta` | `eed08754436750c1a8d616d8ad15610f399e1aace9963db39145f8d35f420815` |
| `28S.dbn` | `f4165a4d9bac1f3398074601a322775eafe47df98861aff06e3f94e1b34fc1c5` |
