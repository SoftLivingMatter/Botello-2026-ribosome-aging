# Botello-2026-ribosome-aging
A collection of scripts for "Ribosome Molecular Aging Shapes Translation Dynamics".
Additional information available upon request.


## Contents
Each script and pipeline in the `scripts` folder is explained below.

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
sites.

### Occupancy dynamics and identification of dropoff codons
Display of relative ribosome occupancy and detection of the dropoff codon is
contained in the `occupancy_dynamics.ipynb` notebook.  pysam is required
to read the genome sequence and alignment files.  Biopython Seq objects are
used for performing codon translation and scipy and sklearn are used for fitting
the occupancy profiles to logistic functions.

### Determination of age-dependent RNA modifications
pseudoU-BIDseq or 2'-O-methylation scores were further analyzed in R to assess
which sites were significantly changed upon ribosome aging, as presented in
`Modification_Analysis.Rmd`

### Permutation test of solvent accessibility at modified sites
Whether modified rRNA sites are more solvent-exposed than other sites of the same
base was tested with `perm_test_sasa_4col.py`.  The input is a csv file with four
columns of SASA values: 28S modified sites, 28S all other sites of that base, 18S
modified sites and 18S all other sites of that base.  Because each background
contains only the modifiable base, random draws are automatically base-matched.
For each rRNA, the mean or median SASA of the modified sites is compared to a null
distribution built from random draws (without replacement) of the same number of
sites from the background.  Mann-Whitney U p-values and AUC are reported as a
cross-check.  The script requires numpy, pandas and matplotlib, and optionally
scipy.  Results are written as csv, png and pdf files.
