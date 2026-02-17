# Botello-2026-ribosome-aging


## Contents
Each script and pipeline in the `scripts` folder is explained below.

### Quantification of nucleolar core/rim enrichment of ribosome signal
Images were analyzed with the `nucleolar_core_rim.cppipe` using CellProfiler
(v4.2.6) and cellpose (v2.3.2).  The distribution of of ribosomes were measured
with the MeasureObjectIntensityDistribution module and saved as csv files
for analysis with `nucleolar_core_rim.ipynb`.

The jupyter notebook additionally requires matplotlib, pandas and seaborn.  The
radial distribution is normalized by the area of each section to generate an
enrichment measurement representing the increase or decrease in intensity relative
to a uniform distribution.  While the code can accommodate any number of bins,
3 bins were used to define the core and rim region of each nucleolus.
