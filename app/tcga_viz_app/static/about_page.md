## About
This web portal accompanies the manuscript "Pan-cancer machine learning
predictors of tissue of origin and molecular subtype", currently available on
[bioRxiv](https://doi.org/10.1101/333914).  Briefly, we constructed classifiers
from annotated expression profiles of cancer samples to predict a cancer's
primary tissue of origin and molecular subtype (if available). The expression
data comes from The Cancer Genome Atlas (TCGA) and is available
[here](https://gdac.broadinstitute.org/).

Abbreviations of primary tumor types are explained
[here](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations).

### This portal
This portal provides interactive visualization of two kinds:

1. 3-dimensional embeddings of the expression profiles of ~10,000 samples
   spanning 33 primary cancer types and matched normals.  This visualization
   aids in understanding the relationships between cancers of different primary
   sites of origin and molecular subtypes.

2. Classification metrics of our primary site and molecular subtype predictors.

### Code availability
The code that runs this visualization and produced the results in the manuscript
is available on
[GitHub](https://github.com/TheJacksonLaboratory/tcga_subtype_classification).

### License
This application is free to use for academic and non-commercial use, and this
application and the code that powers it is subject to [this
license](/static/LICENSE.txt).
