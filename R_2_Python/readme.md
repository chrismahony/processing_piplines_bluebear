## Easily write out that data you need to convert Anndata obj to Seurat


1. Load what you need

```python

import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

amp2 = sc.read("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/AMP2_data/AMPII_annotated.h5ad")


```
