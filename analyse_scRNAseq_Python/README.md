# Analysing scRNAseq data using Python and Scanpy

1. To install Python modules using bluebear you must follow this: https://docs.bear.bham.ac.uk/bluebear/software/self_installs_python/#creating-a-virtual-environment-and-installing-a-python-module
!!! Do not just pip install in a python enivroment!!!

<br>

2. Brefily, open a terminal, navigate to where you are doing your analysis and load the versin of python you want to use:

```Bash

module load bear-apps/2022a
module load Python/3.10.4-GCCcore-11.3.0

```

3. Then create a virtual eniroment and activate it

```Bash

python3 -m venv --system-site-packages my-virtual-env-${BB_CPU}
source my-virtual-env-${BB_CPU}/bin/activate

```


4. Now install packages:

```Python

pip install scanpy
pip install bbknn
pip install scikit-image
pip install git+https://github.com/chrismahony/mcseqy.git

```


5. Next start a Jupitylab session here: https://portal.bear.bham.ac.uk/pun/sys/dashboard/batch_connect/sessions
!!!Make sure you select the same version of Python that you loaded above, in this case it was 3.10.4!!!

6. Now you can load modules that you installed in your virtual eniroment:

```Python

import os
from pathlib import Path
import sys
node_type = os.getenv('BB_CPU')
venv_dir = f'/rds/my_path/my-virtual-env-icelake'  # edit this line to match the venv directory format
venv_site_pkgs = Path(venv_dir) / 'lib' / f'python{sys.version_info.major}.{sys.version_info.minor}' / 'site-packages'
if venv_site_pkgs.exists():
    sys.path.insert(0, str(venv_site_pkgs))
else:
    print(f"Path '{venv_site_pkgs}' not found. Check that it exists and/or that it exists for node-type '{node_type}'.")


import scanpy as sc


```



7. Next defin these function that will make processing easier and reduce the amout of work a bit:


```Python

import scanpy as sc
from pathlib import Path

def load_10x_h5_to_anndata(file_paths, sample_names=None, combine=False):
    """
    Load 10x Genomics .h5 files into AnnData objects and optionally combine them.
    
    Parameters:
    - file_paths: List of strings or Path objects, each path pointing to a 10x .h5 file.
    - sample_names: List of sample names, one per file path. Defaults to the file name if not provided.
    - combine: Boolean, if True, combines all matrices into a single AnnData object.
    
    Returns:
    - If combine is False, returns a list of AnnData objects.
    - If combine is True, returns a single combined AnnData object with sample names.
    """
    
    if sample_names is None:
        sample_names = [Path(fp).stem for fp in file_paths]
    
    if len(sample_names) != len(file_paths):
        raise ValueError("The length of sample_names must match the number of file_paths.")

    adata_list = []
    
    for file_path, sample_name in zip(file_paths, sample_names):
        # Read the .h5 file into an AnnData object
        adata = sc.read_10x_h5(file_path)
        
        # Ensure gene names are unique
        adata.var_names_make_unique()
        
        # Add sample name as metadata
        adata.obs['sample'] = sample_name
        
        # Append AnnData object to list
        adata_list.append(adata)
    
    # Combine into a single AnnData object if specified
    if combine:
        adata_combined = sc.concat(adata_list, join='outer', label="sample", keys=sample_names)
        return adata_combined
    
    return adata_list



```

8. List the files required to read in and sample names and create you object. Note that in this example, I have randomly downsampled the object to 2000 cells for spead and a feasability check:

```Python

import numpy as np
# List of file paths for 10x .h5 files
file_paths = ["/rds/projects/c/croftap-mapjagb10/scRNAseq/count/S2519_021LK/outs/per_sample_outs/S2519_021LK/count/sample_filtered_feature_bc_matrix.h5", "/rds/projects/c/croftap-mapjagb10/scRNAseq/count/S2519_021RK/outs/per_sample_outs/S2519_021RK/count/sample_filtered_feature_bc_matrix.h5"]

sample_names = ["S1921LK", "S1921RK"]

combined_adata = load_10x_h5_to_anndata(file_paths, sample_names=sample_names, combine=True)

combined_adata.var_names_make_unique()
combined_adata.obs_names_make_unique()

current_num_cells = combined_adata.n_obs
random_indices = np.random.choice(current_num_cells, 2000, replace=False)
combined_adata_small = combined_adata[random_indices].copy()
combined_adata_small


```



9. Run QC and plot

```Python


#QC
combined_adata.var["mt"] = combined_adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    combined_adata, qc_vars=["mt"], inplace=True, log1p=True
)

sc.pl.violin(
    combined_adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

```



