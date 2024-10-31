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




import scanpy as sc

def preprocess_and_analyze(adata, n_top_genes=2000, batch_correction=False, batch_key="sample"):
    """
    Preprocesses the AnnData object and performs PCA and UMAP analysis.

    Parameters:
    adata (AnnData): The AnnData object to preprocess.
    n_top_genes (int): Number of highly variable genes to select. Default is 2000.
    batch_correction (bool): Whether to apply batch correction using BBKNN. Default is False.
    batch_key (str): The key in obs to use for batch information. Default is "sample".

    Returns:
    AnnData: The processed AnnData object with PCA and UMAP embeddings.
    """
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
     
    if batch_correction:
        sc.external.pp.bbknn(adata, batch_key=batch_key)
    
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=[batch_key])

    return adata





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

10. Filter based on QC metrics and remove dublets:

```Python

sc.pp.filter_cells(combined_adata_small, min_genes=200)
sc.pp.filter_cells(combined_adata_small, max_genes=7000)
sc.pp.filter_genes(combined_adata_small, min_cells=3)
combined_adata_small = combined_adata_small[combined_adata_small.obs['pct_counts_mt'] < 10, :]
sc.pp.scrublet(combined_adata_small, batch_key="sample")
mask = ~combined_adata_small.obs['predicted_doublet']
combined_adata_small = combined_adata_small[mask].copy()
combined_adata_small

```


11. Process and batch corerect

```Python
combined_adata_small = preprocess_and_analyze(
    combined_adata_small, 
    n_top_genes=2000, 
    batch_correction=True, 
    batch_key="sample"
)

```


12. Find n clusters and then find genes once this has been reached

```Python
resolution = 0.1
max_resolution = 2.0
while resolution <= max_resolution:
    sc.tl.louvain(combined_adata_small, resolution=resolution, key_added=f'louvain_{resolution:.1f}')
    num_clusters = combined_adata_small.obs[f'louvain_{resolution:.1f}'].nunique()
    print(f"Resolution: {resolution:.1f}, Number of clusters: {num_clusters}")
    
    if num_clusters == 8:
        print(f"Running sc.tl.rank_genes_groups for resolution {resolution:.1f}")
        sc.tl.rank_genes_groups(combined_adata_small, groupby=f'louvain_{resolution:.1f}', method='t-test')
        break
    
    
    resolution += 0.1


if resolution > max_resolution:
    print("Did not find exactly 8 clusters within the given resolution range.")

```


13. Visulise clusters and top genes

```Python

sc.pl.umap(
    combined_adata_small,
    color=["louvain_0.2", "sample"],
    legend_loc="on data",
)


sc.pl.rank_genes_groups_dotplot(
    combined_adata_small, groupby="louvain_1.0", standard_scale="var", n_genes=10
)

```

14. Assign cell labels and replot dotplot

```Python

combined_adata_small.obs["cell_types"] = combined_adata_small.obs["louvain_0.2"].map(
    {
        "0": "Tcells",
        "1": "fibs",
        "2": "macs",
        "3": "adipose",
        "4": "Bcells",
        "5": "fibs2",
        "6": "contam",
        "7": "contam2",
    }
)

sc.tl.rank_genes_groups(combined_adata_small, groupby="cell_types", method="wilcoxon")


sc.pl.rank_genes_groups_dotplot(
    combined_adata_small, groupby="cell_types", standard_scale="var", n_genes=10
)

```

