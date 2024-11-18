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

import numpy as np
from pathlib import Path

def load_10x_h5_to_anndata(file_paths, sample_names=None, combine=False, downsample=None):
    """
    Load 10x Genomics .h5 files into AnnData objects and optionally combine them.
    Also supports downsampling to a fixed number of cells per sample.

    Parameters:
    - file_paths: List of strings or Path objects, each path pointing to a 10x .h5 file.
    - sample_names: List of sample names, one per file path. Defaults to the file name if not provided.
    - combine: Boolean, if True, combines all matrices into a single AnnData object.
    - downsample: Integer, the number of cells to retain for each sample. If None, no downsampling is done.

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
        
        # Downsample if specified
        if downsample is not None and downsample < adata.n_obs:
            # Randomly select `downsample` number of cells
            adata = adata[np.random.choice(adata.n_obs, downsample, replace=False), :]
        
        # Append AnnData object to list
        adata_list.append(adata)
    
    # Combine into a single AnnData object if specified
    if combine:
        adata_combined = sc.concat(adata_list, join='outer', label="sample", keys=sample_names)
        return adata_combined
    
    return adata_list


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



import numpy as np
# List of file paths for 10x .h5 files
file_paths = ["/rds/projects/c/croftap-mapjagb10/scRNAseq/count/S2519_021LK/outs/per_sample_outs/S2519_021LK/count/sample_filtered_feature_bc_matrix.h5", "/rds/projects/c/croftap-mapjagb10/scRNAseq/count/S2519_021RK/outs/per_sample_outs/S2519_021RK/count/sample_filtered_feature_bc_matrix.h5"]

downsample = 1000
sample_names = ["S1921LK", "S1921RK"]

combined_adata = load_10x_h5_to_anndata(file_paths, sample_names, combine=True, downsample=downsample)

combined_adata.var_names_make_unique()
combined_adata.obs_names_make_unique()


combined_adata.var["mt"] = combined_adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    combined_adata, qc_vars=["mt"], inplace=True, log1p=True
)

sc.pp.filter_cells(combined_adata, min_genes=200)
sc.pp.filter_cells(combined_adata, max_genes=7000)
sc.pp.filter_genes(combined_adata, min_cells=3)
combined_adata = combined_adata[combined_adata.obs['pct_counts_mt'] < 10, :]
sc.pp.scrublet(combined_adata, batch_key="sample")
mask = ~combined_adata.obs['predicted_doublet']
combined_adata = combined_adata[mask].copy()

combined_adata = preprocess_and_analyze(
    combined_adata, 
    n_top_genes=2000, 
    batch_correction=True, 
    batch_key="sample"
)


resolution = 0.1
max_resolution = 2.0
while resolution <= max_resolution:
    sc.tl.louvain(combined_adata, resolution=resolution, key_added=f'louvain_{resolution:.1f}')
    num_clusters = combined_adata.obs[f'louvain_{resolution:.1f}'].nunique()
    print(f"Resolution: {resolution:.1f}, Number of clusters: {num_clusters}")
    
    if num_clusters == 8:
        print(f"Running sc.tl.rank_genes_groups for resolution {resolution:.1f}")
        sc.tl.rank_genes_groups(combined_adata, groupby=f'louvain_{resolution:.1f}', method='t-test')
        break
    
    
    resolution += 0.1


if resolution > max_resolution:
    print("Did not find exactly 8 clusters within the given resolution range.")

combined_adata.write('/my_path/combined_adata.h5ad')
