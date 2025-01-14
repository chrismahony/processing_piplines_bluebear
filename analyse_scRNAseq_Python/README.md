# Analysing scRNAseq data using Python and Scanpy

Before you start: Optional- Running Cellbender

e.g. here:
https://github.com/chrismahony/processing_piplines_bluebear/blob/main/analyse_scRNAseq_Python/cellbender.txt

And as an array:
https://github.com/chrismahony/processing_piplines_bluebear/blob/main/analyse_scRNAseq_Python/cell_bender_array.txt


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
pip install igraph
pip install louvain
pip intall leidenlg
#pip install git+https://github.com/chrismahony/mcseqy.git #ignore for now

```

# Troubleshooting with louvain clustering!

I commonly go this error messgae:

```Python


---------------------------------------------------------------------------
AttributeError                            Traceback (most recent call last)
Cell In [10], line 4
      2 max_resolution = 2.0
      3 while resolution <= max_resolution:
----> 4     sc.tl.louvain(combined_adata, resolution=resolution, key_added=f'louvain_{resolution:.1f}')
      5     num_clusters = combined_adata.obs[f'louvain_{resolution:.1f}'].nunique()
      6     print(f"Resolution: {resolution:.1f}, Number of clusters: {num_clusters}")

File /rds/projects/c/croftap-mapjag-study/Beth/Python/my-virtual-env-icelake/lib/python3.10/site-packages/legacy_api_wrap/__init__.py:80, in legacy_api.<locals>.wrapper.<locals>.fn_compatible(*args_all, **kw)
     77 @wraps(fn)
     78 def fn_compatible(*args_all: P.args, **kw: P.kwargs) -> R:
     79     if len(args_all) <= n_positional:
---> 80         return fn(*args_all, **kw)
     82     args_pos: P.args
     83     args_pos, args_rest = args_all[:n_positional], args_all[n_positional:]

File /rds/projects/c/croftap-mapjag-study/Beth/Python/my-virtual-env-icelake/lib/python3.10/site-packages/scanpy/tools/_louvain.py:176, in louvain(adata, resolution, random_state, restrict_to, key_added, adjacency, flavor, directed, use_weights, partition_type, partition_kwargs, neighbors_key, obsp, copy)
    174 if use_weights:
    175     partition_kwargs["weights"] = weights
--> 176 if Version(louvain.__version__) < Version("0.7.0"):
    177     louvain.set_rng_seed(random_state)
    178 else:

AttributeError: module 'louvain' has no attribute '__version__'

```


And was solved by running this:


```Python
import louvain
louvain.__version__ = "0.8.2"
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





```

8. List the files required to read in and sample names and create you object. Note that in this example, I have randomly downsampled the object to 2000 cells for spead and a feasability check:

```Python

# List of file paths for 10x .h5 files
file_paths = ["/rds/projects/c/croftap-mapjagb10/scRNAseq/count/S2519_021LK/outs/per_sample_outs/S2519_021LK/count/sample_filtered_feature_bc_matrix.h5", "/rds/projects/c/croftap-mapjagb10/scRNAseq/count/S2519_021RK/outs/per_sample_outs/S2519_021RK/count/sample_filtered_feature_bc_matrix.h5"]

downsample = 1000
sample_names = ["S1921LK", "S1921RK"]

combined_adata = load_10x_h5_to_anndata(file_paths, sample_names, combine=True, downsample=downsample)

combined_adata.var_names_make_unique()
combined_adata.obs_names_make_unique()

#optoinal to downsample further if needed
#current_num_cells = combined_adata.n_obs
#random_indices = np.random.choice(current_num_cells, 2000, replace=False)
#combined_adata_small = combined_adata[random_indices].copy()
#combined_adata_small


```


Read in files from a dir based on a string and generate a list of files names


Generate a list of file names and sample names based on the string: '_cellbender_filtered.h5'

```Python

import glob

directory = Path("/my_path/ambient_RNA_output")
pattern = f"{directory}/*_cellbender_filtered.h5"
file_paths = glob.glob(pattern)
sample_names = [os.path.basename(file).split("_cellbender_filtered.h5")[0] for file in file_paths]

```


Add in further meta.data from csv file

```Python

metadata = pd.read_csv("/my_path/metadata.csv")
adata.obs["sample"] = adata.obs["sample"].astype(str)

for column in metadata.columns:
    if column != "sample":  # Skip the sample column itself
        combined_adata.obs[column] = combined_adata.obs["sample"].map(
            metadata.set_index("sample")[column]
        )

combined_adata
```




An alternaive way to generate a list of files and sample names:

```Pyhton

# Get file names
directory = '/path/to/your/directory'
file_list = [f for f in os.listdir(directory) if f.endswith('.h5') and os.path.isfile(os.path.join(directory, f))]
base_path = Path("/rds/projects/my_path/")
file_paths = [base_path / file for file in file_list]

# Get Sample names
sample_names = [file.split('_')[0] for file in file_list]

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

sc.pp.filter_cells(combined_adata, min_genes=200)
sc.pp.filter_cells(combined_adata, max_genes=7000)
sc.pp.filter_genes(combined_adata, min_cells=3)
combined_adata = combined_adata[combined_adata.obs['pct_counts_mt'] < 10, :]
sc.pp.scrublet(combined_adata, batch_key="sample")
mask = ~combined_adata.obs['predicted_doublet']
combined_adata = combined_adata[mask].copy()
combined_adata

```


11. Process and batch corerect

```Python
combined_adata = preprocess_and_analyze(
    combined_adata, 
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
    sc.tl.louvain(combined_adata, resolution=resolution, key_added=f'louvain_{resolution:.1f}')
    num_clusters = combined_adata.obs[f'louvain_{resolution:.1f}'].nunique()
    print(f"Resolution: {resolution:.1f}, Number of clusters: {num_clusters}")
    
    if num_clusters >= 8:
        print(f"Running sc.tl.rank_genes_groups for resolution {resolution:.1f}")
        sc.tl.rank_genes_groups(combined_adata, groupby=f'louvain_{resolution:.1f}', method='t-test')
        break
    
    
    resolution += 0.1


if resolution > max_resolution:
    print("Did not find exactly 8 clusters within the given resolution range.")

```


13. Visulise clusters and top genes

```Python

sc.pl.umap(
    combined_adata,
    color=["louvain_0.2", "sample"],
    legend_loc="on data",
)


sc.pl.rank_genes_groups_dotplot(
    combined_adata, groupby="louvain_1.0", standard_scale="var", n_genes=10
)

```

14. Assign cell labels and replot dotplot

```Python

combined_adata.obs["cell_types"] = combined_adata.obs["louvain_0.2"].map(
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

sc.tl.rank_genes_groups(combined_adata, groupby="cell_types", method="wilcoxon")


sc.pl.rank_genes_groups_dotplot(
    combined_adata, groupby="cell_types", standard_scale="var", n_genes=10
)

```


15. Save load your anndata object

```Python
#save
combined_adata.write('/my_path/combined_adata.h5ad')

#read
combined_adata = sc.read("/my_path/combined_adata.h5ad")


```


16. For working with larger data sets need to run through sbatch, set up the python script as in:

process.py

That is attahced in this folder


17. The set up the sbatch script to ecxcut this:

run.txt


18. Adding TCR/BCR data once you have a combined adata obj. First create a new dir, in this create a dir for each sample and in each sample dir create a vdj_t and vdj_b dir.

Next you need to copy 'filtered_contig_annotations.csv' and 'clonotypes.csv' from the vdj_t or vdj_b cell ranger outs. Then you should be able to run the script below.


```bash

import os
import pandas as pd

# Path to the parent directory containing sample folders
parent_dir = "/rds/projects/c/croftap-mapjagb10/scRNAseq/count/all_clonotye"

# Initialize DataFrames for TCR and BCR
all_tcr_data = []
all_bcr_data = []

# Iterate over all sample folders
for sample_name in os.listdir(parent_dir):
    sample_folder = os.path.join(parent_dir, sample_name)
    if not os.path.isdir(sample_folder):
        continue  # Skip if not a directory

    # TCR data paths
    tcr_path = os.path.join(sample_folder, "vdj_t/filtered_contig_annotations.csv")
    clono_tcr_path = os.path.join(sample_folder, "vdj_t/clonotypes.csv")

    # BCR data paths
    bcr_path = os.path.join(sample_folder, "vdj_b/filtered_contig_annotations.csv")
    clono_bcr_path = os.path.join(sample_folder, "vdj_b/clonotypes.csv")

    # Process TCR data if files exist
    if os.path.exists(tcr_path) and os.path.exists(clono_tcr_path):
        tcr_data = pd.read_csv(tcr_path)
        clono_tcr_data = pd.read_csv(clono_tcr_path)

        tcr_data = tcr_data[~tcr_data['barcode'].duplicated()][['barcode', 'raw_clonotype_id']]
        tcr_data.rename(columns={'raw_clonotype_id': 'clonotype_id'}, inplace=True)
        tcr_data = tcr_data.merge(clono_tcr_data[['clonotype_id', 'cdr3s_aa']], on='clonotype_id', how='left')
        tcr_data['sample'] = sample_name  # Add sample name for reference
        all_tcr_data.append(tcr_data)

    # Process BCR data if files exist
    if os.path.exists(bcr_path) and os.path.exists(clono_bcr_path):
        bcr_data = pd.read_csv(bcr_path)
        clono_bcr_data = pd.read_csv(clono_bcr_path)

        bcr_data = bcr_data[~bcr_data['barcode'].duplicated()][['barcode', 'raw_clonotype_id']]
        bcr_data.rename(columns={'raw_clonotype_id': 'clonotype_id'}, inplace=True)
        bcr_data = bcr_data.merge(clono_bcr_data[['clonotype_id', 'cdr3s_aa']], on='clonotype_id', how='left')
        bcr_data['sample'] = sample_name  # Add sample name for reference
        all_bcr_data.append(bcr_data)

# Combine all TCR and BCR data into single DataFrames
tcr_combined = pd.concat(all_tcr_data, ignore_index=True) if all_tcr_data else pd.DataFrame()
bcr_combined = pd.concat(all_bcr_data, ignore_index=True) if all_bcr_data else pd.DataFrame()

# Format TCR data for merging with AnnData
if not tcr_combined.empty:
    tcr_combined.set_index('barcode', inplace=True)
    tcr_combined.rename(columns={'clonotype_id': 'T_clonotype_id', 'cdr3s_aa': 'T_cdr3s_aa'}, inplace=True)

# Format BCR data for merging with AnnData
if not bcr_combined.empty:
    bcr_combined.set_index('barcode', inplace=True)
    bcr_combined.rename(columns={'clonotype_id': 'B_clonotype_id', 'cdr3s_aa': 'B_cdr3s_aa'}, inplace=True)


tcr_combined = tcr_combined[tcr_combined.index.isin(combined_adata.obs.index)]
bcr_combined = bcr_combined[bcr_combined.index.isin(combined_adata.obs.index)]

# Resolve duplicates by keeping the first occurrence
tcr_combined = tcr_combined.groupby(tcr_combined.index).first()
bcr_combined = bcr_combined.groupby(bcr_combined.index).first()

# Reindex to match AnnData indices
tcr_combined = tcr_combined.reindex(combined_adata.obs.index)
bcr_combined = bcr_combined.reindex(combined_adata.obs.index)

# Merge into AnnData.obs
if not tcr_combined.empty:
    combined_adata.obs = combined_adata.obs.merge(
        tcr_combined, how='left', left_index=True, right_index=True
    )
if not bcr_combined.empty:
    combined_adata.obs = combined_adata.obs.merge(
        bcr_combined, how='left', left_index=True, right_index=True
    )


combined_adata

```


19. For further downstream analysis, subset a particular cell type:

```python

subset_adata = adata[adata.obs['cell_type'] == 'T_cells'].copy()

```
