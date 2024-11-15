# Plotting using default segmentation from Xenium analyser



1.  Load session, modules and xenium data as spatial data obj

```Python

import os
from pathlib import Path
import sys
node_type = os.getenv('BB_CPU')
venv_dir = f'/rds/projects/g/gilberts-spatial-biology-image-analysis/Chris/my-virtual-env-icelake'  # edit this line to match the venv directory format
venv_site_pkgs = Path(venv_dir) / 'lib' / f'python{sys.version_info.major}.{sys.version_info.minor}' / 'site-packages'
if venv_site_pkgs.exists():
    sys.path.insert(0, str(venv_site_pkgs))
else:
    print(f"Path '{venv_site_pkgs}' not found. Check that it exists and/or that it exists for node-type '{node_type}'.")

import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq



sdata = xenium("/rds/projects/c/croftap-mapjagx1/data/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2518-11__20231107__215839")
sdata

````



2. Plotting morphoogy plots

```Python

from spatialdata import bounding_box_query

fig, ax = plt.subplots(figsize=(10, 10))

# Render the image and display it
sdata.pl.render_images("morphology_focus").pl.show(ax=ax, title="Morphology Image")

# Set the crop boundaries (adjust the values to your desired crop area)
ax.set_xlim(left=15_000, right=25_000)  # Example x-axis limits for cropping
ax.set_ylim(bottom=0, top=10_000)  # Example y-axis limits for cropping

# Show the cropped plot
plt.show()


```


```Python


from spatialdata import bounding_box_query

fig, ax = plt.subplots(figsize=(10, 10))

# Render the image and display it
#sdata.pl.render_images("cell_boundaries").pl.show(ax=ax, title="Morphology Image")

sdata.pl.render_labels("cell_labels").pl.show(ax=ax, title="Cell labels", coordinate_systems="global")

# Set the crop boundaries (adjust the values to your desired crop area)
ax.set_xlim(left=15_000, right=25_000)  # Example x-axis limits for cropping
ax.set_ylim(bottom=0, top=10_000)  # Example y-axis limits for cropping

# Show the cropped plot
plt.show()

```



3. Combining morphology and cells and adjusting colors


```Python


from spatialdata import bounding_box_query

fig, ax = plt.subplots(figsize=(10, 10))

# Render the image and display it
sdata.pl.render_images("morphology_focus").pl.show(ax=ax, title="Morphology Image")
sdata.pl.render_shapes("cell_boundaries", fill_alpha = 0.2, color="red").pl.show(ax=ax, title="Cell labels", coordinate_systems="global")

# Set the crop boundaries (adjust the values to your desired crop area)
ax.set_xlim(left=18_000, right=20_000)  # Example x-axis limits for cropping
ax.set_ylim(bottom=0, top=4000)  # Example y-axis limits for cropping

# Show the cropped plot
plt.show()


```



Another option

```Python

from spatialdata import bounding_box_query

fig, ax = plt.subplots(figsize=(10, 10))

# Render the image and display it
sdata.pl.render_images("morphology_focus").pl.show(ax=ax, title="Morphology Image")
#sdata.pl.render_shapes("cell_boundaries", fill_alpha = 0.2, color="red").pl.show(ax=ax, title="Cell labels", coordinate_systems="global")
sdata.pl.render_labels("cell_labels", fill_alpha = 0.2).pl.show(ax=ax, title="Cell labels", coordinate_systems="global")

# Set the crop boundaries (adjust the values to your desired crop area)
ax.set_xlim(left=18_000, right=20_000)  # Example x-axis limits for cropping
ax.set_ylim(bottom=0, top=4000)  # Example y-axis limits for cropping

# Show the cropped plot
plt.show()

```




4. Adding nuc boundaries


```Python

from spatialdata import bounding_box_query

fig, ax = plt.subplots(figsize=(10, 10))

# Render the image and display it
sdata.pl.render_images("morphology_focus").pl.show(ax=ax, title="Morphology Image")
sdata.pl.render_shapes("cell_boundaries", fill_alpha = 0.2, color="red").pl.show(ax=ax, title="Cell labels", coordinate_systems="global")
sdata.pl.render_shapes("nucleus_boundaries", color="blue").pl.show(ax=ax, title="nucleus_boundaries", coordinate_systems="global")

 

# Set the crop boundaries (adjust the values to your desired crop area)
ax.set_xlim(left=18_000, right=20_000)  # Example x-axis limits for cropping
ax.set_ylim(bottom=0, top=4000)  # Example y-axis limits for cropping

# Show the cropped plot
plt.show()



```

5. And now addig in txs from a specific gene

```Python

from spatialdata import bounding_box_query

fig, ax = plt.subplots(figsize=(10, 10))

# Render the image and display it
sdata.pl.render_images("morphology_focus").pl.show(ax=ax, title="Morphology Image")
sdata.pl.render_shapes("cell_boundaries", fill_alpha = 0.2, color="red").pl.show(ax=ax, title="Cell labels", coordinate_systems="global")
sdata.pl.render_shapes("nucleus_boundaries", color="blue").pl.show(ax=ax, title="nucleus_boundaries", coordinate_systems="global")

sdata.pl.render_points(
    "transcripts",
    color="feature_name",
    groups='VWF',
    palette="orange",
).pl.show(ax=ax, title="Morphology Image")


# Set the crop boundaries (adjust the values to your desired crop area)
ax.set_xlim(left=18_000, right=20_000)  # Example x-axis limits for cropping
ax.set_ylim(bottom=0, top=4000)  # Example y-axis limits for cropping

# Show the cropped plot
plt.show()

```
