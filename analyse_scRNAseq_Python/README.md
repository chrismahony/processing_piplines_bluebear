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
