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

pip install Scanpy

```
