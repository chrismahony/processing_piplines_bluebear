```markdown

# Downloading data from Novogene

Once seq is completed, confirm with Novogene to intiate data release.

Follow email insctructions to login to avalible data.

Select batch.

Select all files and click 'batch download'.

Copy this file into the directory where you wan to download the data.

No copy the script below and modify the wget line to the name of your csv


```

```bash

#!/bin/bash
#SBATCH -n 30
#SBATCH -N 1
#SBATCH --mem 199G
#SBATCH --time 49:0:0
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


set -e

module purge; module load bluebear

wget -i ./X204SC24103158-Z01-F001.csv
for f in *.tar; do tar xf "$f"; done

```

```markdown

Next open a terminal and login

Navigate to the dir where the script is saved.

Submit the script, done!
```


