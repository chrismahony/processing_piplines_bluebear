
# Downloading data from Novogene

1. Once seq is completed, confrim with Novogene to intiate data release.

2. Follow email insctructions to login to avalible data.

3. Select batch to download.

4. Select all files and click 'batch download'.

5. Copy this file into the directory where you wan to download the data.

6. Now copy the script below and modify the wget line to the name of your csv



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

wget -i ./X204SC24101141-Z01-F001.csv
for f in *.tar; do tar xf "$f"; done

```

7. Sav ethe script in the same dir as the csv file

8. Next open a terminal and login

9. Navigate to the dir where the script is saved.

10. Submit the script, have a coffee, done!



