
# Downloading data from Novogene

1. Once seq is completed, confrim with Novogene to intiate data release.

2. Follow email insctructions to login to avalible data.

3. Select batch to download.

4. Select all files and click 'export link'.

5. Copy this file into the directory where you wan to download the data. Probaly create a folder 'fastqs' as this is a logical way to name the data downalded

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

7. Save the script (from above) in the same dir as the csv file, lest assume you call it:

```bash

download_data.txt

```

8. Next open a terminal and login

9. Navigate to the dir where the script is saved. E.g.:

 ```bash

cd /rds/projects/c/croftap-XXXX/my_project/

 ```




11. Submit the script like this

 ```bash

sbatch download_data.txt

 ```

## Troubleshooting

Upon submission you might get this error:

 ```bash

sbatch: error: Batch script contains DOS line breaks (\r\n)
sbatch: error: instead of expected UNIX line breaks (\n).

 ```

Run this and the re-submit the job:

```bash

dos2unix download_data.txt

 ```
