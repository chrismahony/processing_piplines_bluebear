# Processing scATACseq data


1. Complete this step first: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene


2. Then make a directory call something like 'all_fastqs'. Put all the fastq files from all samples that you want to process. Do not make sub directories, just put all files in the one dir 'all_fastqs'

Make a directory called something like 'count', so your project folder would look like:
/croftap-XXXX/
|-- my_project
| |--fastqs #downloaded fatsqs from novogene
| |--all_fastqs   #where you put all fastqs
| |--count   #where you will run the next steps and perform alingment


You may have a load of sub dirs, each with fastq files in them, if so, open a termonal, cd to where all the sub dirs are and run this to move all the fastq files out:


```bash

for dir in */; do
    mv "$dir"*.fastq* .
done

```



4. The copy the script below and save as something like: count.txt
!Important!
Be sure to amend the:

```bash

#SBATCH --account=croftap-XXXX

```

To a RDS folder you have access to.

You will also need to amend:


```bash
FASTQ_DIR #the folder you put all the fastq files
OUTPUT_DIR #where you wan to store the data
REF_PATH #where you ref genome is

```
