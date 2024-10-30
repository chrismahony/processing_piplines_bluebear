# Processing scRNAseq data

This is an approch to process scRNAseq data using cellrnager. 
<br>

The script will find all fastq files in a folder, extract their names an dprocess them with CellRanger (10x Genomics). 
<br>


1. Complete this step: https://github.com/chrismahony/processing_piplines_bluebear/tree/main/dowloand_data_novogene
 <br>

2. Then make a directory call something like 'all_fastqs'. Put all the fastq files from all samples that you want to process.
 <br>
   Do not make sub directories, just put all files in the one dir 'all_fastqs'
<br>
3. Make a directory called something like 'count', so your project folder would look like:
