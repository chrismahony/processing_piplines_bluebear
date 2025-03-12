# Uploading your sequenced data to GEO

1. First login to: https://www.ncbi.nlm.nih.gov/geo/info/submission.html

2. Now download metadata excel sheet

3. Next gather all your fastqs for a specific experiment in a folder.

   Important- make seperate submisison for different experiments. If for a single paper you have done one bulk RNAseq experiment and two scRNAseq experiment, then this is three seperate submissions.

4. Fill out the spead sheet. The easisest way to gather all you file names is this:

```bash

cd /your_dir_with_all_fastqs/

ls > files.txt

```

Now open 'files.txt' and copy and past all the file names into excel


5. Next follow the instruction on the GEO website to login and naviagte to the upload folder.

For me, using BlueBear, the sftp option works best


6. With the sftp login, create a folder that exactly matches the folder name where all the fastqs are stored

7. Now transfer the files!

```bash
#assume that you have created a folder in the geo upload repo called 'all_fastqs'
put -r /rds/projects/my_path/all_fastqs

```

8. Once this is done, submit the metadata via the GEO upload system, they will contact you about missing or corroupt files.
