# Introduction
Downloading and processing SRA files in parallel with input of metadata file downloaded from SRA run selector. This is by no means a professional project, and it is a personal project I embarked to ease my own sufferings from RNA-seq alignments. Always happy to hear feedbacks and what can be done better

# Documentation
To run the entire pipeline, you can configure the necessary parameters in the config file, followed by runnning the fastq_download_unstable.py (I will have to rename this) script. Otherwise, you can run each module individually.

## Server set up
Not all servers are equal, some are harder to deal with than the rest. I am currently working on SLURM compatible servers, however there might be other form of server managements I am not aware of and would look into if I have the time. Beyond that, do note that as STAR aligner is memory intensive, not all the parallelized operations can be done smoothly. As such, there are parameters in the config file to aid in batching the files such that any server SHOULD be able to run the pipeline without dying. Proxy option is included as well, but it is highly not recommended due to instability of downloads from proxies via cURL

## Config file
YAML formatted config file should point filepaths to directories, and if the directories do not exist, utilities will create a new directory according to the path. Do note that this config file is necessary for all scripts, even if you were to run the modules individually, and you have to pass it as a command line argument

## Bulk RNA-seq pipeline
Currently the only available fully built pipeline. Handles from FASTQ downloads to generating counts matrix. Alignment currently hardcoded and to be done on GENCODE v46 genome files. STAR aligner is used, and the features are not very robust as of yet (memory and nThreads are set as of now, unless there is good reason to make it a configuration I will leave it as such. Feel free to ammend to own needs) featureCounts is currently used to generate the counts matrix from the BAM files.
