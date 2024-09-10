import subprocess
import pandas as pd
from tqdm import tqdm 
import logging
import re
import os
import yaml
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    sys.path.append(parent_dir)
    
try:
    from utils.utilities import *
except ModuleNotFoundError:
    sys.path.append(parent_dir)
    from utilities import *


# Generate genome index using STAR
class generate_genome_index:
    def __init__(self, params, logger, dependency=None):
        self.params = params
        self.logger = logger

    def execute(self):
        num_threads = self.params['num_threads']
        genome_dir = self.params['genome_dir']
        mem_limit = self.params['genome_indexing']['mem_limit']
        

        if self.params['proxy']['exists'] == True:
            proxy = params['proxy']['proxy_server']
        else:
            proxy = None


        # Check for existing genome and annotation files
        gtf_re = re.compile(r'.*annotation.*\.gtf.*')
        fa_re = re.compile(r'.*genome.*\.fa.*')
        gz_re = re.compile(r'.*\.gz$')
        files = os.listdir(genome_dir)

        gtf_file = any(gtf_re.match(file) for file in files)
        fa_file = any(fa_re.match(file) for file in files)

        # FTP links for downloading genome files if not present
        gtf_fp = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz'
        fasta_fp = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz'
        
        if not gtf_file or not fa_file:
            self.logger.info('Downloading genome/annotation files')
            if not gtf_file:
                download_ftp(gtf_fp, genome_dir, self.logger, proxy)
            if not fa_file:
                download_ftp(fasta_fp, genome_dir, self.logger, proxy)
        else:
            self.logger.info('Genome and annotation files already present')

        # Unzip .gz files if necessary
        gz_files = [os.path.join(genome_dir, file) for file in files if gz_re.match(file)]
        for gz_file in gz_files:
            self.logger.info(f'Unzipping: {gz_file}')
            gunzip_files(gz_file)

        # Identify unzipped files
        files_new = os.listdir(genome_dir)
        gtf_file_new = [os.path.join(genome_dir, file) for file in files_new if gtf_re.match(file)]
        fa_file_new = [os.path.join(genome_dir, file) for file in files_new if fa_re.match(file)]
        
        if len(gtf_file_new) != 1 or len(fa_file_new) != 1:
            self.logger.error('Error in genome directory, multiple or no files found')
            exit(1)
        
        genome_gtf = gtf_file_new[0]
        genome_fasta = fa_file_new[0]
        
        # STAR command for genome index generation
        self.logger.info('Generating genome index')
        mem_to_use = string_to_bytes(mem_limit)
        star_cmd = f"STAR --runThreadN {num_threads} --runMode genomeGenerate --genomeDir {genome_dir} --genomeFastaFiles {genome_fasta} \
    --sjdbGTFfile {genome_gtf} --sjdbOverhang 100 --limitGenomeGenerateRAM {mem_to_use}"

        if params['slurm']:
            slurm_script_path = generate_slurm_script(star_cmd, "genome_index", num_threads, mem_limit, self.params['slurm_dir'], dependency)
            job_id = submit_slurm_job(slurm_script_path, self.logger)
            return job_id
        else:
            subprocess.run(star_cmd, shell=True, check=True, executable='/bin/bash')
            self.logger.info("Genome index generation completed.")
            return None
    
if __name__ == '__main__':
    args = create_parser().parse_args()
    config_file = args.config_file

    with open(config_file, 'r') as f:
        params = yaml.safe_load(f)

    logger = setup_logger(params['log_file'])
    create_directories(params, logger)
    logger.info('This script is only running Genome Index segment.')
    generate_genome_index(params, logger).execute()
