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


class generate_bam_files:
    def __init__(self, params, logger):
        self.params = params
        self.logger = logger
        self.sra_ids = pd.read_csv(params['pheno_file'], sep=',').iloc[:, 0].tolist()
        self.fastq_files = os.listdir(params['fastq_dir'])
        self.genome_dir = params['genome_dir']
        self.batch_size = params['star_alignment']['batch_size'] # Retrieve batch size from params, default to 0 if not specified
        self.tmp_dir = params['tmp_dir']
        self.rerun_sra_id = params['star_alignment'].get('rerun_sra_id', None)
        self.sra_rerun_unzip_method = params['star_alignment'].get('sra_rerun_unzip_method', 'gunzip -c')
        
    def execute(self):
        if self.rerun_sra_id:
            # If rerun_sra_id is specified, process only that file
            self.logger.info(f"Rerunning alignment for SRA ID: {self.rerun_sra_id}")
            args = [(self.rerun_sra_id, self.params, self.fastq_files, self.genome_dir, os.path.join(self.tmp_dir, self.rerun_sra_id), self.logger)]
            return self.execute_non_batched(args)
        else:
            args = [(sra_id, self.params, self.fastq_files, self.genome_dir, os.path.join(self.tmp_dir, sra_id),self.logger) for sra_id in self.sra_ids]

            # Proceed with batching if batch_size > 0
            if self.batch_size > 0:
                self.logger.info(f'Batched approach: Aligning in batches of {self.batch_size}.')
                return self.execute_batched(args)
            else:
                self.logger.info('Non-batched approach: Aligning all BAM files at once.')
                return self.execute_non_batched(args)

    def execute_batched(self, args):
        job_ids = []
        batch_dependency = None  # Track dependency between batches

        # Process jobs in batches
        for i in range(0, len(args), self.batch_size):
            batch_args = args[i:i + self.batch_size]  # Select batch of jobs
            batch_job_ids = []

            with ProcessPoolExecutor() as executor:
                futures = [executor.submit(self.process_sra_id, arg, batch_dependency) for arg in batch_args]
                for future in as_completed(futures):
                    job_id = future.result()
                    if job_id:
                        batch_job_ids.append(job_id)

            if self.params['slurm'] and len(batch_job_ids) > 0:
                batch_dependency = ":".join(batch_job_ids)  

            job_ids.extend(batch_job_ids)
            self.logger.info(f'Batch {i//self.batch_size + 1} submitted. Dependency: {batch_dependency}')

        return batch_job_ids if batch_job_ids else job_ids

    def execute_non_batched(self, args):
        job_ids = []

        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(self.process_sra_id, arg) for arg in args]
            for future in as_completed(futures):
                job_id = future.result()
                if job_id:
                    job_ids.append(job_id)

        return job_ids


    def process_sra_id(self, args, dependency=None):
        sra_id, params, fastq_files, genome_dir, tmp_dir, logger = args
        matching_files = [os.path.join(params['fastq_dir'], f) for f in fastq_files if sra_id in f]
        if 0 < len(matching_files) < 3:
            read_type = 'single' if len(matching_files) == 1 else 'paired'
            filepaths_string = ' '.join(matching_files)
            out_path = os.path.join(params['bam_dir'], f'{sra_id}_')
            job_id = self.star_alignment(filepaths_string, sra_id, genome_dir, out_path, tmp_dir, logger, params['slurm'], params['slurm_dir'], dependency=dependency, read_type=read_type, unzip_method = self.sra_rerun_unzip_method)
            return job_id
        return None

    def star_alignment(self, filepaths_string, sra_id, genome_dir, out_path, tmp_dir,logger, slurm=False, slurm_dir=None, dependency=None, read_type='paired', unzip_method = 'gunzip -c'):
        star_cmd = (
            f'STAR --runThreadN 4 --genomeDir {genome_dir} '
            f'--readFilesIn {filepaths_string} --outSAMtype BAM SortedByCoordinate '
            f'--quantMode GeneCounts --readFilesCommand {unzip_method} --outFileNamePrefix {out_path}{read_type} --outTmpDir {tmp_dir}'
        )

        if slurm:
            slurm_script_path = generate_slurm_script(star_cmd, f"star_alignment_{sra_id}", 4, slurm_dir, memory="64G", dependency=dependency)
            job_id = submit_slurm_job(slurm_script_path, logger)
            logger.info(f'{sra_id} submitted with job ID {job_id}')
            return job_id
        else:
            subprocess.run(star_cmd, shell=True, check=True, executable='/bin/bash')
            logger.info(f'Alignment for {filepaths_string} done.')
            return None



if __name__ == '__main__':

    def create_script_parser():
        parser = ArgumentParser(description="STAR Alignment Script")
        parser.add_argument("-c", "--config_file", type=str, help="Path to the YAML config file", required=True)
        parser.add_argument('-i',"--rerun_sra_id", type=str, help="Specific SRA ID to rerun (optional)", default=None)
        parser.add_argument('-m',"--sra_rerun_unzip_method", type=str, help="Method to use in unzipping gz files", choices=['zcat', 'gunzip -c'], default='gunzip -c')
        
        return parser

    args = create_script_parser().parse_args()
    config_file = args.config_file

    with open(config_file, 'r') as f:
        params = yaml.safe_load(f)
    if args.rerun_sra_id:
        params['star_alignment']['rerun_sra_id'] = args.rerun_sra_id
    if args.sra_rerun_unzip_method:
        params['star_alignment']['sra_rerun_unzip_method'] = args.sra_rerun_unzip_method

    logger = setup_logger(params['log_file'])
    create_directories(params, logger)
    logger.info('This script is only running Download FASTQ segment.')
    generate_bam_files(params, logger).execute()
