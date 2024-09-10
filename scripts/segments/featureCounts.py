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

class featurecounts_generate_counts_matrix:
    def __init__(self, params, logger, dependency=None):
        self.params = params
        self.logger = logger 
        self.dependency = dependency


    def execute(self):
        genome_gtf = [os.path.join(self.params['genome_dir'], file) for file in os.listdir(self.params['genome_dir']) if file.endswith('.gtf')][0]
        # paired_files = ' '.join([os.path.join(self.params['bam_dir'], i) for i in os.listdir(self.params['bam_dir']) if 'pairedAligned.sortedByCoord.out.bam' in i])
        # single_files = ' '.join([os.path.join(self.params['bam_dir'], i) for i in os.listdir(self.params['bam_dir']) if 'singleAligned.sortedByCoord.out.bam' in i])
        counts_matrix_fp = os.path.join(self.params['counts_dir'], 'counts')
        
        # paired_cmd = f'featureCounts -T 6 -a {genome_gtf} -o {counts_matrix_fp}_paired_reads.txt -p -B -C {paired_files}'
        # single_cmd = f'featureCounts -T 6 -a {genome_gtf} -o {counts_matrix_fp}_single_reads.txt {single_files}'
        cmd = f"""
# Split the paired and single BAM files into batches of 10
paired_files=$(ls {os.path.join(self.params['bam_dir'], '*pairedAligned.sortedByCoord.out.bam')} 2>/dev/null)
single_files=$(ls {os.path.join(self.params['bam_dir'], '*singleAligned.sortedByCoord.out.bam')} 2>/dev/null)

echo "$paired_files" | split -l 1 - bam_batch_paired_

echo "$single_files" | split -l 1 - bam_batch_single_

# Process each batch of paired BAM files
for batch in bam_batch_paired_*; do
    batch_id=$(basename $batch)
    featureCounts -T 6 -a {genome_gtf} -o {counts_matrix_fp}_paired_reads_${{batch_id}}.txt -p -B -C $(cat $batch)
done

# Process each batch of single BAM files
for batch in bam_batch_single_*; do
    batch_id=$(basename $batch)
    featureCounts -T 6 -a {genome_gtf} -o {counts_matrix_fp}_single_reads_${{batch_id}}.txt $(cat $batch)
done

# After processing, concatenate all paired and single CSV files into a single CSV each
cat {counts_matrix_fp}_paired_reads_*.txt > {counts_matrix_fp}_final_paired_reads.txt
cat {counts_matrix_fp}_single_reads_*.txt > {counts_matrix_fp}_final_single_reads.txt

# Optionally, clean up the batch files and intermediate results
rm bam_batch_paired_*
rm bam_batch_single_*
rm {counts_matrix_fp}_paired_reads_*.txt
rm {counts_matrix_fp}_single_reads_*.txt
"""
        if self.params['slurm']:
            job_id = []
            
            # single_slurm_script_path = generate_slurm_script(f'{slurm_script}\n{single_cmd}', "counts_matrix_single", 6,  self.params['slurm_dir'], "16G", self.dependency)
            # paired_slurm_script_path = generate_slurm_script(f'{slurm_script}\n{paired_cmd}', "counts_matrix_paired", 6,  self.params['slurm_dir'], "16G", self.dependency)
            slurm_script_path = generate_slurm_script(cmd, 'featureCounts', 6, self.params['slurm_dir'], "64G", self.dependency)
            # job_id.append(submit_slurm_job(single_slurm_script_path, self.logger, self.dependency))
            # job_id.append(submit_slurm_job(paired_slurm_script_path, self.logger, self.dependency))
            job_id = submit_slurm_job(slurm_script_path, self.logger, self.dependency)
            # for idx in job_id:
            #     self.logger.info(f'{idx} submitted.')
            return job_id
        else:
            result = subprocess.run(cmd, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)
            self.logger.info(f"Command executed successfully: {result.stdout}")
            self.logger.info('Counts matrix generation completed.')
            return None

    def run_cmd(self, cmd):
        '''
            Basic
        '''
        self.logger.info(f'Running {cmd}')
        subprocess.run(cmd, shell=True, check=True)
        self.logger.info('Counts matrix generation completed for: {cmd}')

if __name__ == '__main__':
    args = create_parser().parse_args()
    config_file = args.config_file

    with open(config_file, 'r') as f:
        params = yaml.safe_load(f)

    logger = setup_logger(params['log_file'])
    create_directories(params, logger)
    logger.info('This script is only running Counts Generations segment.')
    featurecounts_generate_counts_matrix(params, logger).execute()
