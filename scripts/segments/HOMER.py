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


def run_cmd(cmd, logger):
    '''
        Basic
    '''
    logger.info(f'Running {cmd}')
    subprocess.run(cmd, shell=True, check=True)
    logger.info('Counts matrix generation completed for: {cmd}')

def generate_counts_matrix(params, logger, dependency=None):
    genome_gtf = [os.path.join(params['genome_dir'], file) for file in os.listdir(params['genome_dir']) if file.endswith('.gtf')][0]
    files = ' '.join([os.path.join(params['bam_dir'], i) for i in os.listdir(params['bam_dir']) if 'Aligned.sortedByCoord.out.bam' in i])
    counts_matrix_fp = os.path.join(params['counts_dir'], 'counts.txt')

    if params['paired_reads']:
        cmd = f'featureCounts -T 6 -a {genome_gtf} -o {counts_matrix_fp} -p -B -C {files}'
    else:
        cmd = f'featureCounts -T 6 -a {genome_gtf} -o {counts_matrix_fp} {files}'
    
    if params['slurm']:
        slurm_script_path = generate_slurm_script(cmd, "counts_matrix", 6, "16G", params['slurm_dir'], dependency)
        job_id = submit_slurm_job(slurm_script_path, logger, dependency)
        return job_id
    else:
        subprocess.run(cmd, shell=True, check=True)
        logger.info('Counts matrix generation completed.')
        return None


# Generate counts matrix
def generate_counts_matrix_homer(params, logger, dependency=None):
    # Find GTF file in genome directory
    genome_gtf = [os.path.join(params['genome_dir'], file) for file in os.listdir(params['genome_dir']) if file.endswith('.gtf')]
    if not genome_gtf:
        logger.error("GTF file not found in genome directory.")
        return None
    genome_gtf = genome_gtf[0]  # Assume there's only one GTF file

        # Find BAM files
    bam_files = [os.path.join(params['bam_dir'], i) for i in os.listdir(params['bam_dir']) if 'Aligned.sortedByCoord.out.bam' in i]
    if not bam_files:
        logger.error("No BAM files found for generating counts.")
        return None
    bam_files_str = ' '.join(bam_files)
    
    counts_matrix_fp = os.path.join(params['counts_dir'], 'counts')

    # HOMER analyzeRepeats.pl command for raw counts, FPKM, TPM, and RPKM
    cmd = (
        f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {bam_files_str} -raw > {counts_matrix_fp}_raw.txt && '
        f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {bam_files_str} -fpkm > {counts_matrix_fp}_fpkm.txt && '
        f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {bam_files_str} -tpm > {counts_matrix_fp}_tpm.txt && '
        f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {bam_files_str} -rpkm > {counts_matrix_fp}_rpkm.txt'
    )
    
    if params['slurm']:
        # Generate SLURM script for HOMER counts matrix generation
        slurm_script_path = generate_slurm_script(cmd, "counts_matrix", 6, "16G", params['slurm_dir'], dependency)
        job_id = submit_slurm_job(slurm_script_path, logger, dependency)
        return job_id
    else:
        # Run the HOMER command locally
        subprocess.run(cmd, shell=True, check=True)
        logger.info('Counts matrix generation with HOMER completed.')
        return None


    # # Find BAM files (paired and single reads)
    # paired_bam_files = [os.path.join(params['bam_dir'], i) for i in os.listdir(params['bam_dir']) if 'Aligned.sortedByCoord.out.bam' in i and 'paired' in i]
    # single_bam_files = [os.path.join(params['bam_dir'], i) for i in os.listdir(params['bam_dir']) if 'Aligned.sortedByCoord.out.bam' in i and 'single' in i]


    # counts_matrix_fp = os.path.join(params['counts_dir'], 'counts')

    # # Function to run a command (used for parallel execution)


    # # List of commands to run
    # cmds = []

    # # Paired-end reads
    # if paired_bam_files:
    #     paired_bam_files_str = ' '.join(paired_bam_files)
    #     paired_cmd = (
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {paired_bam_files_str} -paired -raw > {counts_matrix_fp}_paired_reads_raw.txt && '
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {paired_bam_files_str} -paired -fpkm > {counts_matrix_fp}_paired_reads_fpkm.txt && '
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {paired_bam_files_str} -paired -tpm > {counts_matrix_fp}_paired_reads_tpm.txt && '
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {paired_bam_files_str} -paired -rpkm > {counts_matrix_fp}_paired_reads_rpkm.txt'
    #     )
    #     cmds.append(paired_cmd)
    
    # # Single-end reads
    # if single_bam_files:
    #     single_bam_files_str = ' '.join(single_bam_files)
    #     single_cmd = (
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {single_bam_files_str} -raw > {counts_matrix_fp}_single_reads_raw.txt && '
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {single_bam_files_str} -fpkm > {counts_matrix_fp}_single_reads_fpkm.txt && '
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {single_bam_files_str} -tpm > {counts_matrix_fp}_single_reads_tpm.txt && '
    #         f'analyzeRepeats.pl {genome_gtf} hg38 -count exons -d {single_bam_files_str} -rpkm > {counts_matrix_fp}_single_reads_rpkm.txt'
    #     )
    #     cmds.append(single_cmd)

    # if params['slurm']:
    #     # Generate SLURM script for HOMER counts matrix generation
    #     for cmd in cmds:
    #         slurm_script_path = generate_slurm_script(cmd, "counts_matrix", 6, "16G", params['slurm_dir'], dependency)
    #         submit_slurm_job(slurm_script_path, logger, dependency)
    # else:
    #     # Run commands in parallel using ProcessPoolExecutor
    #     with ProcessPoolExecutor() as executor:
    #         executor.map(run_cmd, cmds)

    # logger.info('Counts matrix generation with HOMER completed for all files.')
    # return None