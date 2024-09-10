import subprocess
from tqdm import tqdm 
import logging
import math
import re
import os
import sys
from argparse import ArgumentParser
from pathlib import Path
import shutil


def create_parser():
    parser = ArgumentParser(
        prog="sra_scraper_prototype",
        description="Test SRA Scraper prototype"
    )
    parser.add_argument("-c", "--config_file", type=Path, help="Path to YAML config file containing all variables.")
    return parser


def setup_logger(log_file_path):
    logging.basicConfig(
        filename=log_file_path,
        encoding='utf-8',
        level=logging.DEBUG,
        format='%(asctime)s %(threadName)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        filemode='w'
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s - %(message)s'))

    logger = logging.getLogger()
    logger.addHandler(console_handler)
    
    logger.info('Script started')

    return logger


def gunzip_files(file_path):
    if file_path.endswith('.gz'):
        subprocess.run(['gunzip', file_path], check=True)
        return file_path[:-3]  # Decompressed file path
    return file_path

def clear_directory_contents(directory, logger):
    # Check if the directory exists and is not empty
    if os.path.exists(directory):
        # Iterate over each item in the directory
        for item in os.listdir(directory):
            item_path = os.path.join(directory, item)
            try:
                # Remove directories and their contents
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                # Remove individual files
                elif os.path.isfile(item_path) or os.path.islink(item_path):
                    os.remove(item_path)
                logger.info(f"Removed {item_path}")
            except Exception as e:
                logger.error(f"Failed to remove {item_path}. Reason: {e}")
    logger.info(f'Emptied {directory}')

def create_directories(params, logger):
    clear_directory_contents(params['tmp_dir'], logger)
    clear_directory_contents(params['slurm_dir'], logger)
   
    directories = [dir_fp for key, dir_fp in params.items() if 'dir' in key]
    
    for directory in directories:
        path = Path(directory)
        if not path.exists():
            logger.info(f"Creating directory: {directory}")
            path.mkdir(parents=True, exist_ok=True)
        else:
            logger.info(f"Directory {directory} already exists.")


def obtaining_ftp_fps(sra_id, logger, proxy = None):
    try:
        # Run the subprocess to get the FTP file paths from the ENA API
        if proxy:
            cmd = ['curl', '-x', proxy,'-X', 'GET', f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sra_id}&fields=fastq_ftp&result=read_run']
        else:
            cmd = ['curl', '-X', 'GET', f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sra_id}&fields=fastq_ftp&result=read_run']
        response = subprocess.run(
            cmd,
            capture_output=True,
            check=True  # Ensure subprocess raises an error for non-zero exit codes
        )
        output = response.stdout.decode()

        # Check if the response contains an invalid SRA ID message
        if 'not valid' in output:
            logger.warning(f'{sra_id} is not valid!')
            return []

        # Split the output to get FTP file paths
        ftp_fps = output.split('\n')[1].split('\t')[0]
        return ftp_fps.split(';') if ';' in ftp_fps else [ftp_fps]

    except subprocess.CalledProcessError as e:
        # Log an error if the subprocess fails (e.g., network error)
        logger.error(f"Failed to retrieve FTP paths for {sra_id}: {e.stderr.decode()}")
        return []

    except IndexError:
        # Log an error if the expected data format is not found
        logger.error(f"Unexpected response format when retrieving FTP paths for {sra_id}: {output}")
        return []

    except Exception as e:
        # Catch any other unexpected exceptions and log them
        logger.error(f"An unexpected error occurred for {sra_id}: {str(e)}")
        return []



def string_to_bytes(size_str):
    size_units = {
        'B': 1,
        'KB': 1024,
        'MB': 1024 ** 2,
        'GB': 1024 ** 3,
        'TB': 1024 ** 4,
        'PB': 1024 ** 5
    }

    size = float(re.findall(r'\\d+\\.?\\d*', size_str)[0])
    unit = re.findall(r'[a-zA-Z]+', size_str)[0].upper()

    if unit not in size_units:
        raise ValueError(f"Unknown unit: {unit}")

    return int(size * size_units[unit])


def bytes_to_string(size_bytes):
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    return f"{round(size_bytes / p, 2)} {size_name[i]}"


def obtain_download_size(ftp_fp, logger, proxy = None):
    try:
        # Simulate obtaining the download size from the response header
        if proxy:
            response = subprocess.run(['curl', '-x', proxy,'-I', ftp_fp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            response = subprocess.run(['curl', '-I', ftp_fp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        file_size = bytes_to_string(int(re.search(r'Content-Length: (\d+)', response.stdout.decode()).group(1)))
        logger.info(f"File size for {ftp_fp}: {file_size}")
        return file_size

    except AttributeError:
        # Handle the case where the Content-Length is missing or malformed
        logger.error(f"Failed to obtain file size for {ftp_fp}. Continuing with download...")
        return None


def download_ftp(ftp_fp, out_dir, logger, proxy = None):
    logger.info(f'Downloading {ftp_fp} into {out_dir}')
    obtain_download_size(ftp_fp, logger)

    file_name = ftp_fp.split("/")[-1]
    dest_file = os.path.join(out_dir, file_name)

    progress_regex = re.compile(r'(?P<percent>\d+(\.\d+)?%?)')
    
    if proxy:
        cmd = ['curl', '--create-dirs', '-x', proxy,'--output', dest_file, ftp_fp]
    else:
        cmd = ['curl', '--create-dirs', '--output', dest_file, ftp_fp]
    
    with subprocess.Popen(cmd,
                          stderr=subprocess.PIPE, text=True) as proc:
        with tqdm(total=100, desc=f'Downloading {file_name}', unit='%') as pbar:
            for line in proc.stderr:
                if match := progress_regex.search(line):
                    percent = float(match.group('percent').strip('%'))
                    pbar.n = percent
                    pbar.refresh()

    logger.info(f'Download complete for {file_name}')
    return dest_file


def generate_slurm_script(command, job_name, num_threads, output_dir, memory = '16G', dependency=None):
    dependency_str = f"#SBATCH --dependency=afterany:{dependency}\n" if dependency else ""
    
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={output_dir}/{job_name}_log.out
#SBATCH --error={output_dir}/{job_name}_err.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={num_threads}
#SBATCH --mem={memory}
{dependency_str}
{command}
"""
    script_path = os.path.join(output_dir, f"{job_name}.slurm")
    with open(script_path, 'w') as f:
        f.write(slurm_script)
    
    return script_path


def submit_slurm_job(script_path, logger, dependency=None):
    try:
        logger.info(f"Submitting SLURM job with script: {script_path}")
        slurm_cmd = f"sbatch {script_path}"
        if dependency:
            slurm_cmd += f" --dependency=afterok:{dependency}"
        result = subprocess.run(slurm_cmd, shell=True, check=True, capture_output=True)
        logger.info(f"Job submitted successfully: {result.stdout.decode()}")
        return result.stdout.decode().split()[-1]  # Returns the SLURM job ID
    except subprocess.CalledProcessError as e:
        logger.error(f"Error submitting SLURM job: {e.stderr.decode()}")
        return None