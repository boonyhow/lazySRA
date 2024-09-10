import pandas as pd
import re
import os
import yaml
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    sys.path.append(parent_dir)
    
try:
    from utils.utilities import *
except ModuleNotFoundError:
    sys.path.append(parent_dir)
    from utilities import *


# Download FASTQ files
class fastqc:
    def __init__(self, params, logger):
        self.params = params
        self.logger = logger

    def fastqc_cmd(self, num_threads, files, outdir, logger):
        # Command to run FastQC on all the files
        cmd = ['fastqc', '-t', str(num_threads), '--outdir', outdir] + files
        logger.info(f'Running command: {" ".join(cmd)}')
        try:
            subprocess.run(cmd, check=True)
            logger.info(f'FASTQC completed for files: {", ".join(files)}')
        except subprocess.CalledProcessError as e:
            logger.error(f'FASTQC failed for files: {", ".join(files)} with error: {e}')

    def execute(self):    
        df = pd.read_csv(self.params['pheno_file'], sep=',')
        sra_ids = df.iloc[:, 0].tolist()  # Assuming the first column contains SRA IDs
        fastq_dir = self.params['fastq_dir']
        out_dir = self.params['fastqc_dir']
        files = os.listdir(fastq_dir)
        
        batch_size = len(files)
        if batch_size > 32:
            batch_size = 32

        # Split into batches based on batch_size
        if len(files) > batch_size:
            ftp_batches = [files[i:i + batch_size] for i in range(0, len(files), batch_size)]
        else:
            ftp_batches = [files]

        for batch_num, batch in enumerate(ftp_batches, start=1):
            if len(batch) > 0:
                self.logger.info(f'Starting FASTQC batch {batch_num}/{len(ftp_batches)} with {len(batch)} files.')

                futures = []
                with ThreadPoolExecutor(max_workers=len(batch)) as executor:
                    # Submit each batch for execution
                    futures.append(
                        executor.submit(self.fastqc_cmd, len(batch), [os.path.join(fastq_dir, file_) for file_ in batch], out_dir, self.logger)
                    )

                # Ensure all tasks are completed
                for future in as_completed(futures):
                    try:
                        future.result()  # This will raise any exception if it occurred
                    except Exception as e:
                        self.logger.error(f'Error during FASTQC execution: {e}')

                self.logger.info(f'Finished FASTQC for batch {batch_num}/{len(ftp_batches)}.')

        self.logger.info('Finished analyzing all FASTQ files.')

if __name__ == '__main__':
    args = create_parser().parse_args()
    config_file = args.config_file

    with open(config_file, 'r') as f:
        params = yaml.safe_load(f)

    logger = setup_logger(params['log_file'])
    create_directories(params, logger)
    logger.info('This script is only running FASTQC segment.')
    fastqc(params, logger).execute()
