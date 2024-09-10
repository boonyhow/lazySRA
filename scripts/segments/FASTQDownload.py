import pandas as pd
import re
import os
import yaml
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
class download_fastq:
    def __init__(self, params, logger):
        self.params = params
        self.logger = logger

    def execute(self):    
        df = pd.read_csv(self.params['pheno_file'], sep=',')
        sra_ids = df.iloc[:, 0].tolist()
        out_dir = self.params['fastq_dir']
        
        if self.params['proxy']['exists'] == True:
            proxy = self.params['proxy']['proxy_server']
        else:
            proxy = None

        ftp_fps = []

        for sra_id in sra_ids:
            if not any(re.search(rf"^{sra_id}.*\.fastq\.gz$", f) for f in os.listdir(out_dir)):
                ftp_fps.extend(obtaining_ftp_fps(sra_id, self.logger, proxy))

        # Split the ftp_fps into batches of 50 if more than 50 files
        batch_size = self.params['fastq_download']['batch_size']
        if batch_size > 50:
            logger.warning('Maximum batch size for FTP requests on the SRA server is 50')
            batch_size = 50

        if len(ftp_fps) > batch_size:
            ftp_batches = [ftp_fps[i:i + batch_size] for i in range(0, len(ftp_fps), batch_size)]
        else:
            ftp_batches = [ftp_fps]

        for batch_num, batch in enumerate(ftp_batches, start=1):
            if len(batch) > 0:
                self.logger.info(f'Starting download batch {batch_num}/{len(ftp_batches)} with {len(batch)} files.')

                with ThreadPoolExecutor(max_workers=len(batch)) as executor:
                    futures = [executor.submit(download_ftp, ftp_fp, out_dir, self.logger, proxy) for ftp_fp in batch]
                    for future in as_completed(futures):
                        future.result()  # Ensures that exceptions are raised if any

                self.logger.info(f'Finished download batch {batch_num}/{len(ftp_batches)}.')

        self.logger.info('Finished downloading all FASTQ files.')

if __name__ == '__main__':
    args = create_parser().parse_args()
    config_file = args.config_file

    with open(config_file, 'r') as f:
        params = yaml.safe_load(f)

    logger = setup_logger(params['log_file'])
    create_directories(params, logger)
    logger.info('This script is only running Download FASTQ segment.')
    download_fastq(params, logger).execute()
