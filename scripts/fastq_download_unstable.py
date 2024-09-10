import yaml
# from src.utils.utilities import *
# from src import download_fastq, generate_genome_index, generate_bam_files, featurecounts_generate_counts_matrix
from utils.utilities import *
from segments.fastqc import fastqc
from segments.FASTQDownload import download_fastq
from segments.GenomeIndex import generate_genome_index
from segments.STARAlignment import generate_bam_files
from segments.featureCounts import featurecounts_generate_counts_matrix

# Main function to handle the process
def main(params, logger):
    genome_job_id = []
    bam_job_ids = []

    create_directories(params, logger)

    # Download FASTQ files
    if params['fastq_download']['initiate']:
        logger.info('Initiating FASTQ download')
        download_fastq(params, logger).execute()

    if params['fastqc']['initiate']:
        fastqc(params, logger).execute()
    
    # Genome indexing (if applicable)
    if params['genome_indexing']['initiate']:
        logger.info('Initiating Genome Index Generation')
        genome_job_id = generate_genome_index(params, logger).execute()
    
    # BAM file generation
    if params['star_alignment']['initiate']:
        logger.info('Starting STAR alignment process.')
        bam_job_ids = generate_bam_files(params, logger).execute()

    # Counts matrix generation, dependent on BAM files
    if params['counts_generations']['initiate']:
        logger.info('Generating counts matrix.')
        featurecounts_generate_counts_matrix(params, logger, dependency=":".join(bam_job_ids)).execute()

# Main entry point
if __name__ == '__main__':
    args = create_parser().parse_args()
    config_file = args.config_file

    with open(config_file, 'r') as f:
        params = yaml.safe_load(f)

    logger = setup_logger(params['log_file'])
    if params['slurm']:
        logger.info('This script is SLURM configured')
    main(params, logger)
