from Bio import SeqIO
import argparse
import utils.utils as utils
import utils.constants as constants
from utils.utils import extract_subsequences
from logging_config import logger
import glob

def extract_monomers(file_path: str, fasta_records: list) -> None:
    df = utils.read_gff_output(file_path,headers=False)
    grouped = df.group_by(["feature"])

    for group_key,group_df in grouped:

        if len(group_df) > 10000:
            logger.info(f"Large number of monomers detected for {group_key[0]}, downsampling to 10000 monomers.")
            group_df = group_df.sample(n=10000)


        subsequence_records = extract_subsequences(fasta_records, group_df)

        output_file = constants.SEQ_SAVE_PATH + group_key[0] + "_monomers.fasta"
        SeqIO.write(subsequence_records, output_file, 'fasta')



def extract_extended_arrays(file_path: str,fasta_records: list) -> None:
    df = utils.read_gff_output(file_path,headers=False)
    grouped = df.group_by(["feature"])

    for group_key,group_df in grouped:

    # Extract subsequences from FASTA records based on polars DataFrame
        subsequence_records = extract_subsequences(fasta_records, group_df,
        flanks="both",flank_size=500,squish_arrays=constants.SQUISH)

# Write the subsequences to a new FASTA file
        output_file = constants.SEQ_SAVE_PATH + group_key[0] + "_extended_arrays.fasta"
        SeqIO.write(subsequence_records, output_file, 'fasta')


def extract_arrays(file_path: str,fasta_records: list) -> None:
    df = utils.read_gff_output(file_path,headers=False)
    grouped = df.group_by(["feature"])

    #extract left flanks
    for group_key,group_df in grouped:

    # Extract subsequences from FASTA records based on polars DataFrame
        subsequence_records = extract_subsequences(fasta_records, group_df,
        flanks="none",flank_size=0)

# Write the subsequences to a new FASTA file
        output_file = constants.SEQ_SAVE_PATH + group_key[0] + "_nonextended_arrays.fasta"
        SeqIO.write(subsequence_records, output_file, 'fasta')
    
def create_monomer_dimers():
    files = glob.glob(constants.SEQ_SAVE_PATH + "*_monomers.fasta")  # Find all files ending with "_monomers.fasta"
    
    for file in files:
        output_file = file.replace("_monomers.fasta", "_monomer_dimers.fasta")  # Generate output filename
        logger.info(f"Creating synthetic dimers for {file}")
        with open(output_file, "w") as out_handle:
            for record in SeqIO.parse(file, "fasta"):
                sequence = record.seq
                dimer_sequence = sequence + sequence  # Concatenate the sequence with itself to create a dimer
                dimer_record = record
                dimer_record.seq = dimer_sequence
                SeqIO.write(dimer_record, out_handle, "fasta")



if __name__=="__main__":

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('genome_path', help='Path to the subject file in FASTA format')
    
    args = parser.parse_args()

    fasta_file = args.genome_path
    logger.info("Reading in the genome")
    fasta_records = list(SeqIO.parse(fasta_file, 'fasta'))
    logger.info(f"Extracting monomers from {constants.BLAST_GFF_PATH}")
    extract_monomers(constants.BLAST_GFF_PATH,fasta_records)

    logger.info(f"Creating monomer-dimers for kmer analysis {constants.BLAST_GFF_PATH}")
    create_monomer_dimers()

    logger.info(f"extracting array sequences from from {constants.ARRAYS_OUT_PATH}")
    extract_arrays(constants.ARRAYS_OUT_PATH,fasta_records)

    logger.info(f"extracting array sequences with flanks from from {constants.ARRAYS_OUT_PATH}")
    extract_extended_arrays(constants.ARRAYS_OUT_PATH,fasta_records)

