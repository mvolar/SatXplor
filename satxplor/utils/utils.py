import polars as pl
from Bio.SeqRecord import SeqRecord
import os
import utils.constants as constants
import shutil

def read_blast_output(file_path):
    # Read the BLAST output file into a DataFrame
    df = pl.read_csv(file_path,skip_rows=0,separator="\t",has_header=False, dtypes={'column_12': float})  # Adjust delimiter based on your file
    df.columns = ["query","subject","perc_id","al_len","mismatches","gap_opens","q_start","q_end","s_start","s_end","evalue","score"]
    return df

def read_gff_output(file_path,headers=True):
    # Read the gff with no headers
    if headers:
        df = pl.read_csv(file_path,skip_rows=0,separator="\t",has_header=True)  # Adjust delimiter based on your file
    else:
        df = pl.read_csv(file_path,skip_rows=0,separator="\t",has_header=False)  # Adjust delimiter based on your file
        df.columns = ["seqnames","source","feature","start","end","score","strand","frame","group"]
    return df

#TODO add filtering step for perc id and qcov hsp to be default at 70
def convert_df_to_gff(df):
    # Read the BLAST output file into a DataFrame
    
    df = df.rename({
        "new_start" : "start",
        "new_end" : "end",
        "query" : "feature",
        "subject" : "seqnames"
    })

    df = df.with_columns(
    [
        (pl.col("feature")).alias("group"),
        (pl.lit("BLAST")).alias("source"),
        (pl.lit(".")).alias("frame"),
    ]
)
    df = df.select(
                   ["seqnames","source","feature","start","end","score","strand","frame","group"]
    )
    return df


# Function to extract subsequences based on start and end positions
def extract_subsequences(records, df:pl.DataFrame,flanks="none",flank_size=0, squish_arrays: bool=False) -> list[SeqRecord]:
    result_records: list = []

    for row in df.rows(named=True):
        header: str = row['seqnames']
        prefix: str = ""

        if flanks =="left":
            start: int = row['start'] - flank_size 
            end: int = row['start']
            prefix = "l_"
        if flanks =="right":
            start: int = row['end']  
            end: int = row['end'] + flank_size
            prefix = "r_"
        if flanks =="both":
            start: int = row['start'] - flank_size 
            end: int = row['end'] + flank_size
        if flanks =="none":
            start: int = row['start'] 
            end: int = row['end'] 

        if squish_arrays:
            for record in records:
                if record.id == header:
                    subsequence = record.seq[start:end]
                    sequence_length = len(subsequence)
                    id_string: str = prefix +  str(row["seqnames"])+ "_" + str(row["feature"])+"_" + str(row["start"]) + "_" +  str(row["end"])
                    if sequence_length > 5000:
                        merged_sequence = subsequence[:constants.SQUISH_LEN] + subsequence[-constants.SQUISH_LEN:]
                        result_records.append(SeqRecord(merged_sequence, id=id_string, description=''))
                    else:
                        result_records.append(SeqRecord(subsequence, id=id_string, description=''))

        
        
        else :
            for record in records:
                if record.id == header:
                    subsequence = record.seq[start:end]
                    id_string: str = prefix +  str(row["seqnames"])+ "_" + str(row["feature"])+"_" + str(row["start"])
                    result_records.append(SeqRecord(subsequence, id=id_string, description=''))
                    
    filtered_seq_records = [record for record in result_records if len(record.seq) > 0]
    return filtered_seq_records


def get_files_containing_string(directory, search_string):
    matching_files = []
    
    # Traverse through the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Check if the search string is present in the file name
            if search_string in file:
                # If found, add the file to the list
                matching_files.append(os.path.join(root, file))

    return matching_files

def move_file(src, dst):
    try:
        shutil.move(src, dst)
    except Exception as e:
        print(f"Error: Failed to move file '{src}' to '{dst}': {e}")
