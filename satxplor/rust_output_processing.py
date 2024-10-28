import os
import polars as pl
import utils.paths as paths
from utils.constant_loader import constants as constants
from utils.utils import read_gff_output, extract_subsequences
from Bio import SeqIO
import argparse

from logging_config import logger

def extract_array_flanks(df: pl.DataFrame,fasta_records: list,
                         flank_size: int =constants.FLANK_SIZE,flank_type: str="flanks") -> None:  


    grouped = df.group_by(["feature"])

    for group_key,group_df in grouped:
    
        #extract left flanks
        left_flanks = extract_subsequences(fasta_records, group_df,flanks="left",flank_size=flank_size)
        #extract right flanks
        right_flanks = extract_subsequences(fasta_records, group_df,flanks="right",flank_size=flank_size)
        #merge them together
        left_flanks.extend(right_flanks)
        
        #Write the subsequences to a new FASTA file
        output_file = paths.FLANKS_SAVE_ROOT + group_key[0] + "_"+ flank_type +".fasta"
        SeqIO.write(left_flanks, output_file, 'fasta')

def load_in_df_from_list_of_files(file_list: list,squish=False):
    
    df_list = []

    if squish:
        for i in file_list:
            name = i.split("/")[-1]
            name = name[:-4]
            name = name.split("_")[0:4]
            name = "_".join(name)
            tmp = pl.read_csv(i,dtypes={"index":int,
            "actual":float,
            "roll_mean":float},separator="\t")
            tmp = tmp.with_columns(
                name = pl.lit(name)
            )
            
            df_list.append(tmp)
    else: 
        for i in file_list:
            name = i.split("/")[-1]
            name = name[:-4]
            name = name.split("_")[0:3]
            name = "_".join(name)
            tmp = pl.read_csv(i,dtypes={"index":int,
            "actual":float,
            "roll_mean":float},separator="\t")
            tmp = tmp.with_columns(
                name = pl.lit(name)
            )
            df_list.append(tmp)
    return df_list

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('genome_path', help='Path to the subject file in FASTA format')

    args = parser.parse_args()
    
    fasta_records = list(SeqIO.parse(args.genome_path, 'fasta'))
    
    logger.info("\n Loaded the fasta")


    file_list = os.listdir(paths.KMER_ANALYSIS_OUT + "data")


    file_list = [paths.KMER_ANALYSIS_OUT + "data/" + string for string in file_list]

    df_list = load_in_df_from_list_of_files(file_list=file_list,squish=constants.SQUISH)

    #find the min max values by unique array name in the filtered data
    df_tot = pl.concat(df_list).filter(
        pl.col("roll_mean")<5.0
    ).group_by('name').agg(
        min_value=pl.col('index').min(),
        max_value=pl.col('index').max()
    )

    df_tot = df_tot.with_columns(
        query = pl.col("name").map_elements(lambda x: x.split("_")[1])
    )
    df_tot =df_tot.rename({
            "min_value" : "start",
            "max_value" : "end",
            "query" : "feature",
            "name" : "seqnames"
        })

    df_tot = df_tot.with_columns(
        [
            (pl.lit("feature")).alias("group"),
            (pl.lit("*")).alias("strand"),
            (pl.lit("EuSatXplor")).alias("source"),
            (pl.lit(".")).alias("frame"),
            (pl.lit("100")).alias("score")
        ]
    )

    df_tot = df_tot.select(
                    ["seqnames","source","feature","start","end","score","strand","frame","group"]
        )

    df_tot.write_csv(paths.KMER_ANALYSIS_OUT + "final_array_table_array_maps.tsv",
                    include_header=False,separator="\t")
    
    logger.info("Created the real edges annotations")
    

    
    #read in the original arrays
    df = read_gff_output(paths.ARRAYS_OUT_PATH,headers=False)

    df_real_edges = df_tot
    if constants.SQUISH:
        df = df.with_columns(
            name = pl.col("seqnames")+"_" + pl.col("group") +"_" + pl.col("start").cast(pl.String)+"_" + pl.col("end").cast(pl.String)
        ).with_columns(
            len= pl.col("end") -pl.col("start")
        )


        df = df.join(df_real_edges,left_on="name",right_on="seqnames")

        df = df.with_columns(
            real_start=pl.col("start") - (500-pl.col("start_right"))).with_columns(
            real_end=pl.when(pl.col("len")<5000)
                        .then(pl.col("end_right") - pl.col("start_right") + pl.col("real_start")) 
                        .otherwise(pl.col("end") - (4500-pl.col("end_right")))
        ).select(
                ["seqnames","source","feature","real_start","real_end","score","strand","frame","group"]
        ).rename(
            {
                "real_start" : "start",
                "real_end" :"end"
            }
        )
        
        df.write_csv(paths.KMER_ANALYSIS_OUT + "real_edges_in_genome.gff",include_header=False,separator="\t")
        
        
        logger.info("Extracted the real flanks")
        
        extract_array_flanks(df=df,fasta_records=fasta_records)

        logger.info("Extracted the microhomology regions")
        extract_array_flanks(df=df,fasta_records=fasta_records,
                            flank_size=20,flank_type="microhomology")
    
    else: 
        df = df.with_columns(
            name = pl.col("seqnames")+"_" + pl.col("group") +"_" + pl.col("start").cast(pl.String)
        )

        df = (df.join(df_real_edges,left_on="name",right_on="seqnames")
        .with_columns(
            real_start=pl.col("start") - (500-pl.col("start_right")))
        .with_columns(
            real_end=pl.col("end_right") - pl.col("start_right") + pl.col("real_start")
        )
        .select(
                ["seqnames","source","feature","real_start","real_end","score","strand","frame","group"]
        )
        .rename(
            {
                "real_start" : "start",
                "real_end" :"end"
            }
        )
        )
        
        df.write_csv(paths.KMER_ANALYSIS_OUT + "real_edges_in_genome.gff",include_header=False,separator="\t")
        
        
        logger.info("Extracted the real flanks")
        
        extract_array_flanks(df=df,fasta_records=fasta_records)

        logger.info("Extracted the microhomology regions")
        extract_array_flanks(df=df,fasta_records=fasta_records,
                            flank_size=20,flank_type="microhomology")






