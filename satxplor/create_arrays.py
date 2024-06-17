import numpy as np
import polars as pl
from plotnine import ggplot, aes, geom_histogram, ggtitle, theme_bw
import pandas as pd
import json
import utils.constants as constants
from utils.utils import read_blast_output, convert_df_to_gff
from logging_config import logger



#temporary function used only in this scripts
def find_and_group_overlapping(group):
    intervals = []
    strand = group["strand"][0]
    current_start, current_end = group['new_start'][0], group['new_end'][0]
    for i in range(1, len(group)):
        if group['new_start'][i] <= current_end:
            # Overlapping intervals, extend the current interval
            current_end = max(current_end, group['new_end'][i])
        else:
            # Non-overlapping interval, start a new interval
            intervals.append((current_start, current_end))
            current_start, current_end = group['new_start'][i], group['new_end'][i]

    intervals.append((current_start, current_end))
    return pl.DataFrame({
        'seqnames': [group['subject'][0]] * len(intervals),
        'source': ["EuSatXplor"] * len(intervals),
        'feature' : [group['query'][0]] * len(intervals),
        'start': [start for start, _ in intervals],
        'end': [end for _, end in intervals],
        'score': [100] * len(intervals),
        'strand' : [strand] * len(intervals),
        'frame' : ["."] * len(intervals),
        'group' : [group['query'][0]] * len(intervals),
    })


if __name__ == '__main__':
    
    logger.info(f"Finding HORs in BLAST output: {constants.BLAST_OUT_PATH}")
    df = read_blast_output(constants.BLAST_OUT_PATH)
    df = df.with_columns(
        pl.when(df['s_start'] > df["s_end"]).then(df['s_end']).otherwise(df['s_start']).alias('new_start'),
        pl.when(df['s_start'] > df["s_end"]).then(df['s_start']).otherwise(df['s_end']).alias('new_end'),
        pl.when(df['s_start'] > df["s_end"]).then(pl.lit("-")).otherwise(pl.lit("+")).alias('strand'),
    )
    
    
    mlendf = (df
          .group_by('query')
          .agg([pl.col('q_end')
                .max()])
                .rename({"q_end":"max_len"}))
    df = df.join(mlendf,on="query")

    df = df.with_columns(
    qcovhsp = pl.col("al_len")*100/pl.col("max_len")
    )

    shifted_start = (
        df
    .shift(-1)
    .select(pl.col("new_start"))
    .to_series()
    )

    df = df.with_columns(
        next_start = shifted_start,
    )

    df = df.with_columns(
        distance = pl.col("next_start") - pl.col("new_end")
    )

    df = df.filter(
        ((pl.col("perc_id") >constants.PERC_ID_FILTER) &  (pl.col("qcovhsp")>constants.QCOVHSP_FILTER))
    )
     
    
    df = df.sort(['query','subject', 'new_start'])

    # Group by 'chromosome' and apply the function to find and group overlapping intervals

    grouped = df.group_by(['query'])

    ef_list = []
    name_list = []
    df_list = []
    pl.Config.set_tbl_rows(500)
    logger.info("Monomer statistics:")
    summary_stats = (
    df.group_by("query")
    .agg(
        pl.col("perc_id").mean().alias("perc_id_mean"),
        pl.col("perc_id").count().alias("perc_id_count"),
        pl.col("perc_id").median().alias("perc_id_median"),
        pl.col("perc_id").quantile(0.25).alias("perc_id_1q"),
        pl.col("perc_id").quantile(0.75).alias("perc_id_3q"),
        pl.col("qcovhsp").mean().alias("aln_len_mean"),
        pl.col("qcovhsp").count().alias("aln_len_count"),
        pl.col("qcovhsp").median().alias("aln_len_median"),
        pl.col("qcovhsp").quantile(0.25).alias("aln_len_1q"),
        pl.col("qcovhsp").quantile(0.75).alias("aln_len_3q"),
    )
)
    print(summary_stats)


    # Iterate through groups
    for group_key, group_df in grouped:
        try:
            logger.info(f"Lengths of group {group_key[0]} is {len(group_df)}")
        

            monomer_len = group_df.select("max_len").row(0)[0]

            tmp_df = group_df.filter(
            ((pl.col("distance")<constants.MAX_PLOT_LEN) & (pl.col("distance")>-10)) 
        )

        # Create a histogram using Plotly Express
            plot = (ggplot(tmp_df, aes(x='distance')) +
                geom_histogram(alpha=0.7,binwidth=constants.HISTOGRAM_BIN_WIDTH,color="black",
                            fill="#e5c8d6") +
                ggtitle(group_key[0]) +
                theme_bw() 
            )
            plot.save(f"results/pictures/{group_key[0]}_HOR.png",verbose=False)



            query_len = group_df.select(pl.col("q_end")).max().to_series()
            query_count = group_df.shape[0]
            # Find peaks

            data = (group_df
                    .filter(pl.col("distance")<constants.MAX_PLOT_LEN)
                    .select(pl.col("distance"))
                    .to_series())
            bin_length = 100
            bins = np.arange(0, np.nanmax(data) + bin_length, bin_length)
            hist, edges = np.histogram(data, bins=bins)


            # Create a Pandas Series with counts
            counts_series = pd.Series(hist, index=bins[1:])
            extension_factor = query_len[0]
            try:
                logger.info(f"Finding peaks for {group_key}")

                x = (counts_series[counts_series>constants.ARRAY_HOR_PERC*query_count]).idxmax()

                extension_factor = int(x) + query_len[0]
            except Exception as e:
                logger.warning(f"No peaks were found for {group_key}, setting default extension factor, reported error {e}")
                
            logger.info(f"Extension factor for {group_key} is {extension_factor}")
            ef_list.append(extension_factor)
            name_list.append(group_key[0])

            group_df = group_df.with_columns(
            new_end = pl.col("new_end") + extension_factor
        )

            df_overlapped = group_df.group_by('subject','query').map_groups(find_and_group_overlapping)
            df_overlapped = df_overlapped.with_columns(
                end=pl.col("end")-extension_factor
            )

            logger.info(f"Number of arrays before filtering: {len(df_overlapped)}")

            df_overlapped = df_overlapped.filter(
            (pl.col("end") - pl.col("start") > constants.MONOMER_NUMBER*monomer_len)
        )
            
            logger.info(f"Number of arrays after filtering: {len(df_overlapped)}")


            df_list.append(df_overlapped)
        except Exception as e:
            logger.warning(f"Error creating arrays for {group_key[0]} error message: {e}")

    data_dict = dict(zip(name_list,
                        ef_list))
    
    file_path = './results/data/extension_factors.json'

    # Write the JSON string to a file
    with open(file_path, 'w') as json_file: 
        json.dump(data_dict, json_file, indent=2) 

    #write the final output table
    out_table = pl.concat(df_list)
    


    df = convert_df_to_gff(df)
    df.write_csv(constants.BLAST_GFF_PATH, separator="\t",include_header=False )
    

    out_table.write_csv(constants.ARRAYS_OUT_PATH, separator="\t",include_header=False)
    logger.info(f"Extenstion factors written {file_path}")
    logger.info(f"Arrays written to {constants.ARRAYS_OUT_PATH}")
    logger.info(f"Unfiltered monomer annots written to {constants.BLAST_GFF_PATH}")

    logger.info("Array statistics:")
    
    print(
        (out_table
        .with_columns(
        width = pl.col("end") - pl.col("start")
        )
        .group_by(["seqnames","feature"])
        .agg(
        pl.col("width").mean().alias("mean"),
        pl.col("width").count().alias("count"),
        pl.col("width").median().alias("median"),
        pl.col("width").quantile(0.25).alias("1q"),
        pl.col("width").quantile(0.75).alias("3q"),
        )
        )
    )
    #logger.info(f"Monomer statistics: {out_table.with_columns(
    #    width = pl.col("end") - pl.col("start")
    #).describe()}")
