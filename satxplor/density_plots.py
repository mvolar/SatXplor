
from utils.utils import read_blast_output
from scipy.stats import gaussian_kde
import polars as pl
import utils.constants as constants
import matplotlib.pyplot as plt
import numpy as np
from logging_config import logger

df = read_blast_output(constants.BLAST_OUT_PATH)



mlendf = (df
          .group_by('query')
          .agg([pl.col('al_len')
                .max()])
                .rename({"al_len":"max_len"}))

df = df.join(mlendf,on="query")

df = df.with_columns(
    qcovhsp = pl.col("al_len")/pl.col("max_len")
)

df = df.filter(pl.col("qcovhsp")>0.25)

grouped = df.group_by(['query'])

for group_key, group_df in grouped:
    try:
        nbins=constants.NKERNEL_BINS
        x,y = group_df["qcovhsp"]*100,group_df["perc_id"]
        k = gaussian_kde([x,y])

        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    #Plot the 2D density plot using Matplotlib
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto',snap=True,rasterized=True)
        plt.title(group_key[0])
        plt.xlabel('Query coverage')
        plt.ylabel('Percentage identity')
        plt.savefig(constants.PIC_SAVE_ROOT + str(group_key[0]) + "_density.png")
        plt.clf()
    except Exception as e:
        logger.warning(f"Error with density plots of {group_key[0]} most likely due to low number of monomers in KDE approximation, this is a non essential error, but we suggest closer examination of the satDNA in question. {e}")
