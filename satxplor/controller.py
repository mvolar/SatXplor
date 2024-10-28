import os
import shutil
import subprocess
from utils.constant_loader import constants as constants
import utils.utils as utils
import utils.paths as paths
import argparse
import sys
from logging_config import logger
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from r_script_runners import (run_pca_script,
                                  run_flank_distances,
                                  run_microhomology,
                                  run_networks_script )
import json
from preprocess import sanitize_and_filter_sequences

def save_constants(data):
    with open('satxplor/utils/constants.json', 'w') as f:
        json.dump(data, f, indent=4)

def update_constants(args):
    # Update constants with new values from command-line arguments (if provided)
    constants = {
        'MAX_PLOT_LEN': args.max_plot_len,
        'HISTOGRAM_BIN_WIDTH': args.histogram_bin_width,
        'NKERNEL_BINS': args.nkernel_bins,
        'FLANK_SIZE': args.flank_size,
        'ARRAY_HOR_PERC': args.array_hor_perc,
        'MONOMER_NUMBER': args.monomer_number,
        'CONTIG_FILTER': args.contig_filter,
        'DIMENSION_RED_MODE': args.dimension_red_mode,
        'PERC_ID_FILTER': args.perc_id_filter,
        'QCOVHSP_FILTER': args.qcovhsp_filter,
        'SQUISH': args.squish,
        'SQUISH_LEN': args.squish_len
    }

    # Save updated constants back to constants.json
    save_constants(constants)

    # Print updated constants
    print(f"Updated Constants:\n{json.dumps(constants, indent=4)}")


with open("run_config.json", "r") as f:
    config = json.load(f)
checkpoints = ["nothing",
               "blast",
                   "arrays",
                   "density",
                   "extract",
                   "filenames",
                   "kmercount",
                   "output",
                   "alignment",
                   "PCA"]

def audit_checkpoints() -> str: 
    checkpoint_file = "./checkpoints.tsv"
    # Check if checkpoint file exists or is empty
    if not os.path.exists(checkpoint_file):
        with open(checkpoint_file, mode='w') as file:
            file.write("Checkpoint\tDatetime\n")
        return("nothing")
    else:
        # Read the content of the last row
        with open(checkpoint_file, mode='r') as file:
            lines = file.readlines()
            last_row = lines[-1].strip()

            # Extract the last checkpoint name and datetime
            last_checkpoint, last_datetime = last_row.split('\t')
            
            if last_checkpoint=="Checkpoint":
                logger.info("Checkpoint file was created, but nothing was in it, starting from scratch.")
                return("nothing")
            else: 
                logger.info(f"Last checkpoint is: {last_checkpoint} at {last_datetime}")
                return(last_checkpoint)
    

def create_checkpoint(checkpoint_name):
    checkpoint_file = "./checkpoints.tsv"
    
    # Append checkpoint and datetime to the checkpoint file
    with open(checkpoint_file, mode='a') as file:
        file.write(f"{checkpoint_name}\t{datetime.now()}\n")


#new blast will always reinitialize folders
        
def initialize_preprocess_blast():
    # Create or clear the 'results' folder
    logger.info('Removing data in existing results and tmp directories.')
    results_folder = './results'
    if os.path.exists(results_folder):
        shutil.rmtree(results_folder)

        
    os.makedirs(results_folder)

    #data saving roots
    os.makedirs(paths.DATA_SAVE_ROOT)
    os.makedirs(paths.TABLE_SAVE_ROOT)

    #sequence saving roots
    os.makedirs(paths.SEQ_SAVE_PATH)
    os.makedirs(paths.FLANKS_SAVE_ROOT)

    #picture saving roots
    os.makedirs(paths.PIC_SAVE_ROOT)
    os.makedirs(paths.PCA_UMAP_SAVE_ROOT)
    os.makedirs(paths.DISTANCE_SAVE_ROOT)
    os.makedirs(paths.NETWORKS_SAVE_ROOT)
    os.makedirs(paths.MICROHOMOLOGY_SAVE_ROOT)
    
    logger.info('Data deleted and folders created.')

    logger.info(
        f"Running preprocessing script for the input.\n  Removing all sequences shorter than {constants.CONTIG_FILTER}.")
    
    input_file = config["INPUT_GENOME_PATH"]
    output_file = config["GENOME_PATH"]
    
    sanitize_and_filter_sequences(input_file, output_file,min_length=constants.CONTIG_FILTER)

    input_file = config["SAT_RAW"]
    output_file = config["SAT_FASTA_PATH"]
    
    
    logger.info("Sanitizing input names")
    
    sanitize_and_filter_sequences(input_file, output_file,min_length=0)


    logger.info( f"Output sequences written to { config['GENOME_PATH'] }")

    command = ['python3', './satxplor/blast.py',
               config["SAT_FASTA_PATH"],
               config["GENOME_PATH"]
               ]
    
    print(f"Blast command:{command}")
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        logger.warning(f"Error. Return code:{ result.returncode}")
        print("STDERR:\n", result.stderr)
        quit()


def run_create_arrays_script():
    # Run the 'create_arrays.py' script

    command = ['python3', './satxplor/create_arrays.py']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()
    

def run_density_plots():
    # Run the 'create_arrays.py' script
    command = ['python3', './satxplor/density_plots.py']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        #this is an optional output, if error, continue normaly
        logger.warning(f"Error with density plots. Return code:{ result.returncode}")
        
    

def run_extract_naive():
    # Run the 'create_arrays.py' script
    command = ['python3', './satxplor/extract_naive.py',
               config["GENOME_PATH"]]
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_filename_mapping():
    # Run the 'create_arrays.py' script

    command = ['python3', './satxplor/filename_mapping.py']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_kmer_edge_finder():
    # Run the 'create_arrays.py' script
    command = ['./satxplor/executables/kmer_edge_finder']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_rust_output_processing():
    # Run the 'create_arrays.py' script
    command = ['python', './satxplor/rust_output_processing.py',
                              config["GENOME_PATH"]]

    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_mafft():
    # Run the 'create_arrays.py' script
    command = ['python3', './satxplor/run_mafft.py']
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        logger.warning(f"Error. Return code:{ result.returncode}. Error in running MAFFT for some sequences. ")

def run_pca():
    search_string = 'monomers_aligned'
    matching_files = utils.get_files_containing_string(paths.SEQ_SAVE_PATH, search_string)

    with ProcessPoolExecutor() as executor:
        executor.map(run_pca_script, matching_files)
    
def run_networks():
    search_string = 'matrix.csv.gz'
    matching_files = utils.get_files_containing_string(paths.PCA_UMAP_SAVE_ROOT, search_string)

    
    with ProcessPoolExecutor() as executor:
        executor.map(run_networks_script, matching_files)

def run_flanks():
    search_string = 'flanks_aligned'
    matching_files = utils.get_files_containing_string(paths.FLANKS_SAVE_ROOT, search_string)

    with ProcessPoolExecutor() as executor:
        executor.map(run_flank_distances, matching_files)


def main():
    logger.info("Python Version:")
    logger.info(sys.version)

    logger.info(config)
    
    if os.path.exists(config["FINAL_RESULTS_DIR"]):
        logger.info("The output folder already exists")
        if not config["OVERWRITE"]:

            logger.error(f"Destination folder '{config[ 'FINAL_RESULTS_DIR' ] }' already exists. but OVERWRITE is off. Exiting the application.")
        else: 
            logger.info("Since OVERWRITE is on, trying to delete the folder!")
            shutil.rmtree(config["FINAL_RESULTS_DIR"])
            logger.info("Folder successfully deleted")



    check_point = audit_checkpoints()

    starting_id = checkpoints.index(check_point)

    parser = argparse.ArgumentParser(description='Update constants in constants.json')

    # Hard-coded default values in argparse
    parser.add_argument('--max_plot_len', type=int, default=5000, help='Maximum plot length')
    parser.add_argument('--histogram_bin_width', type=int, default=50, help='Histogram bin width')
    parser.add_argument('--nkernel_bins', type=int, default=100, help='Number of kernel bins')
    parser.add_argument('--flank_size', type=int, default=500, help='Flank size')
    parser.add_argument('--array_hor_perc', type=float, default=0.01, help='Array horizontal percentage')
    parser.add_argument('--monomer_number', type=int, default=4, help='Monomer number')
    parser.add_argument('--contig_filter', type=int, default=100000, help='Contig filter')
    parser.add_argument('--dimension_red_mode', type=str, default='both', help='Dimension reduction mode')
    parser.add_argument('--perc_id_filter', type=float, default=70.0, help='Percentage identity filter')
    parser.add_argument('--qcovhsp_filter', type=float, default=70.0, help='Query coverage filter')
    parser.add_argument('--squish', type=bool, default=True, help='Squish flag')
    parser.add_argument('--squish_len', type=int, default=2500, help='Squish length')

    args = parser.parse_args()

    update_constants(args)


    logger.info("The output folder already exists")
    fun_list = [initialize_preprocess_blast,
                run_create_arrays_script,
                run_density_plots,
                run_extract_naive,
                run_filename_mapping,
                run_kmer_edge_finder,
                run_rust_output_processing,
                run_mafft,
                run_flanks,
                run_pca,
                run_networks,
                run_microhomology]
    
    for i in range(starting_id,len(checkpoints)+2):
        fun_list[i]()
        try:
            create_checkpoint(checkpoints[i])
        except IndexError:
            pass
        except Exception as e:
            print("Stopping execution due to:", e)

    
    
    logger.info(f"Running finished at {datetime.now()}, deleting the checkpoints file. Copying results file to output directory {config['FINAL_RESULTS_DIR']}")
    os.remove("checkpoints.tsv")


        

# Copy the entire folder from source to destination
    shutil.copytree("./results", config["FINAL_RESULTS_DIR"])
    utils.move_file(config["GENOME_PATH"],config['FINAL_RESULTS_DIR'])
    utils.move_file(config["SAT_FASTA_PATH"],config['FINAL_RESULTS_DIR'])
    

if __name__ == '__main__':

    main()