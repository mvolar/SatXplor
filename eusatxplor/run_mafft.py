import os
import subprocess
import multiprocessing
import utils.constants as constants
from logging_config import logger

def run_mafft(input_file, output_file):
    # Get the number of available threads
    available_threads = multiprocessing.cpu_count()
    
    # Calculate half of the available threads
    half_threads = max(1, available_threads // 2)

    # Build the MAFFT command
    mafft_command = [
        'mafft',
        '--adjustdirection',
        '--reorder',
        '--quiet',
        '--thread',
         str(half_threads),
        input_file
    ]
    try: 
        result = subprocess.run(mafft_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,check=True)
    # Redirect output to the specified file
        if result.stderr:
            print(f"MAFFT STDERR for {input_file}:\n{result.stderr}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT for {input_file}. Error message: {e}")
        raise e
    
    # Save output to the specified file
    with open(output_file, 'w') as output:
        output.write(result.stdout)


def main():
    # List of input files for monomer alignment
    input_files = [file for file in os.listdir(constants.SEQ_SAVE_PATH) if file.endswith("_monomers.fasta")]
    for input_file in input_files:
        # Construct output file name by adding '_aligned' suffix
        output_file = os.path.splitext(input_file)[0] + '_aligned.fasta'
        input_path = constants.SEQ_SAVE_PATH + input_file
        output_path = constants.SEQ_SAVE_PATH + output_file
        

        # Run MAFFT for each file
        run_mafft(input_path, output_path)
        logger.info(f"Alignment completed for {input_path}. Output saved to {output_path}")
    
        #list of input files for flank alignment
    input_files = [file for file in os.listdir(constants.FLANKS_SAVE_ROOT   ) if file.endswith("_flanks.fasta")]
    
    for input_file in input_files:
        # Construct output file name by adding '_aligned' suffix
        output_file = os.path.splitext(input_file)[0] + '_aligned.fasta'
        input_path = constants.FLANKS_SAVE_ROOT + input_file
        output_path = constants.FLANKS_SAVE_ROOT + output_file
    
        # Run MAFFT for each file
        try:
            run_mafft(input_path, output_path)
            logger.info(f"Alignment completed for {input_path}. Output saved to {output_path}")
        except:
            logger.warn("Error in mafft")
    


    input_files = [file for file in os.listdir(constants.FLANKS_SAVE_ROOT   ) if file.endswith("_microhomology.fasta")]
    
    for input_file in input_files:
        # Construct output file name by adding '_aligned' suffix
        output_file = os.path.splitext(input_file)[0] + '_aligned.fasta'
        input_path = constants.FLANKS_SAVE_ROOT + input_file
        output_path = constants.FLANKS_SAVE_ROOT + output_file
    
        # Run MAFFT for each file
        try:
            run_mafft(input_path, output_path)
            logger.info(f"Alignment completed for {input_path}. Output saved to {output_path}")
        except:
            logger.warn("Error in mafft")




if __name__ == "__main__":
    main()