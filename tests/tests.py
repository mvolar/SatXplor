import subprocess

def is_ncbi_blast_installed():
    try:
        # Run the blastn command with the --version option to check if BLAST is installed
        subprocess.run(['blastn', '-version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False



def run_r_script(script_path):
    try:
        # Run the R script and capture both stdout and stderr
        result = subprocess.run(['Rscript', script_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("R script executed successfully.")
        return True
    except subprocess.CalledProcessError as e:
        # Print the error message if the execution fails
        print("Error: Failed to execute the R script.")
        print("Error message:", e.stderr.decode())
        return False


def is_mafft_installed():
    try:
        subprocess.run(['mafft', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False

if __name__ == "__main__":

    print("Running installation tests.")

    script_path = "./tests/rlib_import.R"

    if run_r_script(script_path):
        print("R is setup properly.")
    else:
        print("There was an error. Some libraries missing or not installed.")

    if is_mafft_installed():
        print("MAFFT is installed.")
    else:
        print("MAFFT is not installed.")
    

    if is_ncbi_blast_installed():
        print("NCBI BLAST is installed.")
    else:
        print("NCBI BLAST is not installed.")
