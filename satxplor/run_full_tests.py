import shutil
import tarfile
import subprocess
from logging_config import logger


def extract_tar_gz(archive_path, extract_dir):
    with tarfile.open(archive_path, 'r:gz') as tar:
        for member in tar.getmembers():
            if member.isfile():
                # Extract the file directly to the specified directory
                tar.extract(member, path=extract_dir)

def copy_run_config(run_config_src, run_config_dest):
    shutil.copy2(run_config_src, run_config_dest)

def run_eusatxplor(eusatxplor_script):
    subprocess.run(['python', eusatxplor_script], check=True)

if __name__ == "__main__":
    logger.info("Setting up testing data environment.")


    archive_path = "testing_data/testing_data.tar.gz"
    extract_dir = "testing_data/"
    run_config_src = "testing_data/run_config_test.json"
    run_config_dest = "run_config.json"
    eusatxplor_script = "satxplor/controller.py"

    try:
        logger.info("Extracting data.")
        # Extract from the tar.gz archive
        extract_tar_gz(archive_path, extract_dir)

        # Copy and overwrite run_config.json
        copy_run_config(run_config_src, run_config_dest)

        # Run eusatxplor.py
        logger.info("Starting the test run.")
        run_eusatxplor(eusatxplor_script)
    except Exception as e:
        print("An error occurred:", e)