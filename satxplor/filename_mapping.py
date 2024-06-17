import os
import json
import utils.constants as constants

from logging_config import logger
def match_monomer_array_pairs(folder_path):
    monomer_files = {}
    array_files = {}

    # List files in the given folder
    files = os.listdir(folder_path)

    # Iterate through files
    for file_name in files:
        if file_name.endswith("_monomer_dimers.fasta"):
            monomer_files[file_name] = None
        elif file_name.endswith("_arrays.fasta"):
            array_files[file_name] = None

    # Match monomer and array files
    pairs = {}
    for monomer_file in monomer_files:
        prefix = monomer_file.replace("_monomer_dimers.fasta", "")
        array_file = f"{prefix}_extended_arrays.fasta"
        if array_file in array_files:
            monomer_path = os.path.join(folder_path, monomer_file)
            array_path = os.path.join(folder_path, array_file)
            pairs[monomer_path] = array_path

    return pairs

def write_to_json(data, json_file):
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)

def main():
    folder_path = constants.SEQ_SAVE_PATH
    output_json_file = "pairs.json"

    pairs = match_monomer_array_pairs(folder_path)
    logger.info("Monomer-Array Pairs:")
    for monomer_path, array_path in pairs.items():
        print(f"{monomer_path} -> {array_path}")

    write_to_json(pairs, output_json_file)
    logger.info(f"Pairs written to {output_json_file}")

if __name__ == "__main__":
    main()
