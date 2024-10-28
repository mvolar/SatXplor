import subprocess
import argparse
import os
import glob
import utils.paths as paths




def create_blast_database(subject_file):

    subprocess.run(['makeblastdb', '-in', subject_file, '-dbtype', 'nucl', '-out', subject_file])
    
    return 

def cleanup_database_files(subject_file):
    # Clean up database files with the specified pattern
    pattern = f"./{subject_file}.n*"

    for file in glob.glob(pattern):
        print (f"removing databse file:{file}")
        os.remove(file)

def run_blast(query, subject, output, evalue=10):
    # Define the BLAST command
    blast_cmd = [
        'blastn',  # Replace with the appropriate BLAST command (e.g., blastp, blastx, etc.)
        '-query', query,
        '-db', subject,
        '-out', output,
        '-evalue', str(evalue),
        '-outfmt', str(6),
        '-max_target_seqs', str(10000),
        '-task', "blastn",
        '-num_threads', str(2)
    ]
    # Run the BLAST command
    subprocess.run(blast_cmd)

def blast_sequences():
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('sat_path', help='Path to the query sequence file in FASTA format')
    parser.add_argument('genome_path', help='Path to the subject file in FASTA format')
    
    args = parser.parse_args()
    print(args)

    create_blast_database(args.genome_path)
    run_blast(args.sat_path, args.genome_path, paths.BLAST_OUT_PATH,10)
    cleanup_database_files(args.genome_path)
        

if __name__ == '__main__':
    
    blast_sequences()
