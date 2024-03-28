from Bio import SeqIO

def sanitize_and_filter_sequences(input_file, output_file, min_length=100000):
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence_length = len(record.seq)
            if sequence_length >= min_length:
                sanitized_name = record.id.replace("_", "")  # Remove underscores from the name
                record.id = sanitized_name
                record.description = ""  # Clear the description
                SeqIO.write(record, out_handle, "fasta")

# Replace 'input.fasta' with the actual filename of your input FASTA file
