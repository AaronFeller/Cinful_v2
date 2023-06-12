from Bio import SeqIO
import multiprocessing

def translate_record(record):
    # Translate the DNA sequence into protein sequence
    dna_seq = record.seq
    protein_seq = dna_seq.translate()
    # Create a new record with the same id and description as the original one
    protein_record = SeqIO.SeqRecord(protein_seq, id=record.id, description=record.description)
    return protein_record

def process_records(records):
    # Process a batch of records and return the translated records
    translated_records = []
    for record in records:
        translated_record = translate_record(record)
        translated_records.append(translated_record)
    return translated_records

def write_records(records, output_file):
    # Write the translated records to the output file in fasta format
    with open(output_file, "w") as out_handle:
        SeqIO.write(records, out_handle, "fasta")

def process_file(input_file, output_file, threads):
    # Read all the records from the input file
    with open(input_file, 'r') as input_handle:
        records = list(SeqIO.parse(input_handle, "fasta"))

    # Use multiprocessing to parallelize the translation step
    pool = multiprocessing.Pool(threads)
    # Determine the chunk size based on the number of CPU threads
    chunk_size = len(records) // threads
    # Divide the records into chunks for parallel processing
    chunks = [records[i:i+chunk_size] for i in range(0, len(records), chunk_size)]
    # Process each chunk of records in parallel
    translated_records = pool.map(process_records, chunks)
    pool.close()
    pool.join()

    # Flatten the list of translated records
    translated_records = [record for chunk in translated_records for record in chunk]
    # Write the translated records to the output file
    write_records(translated_records, output_file)

# Get input and output file paths from Snakemake
input_file = snakemake.input[0]
output_file = snakemake.output[0]
# Get the number of threads from Snakemake
threads = snakemake.threads

# Process the file
process_file(input_file, output_file, threads)
