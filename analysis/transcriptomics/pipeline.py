import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess

def load_fasta(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(record)
    return sequences

def remove_seqrecord(sequences, ids_to_remove):
    sequences_to_keep = []
    for sequence in sequences:
        if sequence.id not in [id_to_remove.id for id_to_remove in ids_to_remove]:
            sequences_to_keep.append(sequence)
    return sequences_to_keep

def split_sequence(seqrecord, target_sequence):
    target_sequence = str(target_sequence)
    seq = str(seqrecord.seq)
    
    start = seq.find(target_sequence)
    end = start + len(target_sequence)
    
    if start == -1:
        target_sequence = str(Seq(target_sequence).reverse_complement())
        start = seq.find(target_sequence)
        end = start + len(target_sequence)
    
    if start == -1:
        raise ValueError(f"Target sequence {target_sequence} not found in {seqrecord.id}")
    
    return seq[:start], target_sequence, seq[end:]

def create_seqrecords(split_sequence, prefix):
    seq1, seq2, seq3 = split_sequence
    id1, id2, id3 = f"{prefix}_pre_split", f"{prefix}_target", f"{prefix}_post_split"
    seq1 = SeqRecord(seq=Seq(seq1), id=id1)
    seq2 = SeqRecord(seq=Seq(seq2), id=id2)
    seq3 = SeqRecord(seq=Seq(seq3), id=id3)
    return seq1, seq2, seq3

def create_kallisto_index(genome_fna):
    index_prefix = os.path.splitext(genome_fna)[0]
    if os.path.exists(f"{index_prefix}.idx"):
        print(f"Index file {index_prefix}.idx already exists, skipping index creation")
        return
    command = f"kallisto index -i {index_prefix}.idx {genome_fna}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def run_kallisto_quant(index_prefix, read1, read2, target_sequence):
    outdir = 'IMBMDB_Kallisto_output'
    # Make the main output directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Include a sanitized version of the target sequence in the output directory name
    safe_target_sequence = target_sequence.replace("/", "_").replace(" ", "_")[:50]  # Truncate for safety
    index_one_layer = index_prefix.split("/")[-1]
    output_dir = f"{outdir}/kallisto_output_{index_one_layer}_{safe_target_sequence}"
    i = 1
    while os.path.exists(output_dir):
        output_dir = f"{outdir}/kallisto_output_{index_one_layer}_{safe_target_sequence}_{i}"
        i += 1

    os.makedirs(output_dir, exist_ok=True)
    
    command = f"kallisto quant -i {index_prefix}.idx -o {output_dir} --threads 20 {read1} {read2}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def update_genome_file(genome_fna, new_contigs, original_contig_id):
    """
    Updates the genome file by replacing the original contig with new contigs.
    """
    # Load all sequences from the original genome
    all_sequences = load_fasta(genome_fna)
    
    # Filter out the original contig
    filtered_sequences = remove_seqrecord(all_sequences, [SeqRecord(Seq(""), id=original_contig_id)])
    
    # Convert new_contigs to a list and add them to the filtered sequences
    updated_sequences = filtered_sequences + list(new_contigs)
    
    # Write the updated sequences back to the genome file
    with open(genome_fna, "w") as output_handle:
        SeqIO.write(updated_sequences, output_handle, "fasta")

def main_pipeline(genome_fna, read1, read2, target_sequence):
    sequences = load_fasta(genome_fna)
    original_contig_id = None
    new_contigs = None
    
    for seqrecord in sequences:
        try:
            pre, target, post = split_sequence(seqrecord, target_sequence)
            print(f"Target sequence found and split in {seqrecord.id}")
            original_contig_id = seqrecord.id
            new_contigs = create_seqrecords((pre, target, post), original_contig_id)
            break
        except ValueError:
            continue
    else:
        raise ValueError(f"Target sequence {target_sequence} not found in any sequences in {genome_fna}")
    
    # Update the genome file with the new contigs
    update_genome_file(genome_fna, new_contigs, original_contig_id)
    
    # Create Kallisto index and run quantification
    create_kallisto_index(genome_fna)
    run_kallisto_quant(os.path.splitext(genome_fna)[0], read1, read2, target_sequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Pipeline for genome modification and Kallisto analysis')
    parser.add_argument('--genome_fna', required=True, help='Path to genome FASTA file')
    parser.add_argument('--read1', required=True, help='Path to R1 FASTQ file')
    parser.add_argument('--read2', required=True, help='Path to R2 FASTQ file')
    parser.add_argument('--target_sequence', required=True, help='Target sequence to split and analyze')
    args = parser.parse_args()
    
    main_pipeline(args.genome_fna, args.read1, args.read2, args.target_sequence)
