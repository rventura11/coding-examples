import subprocess

# Variables
reference_genome = "reference_genome.fasta"  # reference genome file path
reads_file = "reads.fastq"                   # FASTQ reads file path
output_file = "alignment.sam"                # output SAM file

# Functions - first one to index reference genome with bwa, second one to conduct alignment with mem or aln

def index_reference(reference_genome):
    command = ["bwa", "index", reference_genome]
    subprocess.run(command, check=True)
    print(f"Indexing completed for {reference_genome}")

def align_reads(reference_genome, reads_file, output_file, algorithm="mem"):
    if algorithm == "mem":
        command = ["bwa", "mem", reference_genome, reads_file]
    elif algorithm == "aln":
        # For aln your need a .sai file
        sai_file = reads_file + ".sai"
        command_aln = ["bwa", "aln", reference_genome, reads_file]
        subprocess.run(command_aln, check=True)
        command = ["bwa", "samse", reference_genome, sai_file, reads_file]
    else:
        raise ValueError("Unsupported algorithm. Use 'mem' for long reads or 'aln' for short reads.")

    with open(output_file, "w") as out:
        subprocess.run(command, stdout=out, check=True)
    
    print(f"Alignment completed. Output saved as {output_file}")


# Step 1: Index reference genome
index_reference(reference_genome)

# Step 2: Align reads
align_reads(reference_genome, reads_file, output_file, algorithm="mem") # or algorithm="aln" for short reads
