import logging
import pandas as pd
from Bio import SeqIO

def filterUsedGenomes(fasta, orf_scores):
    """Gets the genomes that are used for further steps (i.e. Phylogeny).

    Args:
        fasta (string): Path to the fasta file containing the genomes.
        orf_scores (string): Path to the file containing the ORF scores.

    Returns:
        list: List of the top scoring genomes.
    """
    orf_scores = pd.read_csv(orf_scores, sep=",", header=None, names=["genbang_accession", "orf", "strain_isolate", "score"])
    genbank_accessions = orf_scores["genbang_accession"].tolist()

    used_genomes = []
    num_records = 0
    for record in SeqIO.parse(fasta, "fasta"):
        num_records += 1
        genome_id = record.id.split(" ")[0]
        genome_id = genome_id.split(".")[0]
        if genome_id in genbank_accessions:
            used_genomes.append(record)
    logging.info(f"{len(used_genomes)} genomes out of {num_records} fetched genomes are used for the phylogeny.")
    return used_genomes

def writeFilteredGenomes(used_genomes, out_fasta):
    """Writes the top scoring genomes to a fasta file.

    Args:
        used_genomes (list): List of the top scoring genomes.
        out_fasta (string): Path where the fasta file is stored
    """
    SeqIO.write(used_genomes, out_fasta, "fasta")

def main():
    genomes_fasta = snakemake.input[0]
    orf_scores = snakemake.input[1]
    used_genomes_fasta = snakemake.output[0]
    std_out = snakemake.log[0]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    used_genomes = filterUsedGenomes(genomes_fasta, orf_scores)
    writeFilteredGenomes(used_genomes, used_genomes_fasta)

if __name__ == "__main__":
    main()