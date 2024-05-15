import logging
import multiprocessing
import os
import pandas as pd
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def readFasta(fasta_file):
    """Reads a fasta file and returns the sequence records.

    Args:
        fasta_file (string): Path to the fasta file.

    Returns:
        list: List of all sequence records from the fasta file.
    """
    all_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        all_records.append(record)
    return all_records


def getORFsToInspect(all_records, orfs):
    """Gets all ORFs that correspond to the top scoring genomes.

    Args:
        all_records (list): List of sequence records.
        orfs (list): List of ORF IDs.

    Returns:
        list: List of lists of ORFs.
    """
    orfs_to_inspect = {}
    for accession in orfs:
        accession_orfs = []
        for record in all_records:
            if record.id.split(".")[0] == accession:
                accession_orfs.append(record)
        orfs_to_inspect[accession] = accession_orfs
    return orfs_to_inspect


def createHeatmap(similarity_matrix, orfs, output):
    """Creates the heatmap of the pairwise similarities.

    Args:
        similarity_matrix (numpy array): Matrix containing the pairwise similarities.
        orfs (list): List of genbank accessions.
        output (string): Path where the matrix png will be stored.
    """

    plt.figure(figsize=(len(similarity_matrix), round((len(similarity_matrix) / 5) * 4))) # 5:4 ratio
    sns.heatmap(similarity_matrix, annot=True, cmap='coolwarm', fmt=".2f")

    plt.title('Pairwise Similarities Heatmap')
    plt.xlabel('Accessions')
    plt.ylabel('Accessions')

    tick_labels = orfs 
    plt.xticks(np.arange(len(tick_labels)) + 0.5, tick_labels, rotation=45)
    plt.yticks(np.arange(len(tick_labels)) + 0.5, tick_labels, rotation=0)

    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


def compute_task(args):
    """Calls the computing function (used for multithreading).

    Args:
        args (tuple): tuple containing the accessions and ORFs to compare.

    Returns:
        int: index of the first accession.
        int: index of the second accession.
        float: similarity between the two proteomes.
    """
    accession1, accession2, orfs1, orfs2 = args
    pid = os.getpid()
    logging.info(f"Process ID: {pid}, computing similarity for indices ({accession1}, {accession2}), {orfs1[0].id}, {orfs2[0].id}")
    similarity = compute_proteome_similarity(orfs1, orfs2)
    return accession1, accession2, similarity


def compute_proteome_similarity(orfs1, orfs2):
    """Computes the pairwise similarity between two proteomes.

    Args:
        orfs1 (list): List of the ORFs of the first proteome.
        orfs2 (list): List of the ORFs of the second proteome.

    Returns:
        float: Pairwise similarity between the two proteomes.
    """
    sorted_orfs1 = sorted([orf.seq for orf in orfs1])
    sorted_orfs2 = sorted([orf.seq for orf in orfs2])
    count = 0
    index = 0
    # count identical ORFs, each ORF can only be counted once as a duplicate
    for orf1 in range(len(sorted_orfs1)):
        for orf2 in range(len(sorted_orfs2))[index:]:
            if sorted_orfs1[orf1] == sorted_orfs2[orf2]:
                count += 1
                index = orf2 + 1
                break
            elif sorted_orfs1[orf1] < sorted_orfs2[orf2]:
                index = orf2
                break

    similarity = count /  max(len(orfs1), len(orfs2))
    return similarity


def compute_similarity_matrix(orfs_to_inspect, orfs, num_threads):
    """Computes the pairwise similarity matrix of the proteomes (Multithreaded).

    Args:
        orfs_to_inspect (list): List of ORFs to inspect.
        orfs (list): List of genbank accessions.
        num_threads (int): Number of available threads.

    Returns:
        numpy array: Matrix containing the pairwise similarities.
    """
    tasks = [(accession, accession2, orfs_to_inspect[orfs[accession]], orfs_to_inspect[orfs[accession2]]) 
             for accession in range(len(orfs)) for accession2 in range(accession+1, len(orfs))]
    similarity_matrix = np.zeros((len(orfs), len(orfs)))

    with multiprocessing.Pool(processes=num_threads) as pool:
        results = pool.map(compute_task, tasks)

    for accession1, accession2, similarity in results:
        similarity_matrix[accession1][accession2] = similarity
        similarity_matrix[accession2][accession1] = similarity
    
    np.fill_diagonal(similarity_matrix, 1)

    return similarity_matrix


def main():
    orfs_csv = snakemake.input[0]
    proteome_fasta = snakemake.input[1]
    num_threads = snakemake.params[0]
    heatmap_png = snakemake.output[0]
    std_out = snakemake.log[1]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    orfs_df = pd.read_csv(orfs_csv)
    orfs = orfs_df["genbank_accession"].to_list()

    all_records = readFasta(proteome_fasta)
    orfs_to_inspect =  getORFsToInspect(all_records, orfs)
    similarity_matrix = compute_similarity_matrix(orfs_to_inspect, orfs, num_threads)   
    createHeatmap(similarity_matrix, orfs, heatmap_png)


if __name__ == '__main__':
    main()