import logging
import pandas as pd
import numpy as np
from Bio import SeqIO


def weightORFs(peptide_shaker_report):
    """Reads in the PSM Report and the weights the ORFs.

    Args:
        peptide_shaker_report (string): Path where the PSM Report (Output of PeptideShaker from the final search steps) is stored

    Returns:
        pandas DataFrame: DataFrame containing the Genbank accessions and their according psmid, weight and confidence of the PSM
    """
    
    df = pd.read_csv(peptide_shaker_report, sep="\t", on_bad_lines="skip", usecols=["Protein(s)", "Confidence [%]"])
    df.columns = ["ORF", "Confidence"]
    df["ORF"] = df["ORF"].astype(str)
    df["weight"] = 1 / (df.ORF.str.count(",") + 1)
    df["ORF"] = df.ORF.apply(lambda x: x.split(","))
    df["psmid"] = df.index + 1
    report_df = df.explode("ORF", ignore_index=True)
    report_df["accession"] = report_df.ORF.apply(lambda x: x.split(".")[0])
    accessions = report_df["accession"].unique().tolist()

    return report_df, accessions


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


def getORFsToInspect(all_records, accessions):
    """Gets all ORFs that correspond to the top scoring genomes.

    Args:
        all_records (list): List of sequence records.
        orfs (list): List of ORF IDs.

    Returns:
        dict: dict of lists of ORFs.
    """
    orfs_to_inspect = {}
    for accession in accessions:
        accession_orfs = []
        for record in all_records:
            if record.id.split(".")[0] == accession:
                accession_orfs.append(record)
        orfs_to_inspect[accession] = accession_orfs
    return orfs_to_inspect


def createGenomeStats(report_df):
    """Counts the PSMs of each genome and compute the mean confidence.

    Args:
        report_df (pandas DataFrame): output of weightORFs()

    Returns:
        pandas DataFrame: DataFrame containing different statistics of the genome
    """

    genome_stats = pd.DataFrame({
        "accession": report_df["accession"],
        "counts": report_df.groupby("accession")["accession"].transform("count")
    })

    mean_conf_df = report_df.groupby("accession")["Confidence"].mean().reset_index()
    mean_conf_df.rename(columns={"Confidence": "mean_confidence"}, inplace=True)
    genome_stats = genome_stats.merge(mean_conf_df, on="accession")
    genome_stats.drop_duplicates(inplace=True)
    genome_stats.reset_index(drop=True, inplace=True)

    return genome_stats


def getStrainAccessionsDF(strain_accessions):
    """Reads in a csv file containing the  genbank_accession, species, strain, isolate, taxa and HigherTaxa
       and converts the ORF ID to get the real genbank accession.

    Args:
        strain_accessions (string): Path where the csv file is stored.

    Returns:
        pandas DataFrame: DataFrame where the ORF ID is converted to the real genbank accession.
    """

    accessions_df = pd.read_csv(strain_accessions, sep=",")
    accessions_df["genbank_accession"] = accessions_df.genbank_accession.apply(lambda x: x.split(".")[0])

    return accessions_df


def mergeDataFrames(accessions_df, genome_stats, relevant_genomes, orfs_dict, out_file):
    """Merges the DataFrames and also computes a confidence scoring.
       This function also writes the merged DataFrame into a tsv file.

    Args:
        accessions_df (pandas DataFrame): Output of getStrainAccessionsDF()
        genome_stats (pandas DataFrame): Output of createGenomeStats()
        relevant_genomes (integer): Maximum number of genomes considered for the merged DataFrame
        orfs_dict (dict): Dict of lists of ORFs.
        out_file (string): Path where the tsv file is stored.
    """
    
    merged_df = pd.merge(accessions_df, genome_stats, how="inner", left_on="genbank_accession", right_on="accession")
    merged_df.replace(np.nan, "NaN", regex=True, inplace=True)
    merged_df.sort_values("counts", ascending=False, inplace=True)
    merged_df.reset_index(drop=True, inplace=True)
    merged_df.drop(columns=["accession"])
    merged_df = merged_df[:relevant_genomes]
    merged_df["confidence_scoring"] = merged_df["counts"] * merged_df["mean_confidence"]
    merged_df["length_scoring"] = merged_df.apply(lambda row: row["counts"] / len(orfs_dict[row["genbank_accession"]]), axis=1)

    merged_df.to_csv(out_file, index=False, sep="\t")


def main():
    peptide_shaker_report = snakemake.input[0]
    strain_accessions = snakemake.input[1]
    fasta = snakemake.input[2]
    out_file = snakemake.output[0]
    relevant_genomes = snakemake.params[0]

    std_out = snakemake.log[0]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    report_df, accessions = weightORFs(peptide_shaker_report)
    logging.info(f"weightORFs done")
    all_records = readFasta(fasta)
    logging.info(f"readFasta done")
    orfs_dict = getORFsToInspect(all_records, accessions)
    logging.info(f"getORFsToInspect done")
    genome_stats_df = createGenomeStats(report_df)
    logging.info(f"createGenomeStats done")
    accessions_df = getStrainAccessionsDF(strain_accessions)
    logging.info(f"getStrainAccessionsDF done")
    mergeDataFrames(accessions_df, genome_stats_df, relevant_genomes, orfs_dict, out_file)
    logging.info(f"mergeDataFrames done")

if __name__ == "__main__":
    main()
