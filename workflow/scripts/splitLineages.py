import os
import logging
import pandas as pd
from Bio import SeqIO


def processFetchedData(fetched_data, n):
    """Reads in the strain_accessions.csv and splits the lineages into categories.
    Each lineage with at least n entries in it, will be separately searched. All lineages
    with less than n will be categorized as "others". They are search together. 
    Each lineage is splitted at the dot. That means that "B" and "B.1" are both 
    considered as "B" for the categorization.

    Args:
        fetched_data (string): Path to the strain_accessions.csv
        n (integer): Threshold of minunum entries for one lineage to not be 
            categorized as "others"

    Returns:
        pandas DataFrame: pandas DataFrame containing the categorization of the lineages.
    """
    fetched_df = pd.read_csv(fetched_data, sep=",")
    fetched_df["major_lineage"] = fetched_df["strain"].apply(lambda x: x.split(".")[0])

    lineages_counts = fetched_df["major_lineage"].value_counts()
    abundant_lineages = [lineages_counts.index[lineage] for lineage in range(len(lineages_counts)) if lineages_counts[lineage] >= n] # get all lineages with at least n entries

    def classifyLineage(row):
        return row['major_lineage'] if row['major_lineage'] in abundant_lineages else "others" # categorize every entry

    fetched_df['classification'] = fetched_df.apply(classifyLineage, axis=1)

    return fetched_df


def splitLineages(processed_df, fasta, dir_path):
    """Splits the lineages based on the categorization. For each,
      a subdirectory is created with a concatenated fasta file 
      containing the genomes for this lineage.

    Args:
        processed_df (pandas DataFrame): pandas DataFrame containing the categorization of the lineages.
        fasta (string): Path to the concatenated fasta file containing all downloaded genomes.
        dir_path (string): Path where the subfolder will be created.
    """
    all_records = list(SeqIO.parse(fasta, "fasta"))
    for classifiction in processed_df["classification"].unique():
        logging.info(classifiction)
        classifiction_df = processed_df[processed_df["classification"] == classifiction]
        accessions = classifiction_df["genbank_accession"].tolist()
        records = [record for record in all_records if record.id in accessions]
        try:
            os.makedirs(f"{dir_path}/{classifiction}", exist_ok=True)
            logging.info(f"Directory '{dir_path}/{classifiction}' was created.")
            SeqIO.write(records, f"{dir_path}/{classifiction}/{classifiction}.fasta", "fasta")
        except OSError as error:
            logging.info(f"Creation of the directory '{dir_path}/{classifiction}' failed due to: {error}")
    

def main():
    fetched_data = snakemake.input[0]
    fasta = snakemake.input[1]
    dir_path = snakemake.params[0]
    n = snakemake.params[1]
    std_out = snakemake.log[0]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    processed_df = processFetchedData(fetched_data, n)
    splitLineages(processed_df, fasta, dir_path)

if __name__ == "__main__":
    main()