import logging
import pandas as pd
from Bio import SeqIO
import re


def getAccessionsWithDuplicateNames(accessions_csv, words_blacklist):
    """Gets all genbank accessions with duplicate strain names 
        and all genbank accessions with blacklisted words in the name.

    Args:
        accessions_csv (string): Path to the csv containing the genbank accessions.
        words_blacklist (list): List of blacklisted words. 
                                They are removed from the strain_isolate names.

    Returns:
        list: List of genbank accessions with duplicate strain names.
        list: List of genbank accessions with blacklisted words.
    """
    accessions_df = pd.read_csv(accessions_csv)
    accessions_df["strain"] = "strain " + accessions_df["strain"].astype(str)
    accessions_df["isolate"] = "isolate " + accessions_df["isolate"].astype(str)

    accessions_df["strain_isolate"] = accessions_df["species"] + " "
    accessions_df.loc[accessions_df["strain"] != "strain nan", "strain_isolate"] += accessions_df["strain"] + " "
    accessions_df.loc[accessions_df["isolate"] != "isolate nan", "strain_isolate"] += accessions_df["isolate"] + " "
    accessions_df["strain_isolate"] = accessions_df["strain_isolate"].str.lower() # lower string for comparison
    entries_with_blacklisted_words = []
    for blacklisted_word in words_blacklist:
        regex = r"\b{}\b".format(re.escape(blacklisted_word.lower())) # "\b" = word boundary
        current_entries_with_blacklisted_words = accessions_df[accessions_df['strain_isolate'].str.contains(regex, regex=True, na=False) == True]["genbank_accession"].to_list() # get all genbank accessions that contains a blacklisted word
        for entry in current_entries_with_blacklisted_words:
            entries_with_blacklisted_words.append(entry) # collect all genbank accessions with blacklisted words 

        accessions_df["strain_isolate"] = accessions_df["strain_isolate"].str.replace(regex, "", regex=True) # remove blacklisted word
        accessions_df["strain_isolate"] = accessions_df["strain_isolate"].str.replace("  ", " ", regex=True) # if the word was in the middle of the string there could be to spaces, replace them by one

    entries_with_blacklisted_words = list(set(entries_with_blacklisted_words)) # if a strain contains multiple blacklisted words it appears multiple times in the list, drop all duplicates
    accessions_df["strain_isolate"] = accessions_df["strain_isolate"].str.strip(" ")


    duplicates = accessions_df[accessions_df["strain_isolate"].duplicated(keep=False)]
    
    duplicate_accessions = []
    for duplicate in duplicates["strain_isolate"].unique():
        current_duplicate_accessions = duplicates.loc[duplicates["strain_isolate"] == duplicate]["genbank_accession"].to_list()
        duplicate_accessions.append(current_duplicate_accessions)

    return duplicate_accessions, entries_with_blacklisted_words


def readFasta(fasta_file):
    """Reads a fasta file and returns all sequence records.

    Args:
        fasta_file (string): Path to the fasta file.

    Returns:
        list: List of all sequence records in the fasta file.
    """
    all_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        all_records.append(record)
    return all_records


def getORFsToInspect(duplicate_accessions, all_records):
    """Gets all ORFs for the duplicate accessions.

    Args:
        duplicate_accessions (list): List of genbank accessions with duplicate strain names.
        all_records (list): List of sequence records.

    Returns:
        dict: Dictionary with the number of ORFs for each accession.
        dict: Dictionary with all ORFs for each accession.
    """
    num_accessions = {}
    orfs_to_inspect = {}
    for duplicate_strains in duplicate_accessions:
        for accession in duplicate_strains:
            accession_orfs = []
            count = 0
            for record in all_records:
                record_id_until_point = record.id.split(".")[0]
                underscore_count = record_id_until_point.count("_")
                if underscore_count > 0:
                    split_position = record.id.find("_", len(record_id_until_point) + 1)
                    if record.id[0:split_position] == accession:
                        count += 1
                        accession_orfs.append(record)
                elif record.id.split("_")[0] == accession:
                    count += 1
                    accession_orfs.append(record)
            num_accessions[accession] = count
            orfs_to_inspect[accession] = accession_orfs
    return num_accessions, orfs_to_inspect


def filterDuplicateProteomes(duplicate_accessions, 
                             entries_with_blacklisted_words, 
                             num_accessions, 
                             orfs_to_inspect, 
                             similarity_threshold, 
                             all_records):
    """_summary_

    Args:
        duplicate_accessions (list): List of genbank accessions with duplicate strain names.
        entries_with_blacklisted_words (list): List of genbank accessions with blacklisted words.
        num_accessions (dict): Dictionary containing the number of ORFs for each accession.
        orfs_to_inspect (dict): Dictionary containing all ORFs for each accession.
        similarity_threshold (float): Threshhold for the pairwise similarity of Proteomes.
        all_records (list): List of sequence records.

    Returns:
        list: List of sequence records without the duplicate proteomes.
    """
    for duplicate_strains in duplicate_accessions:
        for accession in duplicate_strains:
            orfs_of_accession1 = sorted([orf.seq for orf in orfs_to_inspect[accession]])
            for accession2 in duplicate_strains[duplicate_strains.index(accession) + 1:]:
                orfs_of_accession2 = sorted([orf.seq for orf in orfs_to_inspect[accession2]])
                count = 0
                index = 0
                # count identical ORFs, each ORF can only be counted once as a duplicate
                for orf1 in range(len(orfs_of_accession1)):
                    for orf2 in range(len(orfs_of_accession2))[index:]:
                        if orfs_of_accession1[orf1] == orfs_of_accession2[orf2]:
                            count += 1
                            index = orf2 + 1
                            break
                        elif orfs_of_accession1[orf1] < orfs_of_accession2[orf2]:
                            index = orf2
                            break
                try:
                    similarity = count / max(num_accessions[accession], num_accessions[accession2])
                except ZeroDivisionError:
                    logging.info(f"There was a ZeroDivisionError. The Similarity of the ORFS  of {accession} and {accession2} is set to 0!")
                    similarity = 0.0
                logging.info(f"Similarity of ORFs of {accession} and {accession2} = {similarity}")

                old_len = len(all_records)
                if similarity > similarity_threshold:
                    logging.info(f"Similarity exceeds threshold of {similarity_threshold * 100}%!")
                    if accession in entries_with_blacklisted_words:
                        logging.info(f"Removing ORFs for Genbank accession {accession}!")
                        orf_ids_to_remove = [orf.id for orf in orfs_to_inspect[accession]]
                    elif accession2 in entries_with_blacklisted_words:
                        logging.info(f"Removing ORFs for Genbank accession {accession2}!")
                        orf_ids_to_remove = [orf.id for orf in orfs_to_inspect[accession2]]
                    elif num_accessions[accession] < num_accessions[accession2]:
                        logging.info(f"Removing ORFs for Genbank accession {accession}!")
                        orf_ids_to_remove = [orf.id for orf in orfs_to_inspect[accession]]
                    else:
                        logging.info(f"Removing ORFs for Genbank accession {accession2}!")
                        orf_ids_to_remove = [orf.id for orf in orfs_to_inspect[accession2]]
                    all_records = [record for record in all_records if record.id not in orf_ids_to_remove]
                    logging.info(f"Removed {old_len - len(all_records)} ORFs.\n")
                else:
                    logging.info(f"Similarity threshold is not met. Keeping ORFs!\n")  
    return all_records


def writeFilteredFasta(filtered_records, filtered_output_fasta):
    """Writes the filtered records to a fasta file.

    Args:
        filtered_records (list): List of sequence records without the duplicate proteomes.
        filtered_output_fasta (string): Path where the filtered fasta is stored.
    """
    with open(filtered_output_fasta, "w") as f:
        SeqIO.write(filtered_records, f, "fasta")

def main():
    accessions_csv = snakemake.input[0]
    proteomes_fasta = snakemake.input[1]
    similariy_threshold = snakemake.params[0]
    words_blacklist = snakemake.params[1]
    filtered_proteomes_fasta = snakemake.output[0]
    std_out = snakemake.log[1]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    duplicate_accessions, entries_with_blacklisted_words = getAccessionsWithDuplicateNames(accessions_csv, words_blacklist)
    all_records = readFasta(proteomes_fasta)
    num_accessions, orfs_to_inspect = getORFsToInspect(duplicate_accessions, all_records)
    filtered_records = filterDuplicateProteomes(duplicate_accessions, entries_with_blacklisted_words, num_accessions, orfs_to_inspect, similariy_threshold, all_records)
    writeFilteredFasta(filtered_records, filtered_proteomes_fasta)



if __name__ == "__main__":
    main()