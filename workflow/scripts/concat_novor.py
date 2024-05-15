import re
import glob 
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def removeModificationNotation(peptide):
    """Removes the modifications of the sequence.

    Args:
        peptide (string): Sequence containing the modifications.

    Returns:
        string: Sequence without modifications.
    """
    aminoacid_sequence = re.sub(r"\(.*?\)|N-term\|\d+", "", peptide)
    return aminoacid_sequence


def concatNovorRecords(novor_csv, fasta, concat_fasta):
    """Parses the NOVOR csv for the peptides, giving them unique identifier 
       and concatenes the sequence database and the NOVOR records.

    Args:
        novor_csv (string): Path to the csv file created by NOVOR.
        fasta (string): Path to the fasta file containing the database sequences.
        concat_fasta (string): Path to the concatenated fasta file.
    """
    all_records = []
    for record in SeqIO.parse(fasta, "fasta"):
        all_records.append(record)

    novor_df = pd.read_csv(novor_csv, sep=", ", comment="#", header=None) # ignore comments in csv
    novor_df[9] = novor_df[9].apply(removeModificationNotation) # novor_df[9] = peptide sequence
    for index, row in novor_df.iterrows():
        record_id = "NOVOR_PEPTIDE_" + str(row[0]) # row[0] = id
        record_sequence = Seq(row[9])
        record = SeqRecord(record_sequence, id=record_id, description="")
        all_records.append(record)

    SeqIO.write(all_records, concat_fasta, "fasta")


def main():
    novor_dir = snakemake.params[0]
    fasta = snakemake.input[1]
    concat_fasta = snakemake.output[0]

    novor_csv = glob.glob(f"{novor_dir}/*.novor.csv", recursive=True)[0]
    concatNovorRecords(novor_csv, fasta, concat_fasta)


if __name__ == "__main__":
    main()
