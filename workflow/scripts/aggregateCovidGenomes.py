import pandas as pd
from Bio import SeqIO


def getTopScoringAccessions(score_csvs, n):
    """Reads in all TaxID scores for the SARS-CoV-2 lineages and gets the top n from each.

    Args:
        score_csvs (list): List of file Paths to the TaxID scores of the lineages.

    Returns:
        list: Combined list of the n top scoring Taxon IDs per lineage.
    """
    all_top_accessions = []
    for csv in score_csvs:
        scores_df = pd.read_csv(csv, sep=",")
        scores_df = scores_df.sort_values(by="weight", ascending=False) # sort again, just to be sure
        top_rows = scores_df.head(n)
        top_accessions = top_rows["genbank_accession"].to_list()
        for accssion in top_accessions:
            all_top_accessions.append(accssion)

    return all_top_accessions


def parseAndWriteFasta(concat_fasta, filtered_fasta, all_top_accessions):
    """Reads in the concatenated fasta file and filters it for the entries returned by getTopScoringAccessions().

    Args:
        concat_fasta (string): Path to the concatenated fasta file.
        filtered_fasta (string): Path where the filtered concatenated fasta file should be stored.
        all_top_accessions (list): Combined list of the n top scoring Taxon IDs per lineage.
    """
    all_records = []
    for record in SeqIO.parse(concat_fasta, "fasta"):
        all_records.append(record)

    records_to_use = []

    for record in all_records:
        if record.id.split(".")[0] in all_top_accessions:
            records_to_use.append(record)
    
    SeqIO.write(records_to_use, filtered_fasta, "fasta")


def main():
    orf_scores_csvs = snakemake.input[:-1]
    concat_fasta = snakemake.input[-1]
    filtered_fasta = snakemake.output[0]
    n = snakemake.params[0]

    top_accessions = getTopScoringAccessions(orf_scores_csvs, n)
    parseAndWriteFasta(concat_fasta, filtered_fasta, top_accessions)

if __name__ == "__main__":
    main()
