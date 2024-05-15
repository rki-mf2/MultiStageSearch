import pandas as pd
from copy import deepcopy


def mapTaxIDs(in_file, taxids, out_file):
    """Reads in the PSM Report and the peptide taxid

    Args:
        in_file (string): Path to the PSM Report (PeptideShaker output of the first search step).
        taxids (string): Path to the tsv file that contains the mappings of peptide to taxon id.
        out_file (string): Path where the csv file is stored.

    Returns:
        pandas DataFrame: pandas DataFrame containing the matched peptides and their according weights, psmids and taxids.
    """

    in_file = pd.read_csv(in_file, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
    id_df = pd.read_csv(taxids, header=None, names=["protein", "taxid"], sep="\t")

    in_file.columns = ["accession"]
    in_file["weight"] = 1 / (in_file.accession.str.count(",") + 1)
    in_file["accession"] = in_file.accession.apply(lambda x: x.split(","))
    in_file["psmid"] = in_file.index + 1
    report_df = in_file.explode("accession", ignore_index=True)
    report_df_wo_version = deepcopy(report_df)
    report_df_wo_version["accession"] = report_df_wo_version["accession"].apply(lambda x: x.split(".")[0])

    merged_df = pd.merge(report_df, id_df, how="inner", left_on="accession", right_on="protein")
    merged_df_wo_version = pd.merge(report_df_wo_version, id_df, how="inner", left_on="accession", right_on="protein")
    concat_df = pd.concat([merged_df, merged_df_wo_version])
    concat_df = concat_df.drop("protein", axis=1)

    concat_df.to_csv(out_file, sep="\t", header=["accession", "weight", "psmid", "taxid"])
    return concat_df


def score(df, out_file):
    """Creates a csv file that contains the the taxon ids with the according weights

    Args:
        df (pandas DataFrame): DataFrame (output of mapTaxIDs).
        out_file (string): Path where the csv file is stored.
    """
    
    df_score = df.groupby('taxid')['weight'].sum().reset_index()
    df_score = df_score.sort_values(by=['weight'], ascending=False)
    df_score.to_csv(out_file, index=False, sep="\t")


def main():
    raw_report = snakemake.input[0]
    taxids = snakemake.input[1]
    mapped_taxids = snakemake.output[0]
    top_scoring = snakemake.output[1]

    merged_df = mapTaxIDs(in_file=raw_report, taxids=taxids, out_file=mapped_taxids)
    score(df=merged_df, out_file=top_scoring)


if __name__ == "__main__":
    main()
