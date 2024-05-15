import pandas as pd
from Bio import SeqIO

def getAccessions(in_file):
    """Truncates the ORF identifier to the genbank ID.

    Args:
        in_file (string): Path to the fasta file containing the ORFs.
    
    Returns:
        pandas DataFrame: pandas DataFrame containing the mapping of ORFs and genbank IDs.

    """
    accessions_dict = {"orf": [], "genbank_accession": []}
    records= SeqIO.parse(in_file, "fasta")
    for record in records:
        orfs = (str(record.id)).split(",")
        for orf in orfs:
            accession = orf.split(".")[0]
            accessions_dict["orf"].append(orf)
            accessions_dict["genbank_accession"].append(accession)
            
    accessions_df = pd.DataFrame(data=accessions_dict)
    return accessions_df

def mapORFs(in_file, accessions_df, out_file):
    """Reads in the final PSM Report and the ORF-Genbank ID mapping and weights the ORFs.

    Args:
        in_file (string): Path to the PSM Report (PeptideShaker output of the final search step).
        accessions_df (pandas DataFrame): pandas DataFrame containing the mapping of ORFs and genbank IDs.
        out_file (string): Path where the csv file is stored.

    Returns:
        pandas DataFrame: pandas DataFrame containing the matched ORFs and their according weights, psmids and genbank accesions.
    """

    in_file = pd.read_csv(in_file, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])

    in_file.columns = ["accession"]
    in_file["accession"] = in_file["accession"].astype(str)
    in_file["weight"] = 1 / (in_file.accession.str.count(",") + 1)
    in_file["accession"] = in_file.accession.apply(lambda x: x.split(","))
    in_file["psmid"] = in_file.index + 1
    report_df = in_file.explode("accession", ignore_index=True)

    merged_df = pd.merge(report_df, accessions_df, how="inner", left_on="accession", right_on="orf")
    merged_df = merged_df.drop("orf", axis=1)

    merged_df.to_csv(out_file, sep="\t", header=["accession", "weight", "psmid", "genbank_accession"])
    return merged_df


def scoreGenbankAccessions(df):
    """Create the scores for the genbank accessions.

    Args:
        df (pandas DataFrame): pandas DataFrame containing the weights for the ORFs

    Returns:
        pandas DataFrame: pandas DataFrame containing the genbank accessions with their according weights.
    """

    df_score = df.groupby('genbank_accession')['weight'].sum().reset_index()
    df_score = df_score.sort_values(by=['weight'], ascending=False)
    return df_score


def scoreTaxIDs(df_score, strain_taxids, out_file):
    """Maps the genbank accessions to taxon IDs and creates the scoring for taxon IDs.

    Args:
        df_score (pandas DataFrame): pandas DataFrame containing the genbank accessions with their according weights.
        strain_taxids (string): Path to the tsv file containing among others the mapping of genbank accession to taxon ID.
        out_file (string): Path where the resulting csv file containing the scoring for the taxon IDs should be stored.
    """
    strain_accessions = pd.read_csv(strain_taxids, sep = "\t")
    strain_accessions = strain_accessions.drop_duplicates(subset= ["species", "strain", "isolate", "taxa", "HigherTaxa"], keep="first").reset_index(drop=True)

    strain_accessions["strain"] = strain_accessions["strain"].astype(str)
    strain_accessions["isolate"] = strain_accessions["isolate"].astype(str)
    strain_accessions["strain_isolate"] = strain_accessions["strain"].replace("nan", "") + "_" + strain_accessions["isolate"].replace("nan", "")
    strain_accessions["strain_isolate"] = strain_accessions["strain_isolate"].str.strip("_")
    strain_accessions['strain_isolate'] = strain_accessions["strain_isolate"].replace("", "NaN")
    
    taxid_weights_df = df_score.merge(strain_accessions[["genbank_accession", "strain_isolate", "taxa"]], on="genbank_accession", how="right")
    taxid_weights_df = taxid_weights_df.rename(columns={"taxa": "taxid"})
    taxid_weights_df = taxid_weights_df.sort_values(by="weight", ascending=False)
    taxid_weights_df = taxid_weights_df[["genbank_accession", "taxid", "strain_isolate", "weight"]]
    taxid_weights_df = taxid_weights_df.reset_index(drop=True)
    taxid_weights_df["weight"] = taxid_weights_df["weight"].astype(float)
    taxid_weights_df = taxid_weights_df.drop_duplicates(subset= ["taxid"], keep="first").reset_index(drop=True)
    taxid_weights_df.to_csv(out_file, sep=",", index=False)


def main():
    fasta_file = snakemake.input[0]
    raw_report = snakemake.input[1]
    strain_taxids = snakemake.input[2]
    mapped_OFRs = snakemake.output[0]
    top_scoring = snakemake.output[1]

    accessions_df = getAccessions(fasta_file)
    merged_df = mapORFs(in_file=raw_report, accessions_df=accessions_df, out_file=mapped_OFRs)
    df_score = scoreGenbankAccessions(df=merged_df)
    scoreTaxIDs(df_score=df_score, strain_taxids=strain_taxids, out_file=top_scoring)


if __name__ == "__main__":
    main()
