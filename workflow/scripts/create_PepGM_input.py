import shutil
import pandas as pd


def createPepGMCsv(psm_report, strain_accessions, out_file):
    """Creates the csv file used as input in PepGM.

    Args:
        psm_report (string): path to the PSM report of the final search step
        strain_accessions (string): path to the 
        out_file (string): path to the csv file that is used as an input for PepGM
    """

    psm_report_df = pd.read_csv(psm_report, sep="\t")
    psm_report_df["Protein(s)"] =  psm_report_df["Protein(s)"].astype(str)
    psm_report_df["Protein(s)"] = psm_report_df["Protein(s)"].apply(lambda x: x.split(","))

    psm_report_df = psm_report_df.explode("Protein(s)", ignore_index=True)
    psm_report_df["Protein(s)"] = psm_report_df["Protein(s)"].apply(lambda x: x.split(".")[0])
    psm_report_df = psm_report_df[["Protein(s)", "Sequence", "Confidence [%]"]]
    psm_report_df = psm_report_df.rename(columns={"Confidence [%]": "score"})
    psm_report_df["score"] = psm_report_df["score"].replace(100.0, 99.9)
    psm_report_df["score"] = psm_report_df["score"] / 100

    strain_accessions_df = pd.read_csv(strain_accessions, sep=",")
    strain_accessions_df["genbank_accession"] = strain_accessions_df["genbank_accession"].apply(lambda x: x.split(".")[0])

    merged_df = psm_report_df.merge(strain_accessions_df[["genbank_accession", "taxa", "HigherTaxa"]], left_on="Protein(s)", right_on="genbank_accession", how="left")
    merged_df = merged_df.drop("genbank_accession", axis=1)
    merged_df = merged_df.rename(columns={"Sequence": "sequence"})
    # type casting
    merged_df["taxa"] = merged_df["taxa"].astype(int)
    merged_df["HigherTaxa"] = merged_df["HigherTaxa"].astype(int)
    merged_df.to_csv(out_file, sep=",", header=True, float_format='%.14f')


def createPepGMInputScores(orig_scores, out_file):
    """Creates the csv file used as input in PepGM.

    Args:
        orig_scores (string): path to the csv file containing the scores with the genbank accessions
        out_file (string): path to the csv file that is used as an input for PepGM
    """

    orig_scores_df = pd.read_csv(orig_scores, sep=",")
    pepGM_input_scores_df = orig_scores_df.drop("genbank_accession", axis=1)
    pepGM_input_scores_df = pepGM_input_scores_df[["taxid", "weight"]]
    pepGM_input_scores_df.to_csv(out_file, sep=",", header=True, index=False)

def copyStrainAccessions(strain_accessions, out_file):
    """Copies the strain accessions file to the output file.

    Args:
        strain_accessions (string): path to the csv file containing the strain accessions
        out_file (string): path to the csv file that is used as an input for PepGM
    """

    shutil.copyfile(strain_accessions, out_file)

def main():
    psm_report = snakemake.input[0]
    strain_accessions = snakemake.input[1]
    orig_scores = snakemake.input[2]
    pepgm_csv = snakemake.output[0]
    pepgm_input_scores = snakemake.output[1]
    pepgm_strain_accessions = snakemake.output[2]
    createPepGMCsv(psm_report, strain_accessions, pepgm_csv)
    createPepGMInputScores(orig_scores, pepgm_input_scores)
    copyStrainAccessions(strain_accessions, pepgm_strain_accessions)


if __name__ == "__main__":
    main()
