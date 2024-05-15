import pandas as pd

def convertReport(in_report, mapping_file, out_report):
    """Maps the ORF lists to the ORFs and replace them in the PeptideShaker report.

    Args:
        in_report (string): Path to the PeptideShaker report.
        mapping_file (string): Path to the file containing the ORF mapping.
        out_report (string): Path where the converted PeptideShaker report should be stored.
    """
    peptide_shaker_df = pd.read_csv(in_report, sep="\t")
    peptide_shaker_df["Protein(s)"] = peptide_shaker_df["Protein(s)"].apply(lambda x: x.split(","))
    peptide_shaker_df = peptide_shaker_df.explode("Protein(s)", ignore_index=True)
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=None, names=["ID_Mapping", "ORF_List"])
    id_mapping = mapping_df.set_index("ID_Mapping")["ORF_List"].to_dict()
    peptide_shaker_df["Protein(s)"] = peptide_shaker_df["Protein(s)"].map(id_mapping) # map lists to ORFs
    # aggregate the rows again
    peptide_shaker_df = peptide_shaker_df.groupby("Spectrum Title").agg({
            "Protein(s)": lambda x: ",".join(x),
            "Sequence": "first",
            "AAs Before": "first",
            "AAs After": "first",
            "Position": "first",
            "Modified Sequence": "first",
            "Variable Modifications": "first",
            "Fixed Modifications": "first",
            "Spectrum File": "first",
            "Spectrum Scan Number": "first",
            "RT": "first",
            "m/z": "first",
            "Measured Charge": "first",
            "Identification Charge": "first",
            "Theoretical Mass": "first",
            "Isotope Number": "first",
            "Precursor m/z Error [ppm]": "first",
            "Localization Confidence": "first",
            "Probabilistic PTM score": "first",
            "D-score": "first",
            "Confidence [%]": "first",
            "Validation": "first",
        }).reset_index()

    # sort columns
    peptide_shaker_df = peptide_shaker_df[["Protein(s)", "Sequence", "AAs Before", "AAs After", "Position", "Modified Sequence", "Variable Modifications", "Fixed Modifications", "Spectrum File", 
                                        "Spectrum Title", "Spectrum Scan Number", "RT", "m/z", "Measured Charge", "Identification Charge", "Theoretical Mass", "Isotope Number", "Precursor m/z Error [ppm]",
                                            "Localization Confidence", "Probabilistic PTM score", "D-score", "Confidence [%]", "Validation"]]


    peptide_shaker_df.to_csv(out_report, sep="\t")


def main():
    peptide_shaker_report = snakemake.input[0]
    mapping_file = snakemake.input[1]
    converted_report = snakemake.output[0]
    convertReport(peptide_shaker_report, mapping_file, converted_report)


if __name__ == "__main__":
    main()
