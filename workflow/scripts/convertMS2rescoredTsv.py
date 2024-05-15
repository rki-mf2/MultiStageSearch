import pandas as pd
import xml.etree.ElementTree as ET
import ast
import re
import logging

pd.set_option('display.precision', 20)


def removeModsFromSequence(sequence):
    """Removes the modifications of the sequence.

    Args:
        sequence (string): Sequence containing the modifications.

    Returns:
        string: Sequence without modifications.
    """
    return re.sub(r"\[.*?\]|/.*", "", sequence)


def extractIdentifier(seq_id):
    """Extracts the protein ID from a sequence identifier, which could contain other symbols.

    Args:
        seq_id (string): Sequence identifier of a fasta sequence containing the protein ID.

    Returns:
        string: extracted protein ID.
    """
    # split Cont entries
    if "|" in seq_id:
        id_list = seq_id.split("|")
        for entry in id_list:
            if entry.startswith("Cont_"):
                return entry
        return seq_id
    # split every other entries
    else:
        id_list = seq_id.split(" ")
        return id_list[0]


def parseXTandemXML(tandem_xml):
    """Parses the XTandem XML and returns the labels of the entries.

    Args:
        tandem_xml (string): Path to the XTandem results XML.

    Returns:
        dict: Dictionary containing the labels of the entries. 
            This corresponds to the spectrum title of the MGF file.
    """
    tree = ET.parse(tandem_xml)
    root = tree.getroot()

    label_dict = {}
    for group in root.findall("group"):
        try:
            group_id = group.get("id")
            inner_group = group.find("group")
            note = inner_group.find("note")
            label = note.text
            label_split = label.split(" RTINSECONDS")[0]
            label_dict[group_id]= label_split
        except:
            continue
    
    return label_dict

def convert_tsv(ms2rescore_tsv, converted_tsv, fdr, label_mapping, id_mapping=None):
    """Converts the MS2Rescore output to match the structure of a PeptideShaker report.

    Args:
        ms2rescore_tsv (string): Path to the MS2Rescore output.
        converted_tsv (string): Path where the converted report should be stored.
        fdr (float): False discovery rate in percent.
        label_mapping (dict): Dictionary containing the label of the entry. 
        id_mapping (string, optional): Path to the file containing the mapping of ORF lists and ORFs. Defaults to None.
    """
    fdr = float(fdr / 100)
    ms2rescore_df = pd.read_csv(ms2rescore_tsv, sep="\t", usecols=["peptidoform", "spectrum_id", "is_decoy", "qvalue", "pep", "precursor_mz", "protein_list", "provenance:xtandem_id"], dtype={'pep': float})
    filtered_ms2rescore_df = ms2rescore_df.loc[ms2rescore_df["is_decoy"] == False] # remove decoys
    filtered_ms2rescore_df = filtered_ms2rescore_df.loc[filtered_ms2rescore_df["qvalue"] < fdr] # filter for FDR
    filtered_ms2rescore_df["Confidence [%]"] = (1 - filtered_ms2rescore_df["pep"]) * 100 # compute confidence

    filtered_ms2rescore_df["peptidoform"] = filtered_ms2rescore_df["peptidoform"].apply(lambda x: removeModsFromSequence(x)) # remove mods from sequence
    

    filtered_ms2rescore_df.rename(columns={"peptidoform": "Sequence", "protein_list": "Protein(s)", "provenance:xtandem_id": "Spectrum Title"}, inplace=True)
    filtered_ms2rescore_df = filtered_ms2rescore_df[["spectrum_id", "Protein(s)", "Sequence", "precursor_mz", "Confidence [%]", "Spectrum Title"]]
    filtered_ms2rescore_df["Spectrum Title"] = filtered_ms2rescore_df["Spectrum Title"].astype(str)
    filtered_ms2rescore_df["Protein(s)"] = filtered_ms2rescore_df["Protein(s)"].apply(ast.literal_eval)
    filtered_ms2rescore_df["Spectrum Title"] = filtered_ms2rescore_df["Spectrum Title"].map(label_mapping)
    filtered_ms2rescore_df = filtered_ms2rescore_df.explode("Protein(s)", ignore_index=True)
    if id_mapping: # ORF Mapping, for the genomic and top scoring search, otherwise there is no mapping to be done here.
        mapping_df = pd.read_csv(id_mapping, sep="\t", header=None, names=["ID_Mapping", "ORF_List"])
        id_mapping = mapping_df.set_index("ID_Mapping")["ORF_List"].to_dict()
        filtered_ms2rescore_df["Protein(s)"] = filtered_ms2rescore_df["Protein(s)"].map(id_mapping)
    else:
        filtered_ms2rescore_df["Protein(s)"] = filtered_ms2rescore_df["Protein(s)"].apply(lambda x: extractIdentifier(x)) # extract ID
    
    # aggregate the rows again
    filtered_ms2rescore_df = filtered_ms2rescore_df.groupby("spectrum_id").agg({
                "Protein(s)": lambda x: ",".join(x),
                "Sequence": 'first',
                "precursor_mz": 'first',
                "Confidence [%]": 'first',
                "Spectrum Title": 'first',
            }).reset_index()

    # write DataFrame to file
    filtered_ms2rescore_df.to_csv(converted_tsv, sep="\t", header=True, index=False, quotechar="'")


def main():
    ms2rescore_tsv = snakemake.input[0]
    xtandem_xml = snakemake.input[1]
    converted_tsv = snakemake.output[0]
    fdr = snakemake.params[0]
    std_out = snakemake.log[1]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    label_mapping = parseXTandemXML(xtandem_xml)

    try: # genomic db search and top scoring db search
        orf_list_mapping = snakemake.input[2]
        logging.info(f"Using provided mapping for ORF lists.")
        convert_tsv(ms2rescore_tsv, converted_tsv, fdr, label_mapping, orf_list_mapping)
    except: # host filtering and reference db search
        logging.info(f"No ORF lists mapping provided. The protein list will not be changed.")
        convert_tsv(ms2rescore_tsv, converted_tsv, fdr, label_mapping)
    


if __name__ == "__main__":
    main()

