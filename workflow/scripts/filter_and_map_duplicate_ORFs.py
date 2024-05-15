import logging
from Bio import SeqIO
import pandas as pd


def filterUniqueORFs(in_file):
    """Filters the concatenated fasta file containing the proteomes and filters for duplicate sequences.
       If there is a duplicate sequence, the ORF ID is concatenated to the header of the "original" sequence.

    Args:
        in_file (string): Path to the concatenated fasta file.

    Returns:
        list: List of sequence records where duplicate are filtered.
    """
    
    sequences = {}
    sequences["orignial"] = {}
    sequences["decoys"] = {}
    id_mapping = {}
    counter = 0
    for record in SeqIO.parse(in_file, "fasta"):
        counter += 1
        current_id = record.id.split(" ")[0]
        if current_id.endswith("_REVERSED"):
            if record.seq not in sequences["decoys"]:
                sequences["decoys"][record.seq] = current_id
            else:
                r = current_id
                sequences["decoys"][record.seq] += f",{r}"
        else:
            if record.seq not in sequences["orignial"]:
                sequences["orignial"][record.seq] = current_id
            else:
                r = current_id
                sequences["orignial"][record.seq] += f",{r}"

    mapping_counter = 0
    for key in sequences["orignial"]:
        id_mapping[f"ORFs_{mapping_counter}"] = sequences["orignial"][key]
        sequences["orignial"][key] = f"ORFs_{mapping_counter}"
        mapping_counter += 1
    mapping_counter = 0
    for key in sequences["decoys"]:
        id_mapping[f"ORFs_{mapping_counter}_REVERSED"] = sequences["decoys"][key]
        sequences["decoys"][key] = f"ORFs_{mapping_counter}_REVERSED"
        mapping_counter += 1

    num_seqs = len(sequences["orignial"]) + len(sequences["decoys"])
    logging.info(f"Number of sequences: {counter}")
    logging.info(f"Number of unique sequences: {num_seqs}")
    return sequences, id_mapping


def writeFastaFile(sequences, out_file):
    """Writes the sequences with the modified sequence headers to a fasta file.

    Args:
        sequences (list): List of sequence records where duplicate are filtered.
        out_file (string): Path where the fasta file should be stored.
    """

    with open(out_file, "w") as f:
        for key in sequences:
            for sequence in sequences[key]:
                seq = str(sequence)
                f.write(">" + str(sequences[key][sequence]) + "\n")
                while seq:
                    f.write(seq[:60] + "\n")
                    seq = seq[60:]
                f.write("\n")


def writeMappingTsv(mapping, out_file):
    id_mapping_df = pd.DataFrame(list(mapping.items()), columns=['ID_Mapping', 'ORF_List'])
    id_mapping_df.to_csv(out_file, sep="\t", header=True, index=False, quotechar="'")


def main():
    in_file = snakemake.input[0]
    out_fasta = snakemake.output[0]
    out_tsv = snakemake.output[1]
    std_out = snakemake.log[1]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    unique_ORFs, id_mapping = filterUniqueORFs(in_file)
    writeFastaFile(unique_ORFs, out_fasta)
    writeMappingTsv(id_mapping, out_tsv)

if __name__ == "__main__":
    main()
