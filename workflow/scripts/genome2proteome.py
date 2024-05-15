import subprocess
import time
import shutil
import os


def getCorrectFasta(fastas, lineage):
    """If using the SARS-CoV-2 mode, find the correct fasta for the current lineage.

    Args:
        fastas (list): List of all paths to the fasta files for the SARS-CoV-2 lineages.
        lineage (string): Current lineage.

    Returns:
        string: Path to the correct fasta for the current lineage.
    """
    for fasta in fastas:
        if f"/{lineage}/" in fasta:
            return fasta
        

def getSequences(sixpack_temp_fasta):
    """Gets a list of all sequences in the concatenated fasta file.

    Args:
        sixpack_temp_fasta (string): Path to the concatenated fasta file.

    Returns:
        list: List of all sequences
    """
    
    with open(sixpack_temp_fasta, "r") as fasta:
        sequences = fasta.read()
        d = ">"
        sequences = [d+e for e in sequences.split(d)][1:]
    return sequences


def copyFasta(concat_fasta, sixpack_temp_fasta, sleep_time):
    shutil.copy(concat_fasta, sixpack_temp_fasta)
    time.sleep(sleep_time)


def deleteTempFasta(sixpack_temp_fasta):
    os.remove(sixpack_temp_fasta)


def runSixpack(sixpack_out, sequences, additional_sixpack_params, sixpack_temp_fasta, sample_path, orfminsize, sleep_time):
    """Runs sixpack for each sequence in the fasta file to perform a sixframe translation.

    Args:
        sixpack_out (string): Path where the standard output of sixpack is stored.
        sequences (list): List of all sequences
        additional_sixpack_params (string): Additional sixpack parameters if provided.
        sixpack_temp_fasta (string): Path to the concatenated fasta file.
        sample_path (string): Path where sixpack stores the resulting fasta files.
        orfminsize (integer): The minimum number of amino acids that an ORF should have.
    """
    # workaround since sixpack is not able to read the first sequence otherwise
    with open(sixpack_temp_fasta, "w") as fasta:
        for seq in sequences:
            fasta.write(seq + "\n")
    time.sleep(sleep_time)  # because file system latency

    with open(sixpack_out, "w") as f:
        for i, sequence in enumerate(sequences):
            if additional_sixpack_params:
                process = subprocess.Popen(
                    [
                        "sixpack",
                        "-sequence",
                        sixpack_temp_fasta,
                        "-outfile",
                        f"{sample_path}{sequence[1:11]}.orfs",
                        "-outseq",
                        f"{sample_path}{sequence[1:11]}.fasta",
                        "-orfminsize",
                        f"{orfminsize}",
                        f"{additional_sixpack_params}"
                    ],
                    stdout=f,
                    stderr=subprocess.STDOUT
                )
                process.wait()
            else:
                process = subprocess.Popen(
                    [
                        "sixpack",
                        "-sequence",
                        sixpack_temp_fasta,
                        "-outfile",
                        f"{sample_path}{sequence[1:11]}.orfs",
                        "-outseq",
                        f"{sample_path}{sequence[1:11]}.fasta",
                        "-orfminsize",
                        f"{orfminsize}",
                    ],
                    stdout=f,
                    stderr=subprocess.STDOUT
                )
                process.wait()
            
            # remove the sequence from the concatenated fasta file, that was used in this iteration.
            with open(sixpack_temp_fasta, "w") as fasta:
                for seq in sequences[i+1:]:
                    fasta.write(seq + "\n")
            time.sleep(sleep_time)  # because file system latency


def main():
    concat_fastas = snakemake.input
    sixpack_out = snakemake.output[1]

    sample_path = snakemake.params[0] + "/" + snakemake.params[3] + "/"
    orfminsize = snakemake.params[1]
    additional_sixpack_params = snakemake.params[2]
    sixpack_temp_fasta = snakemake.params[4]
    sleep_time = int(snakemake.params[5])
    
    try: # if the SARS-CoV-2 mode is used, there will be a second wildcard
        lineage = snakemake.wildcards[1]
        concat_fasta = getCorrectFasta(concat_fastas, lineage)
    except: # normal mode
        concat_fasta = snakemake.input[0]

    copyFasta(concat_fasta, sixpack_temp_fasta, sleep_time)
    sequences = getSequences(concat_fasta) 
    deleteTempFasta(sixpack_temp_fasta)
    runSixpack(sixpack_out, sequences, additional_sixpack_params,
                sixpack_temp_fasta, sample_path, orfminsize, sleep_time)


if __name__ == "__main__":
    main()
