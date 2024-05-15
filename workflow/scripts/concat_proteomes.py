import os
import shutil


def concatProteomes(concat_proteomes_file, result_dir, sample_name, sixpack_dir_name):
    """Reads in all proteome fasta files created by sixpack and concatenates them into one file.

    Args:
        concat_proteomes_file (string): Path where the concatened proteomes file should be stored.
        result_dir (string): Path to the directory containing the sixpack proteome files.
        sample_name (string): Name of the sample.
    """
    
    cwd = str(os.getcwd())
    full_in_path = f"{cwd}/{result_dir}/{sample_name}/{sixpack_dir_name}"

    included_extensions = ["fasta"]
    file_names = [
        fn for fn in os.listdir(full_in_path) if any(fn.endswith(ext) for ext in included_extensions)
    ]

    with open(concat_proteomes_file, "w") as out_file:
        for file in file_names:
            with open(f"{full_in_path}/{file}", "r") as f:
                shutil.copyfileobj(f, out_file)


def main():
    concat_proteomes = snakemake.output[0]
    res_dir = str(snakemake.params[0])
    sample_name = str(snakemake.params[1])
    sixpack_dir_name = str(snakemake.params[2])
    concatProteomes(concat_proteomes, res_dir, sample_name, sixpack_dir_name)


if __name__ == "__main__":
    main()
