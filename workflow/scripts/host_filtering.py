# from PepGM

import pandas as pd


def getSpectraNames(report):

    """
    get spectrum titles from PeptideShaker output
    input: path to spectrum file
    returns: list of spectra tiltes
    """

    raw = pd.read_csv(report, sep="\t", on_bad_lines="skip", usecols=["Spectrum Title"])
    spectrumNames = raw["Spectrum Title"].to_list()
    return spectrumNames


def FilterMGF(SpectraToFilter, SpectrumFile, out_file, log_file):

    """
    filter spectra with title provided in a list from .mgf file
    input: list of spectrum titles, .mgf file of spectra
    output: filtered .mgf file
    """

    SpectraToFilter = set(SpectraToFilter)  # remove duplicates
    with open(SpectrumFile, "r", newline="\n") as mgf:
        LinesToWrite = []
        writeLine = True

        for line in mgf:
            lineNoN = line.rstrip()
            if lineNoN.startswith("TITLE"):
                if lineNoN[6:] in SpectraToFilter:
                    writeLine = False
                else:
                    writeLine = True
            if writeLine:
                LinesToWrite.append(line)

    with open(out_file, "w") as f:
        f.writelines(LinesToWrite)

    with open(log_file, "w") as f:
        f.write("filtered " + str(len(SpectraToFilter)) + " host spectra")


def main():
    mgf = snakemake.input[0]
    psm_report = snakemake.input[1]
    out_file = snakemake.output[0]
    log_file = snakemake.log[0]

    FilterList = getSpectraNames(psm_report)
    FilterMGF(FilterList, mgf, out_file, log_file)


if __name__ == "__main__":
    main()
