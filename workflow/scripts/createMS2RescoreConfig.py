# adapted from peptonizer 2000
def createConfig(outfile, fragmentation_method, fragment_tolerance, spectrum_pattern, threads):
    """Creates the configuration file for MS2Rescore.

    Args:
        outfile (string): Path where the MS2Rescore config is stored.
        fragmentation_method (string): Fragmentation method used in the experiment for this sample.
        fragment_tolerance (float): Fragment Ion MZ Tolerance used in the experiment for this sample.
        spectrum_pattern (string): Regex pattern to identify the spectrum titles.
        threads (int): Number of threads for MS2Rescore.
    """
    with open(outfile, "w") as f:
        lines = [' {\"$schema\":\"./config_schema.json\",']
        lines.append('"ms2rescore":{')
        lines.append('"feature_generators":{')
        lines.append('"basic":{},')
        lines.append('"ms2pip":{')
        lines.append(f'"model":"{fragmentation_method}",')
        lines.append(f'"ms2tolerance":{fragment_tolerance}')
        lines.append('}')
        lines.append('},')
        lines.append(f'"psm_id_pattern":"{spectrum_pattern}",')
        lines.append(f'"spectrum_id_pattern":"{spectrum_pattern}",')
        lines.append('"modification_mapping":{')
        lines.append('"+57":"Carbamidomethyl",')
        lines.append('"+15.99":"Oxidation",')
        lines.append('"+39.9945":"Pyro-Carbamidomethyl",')
        lines.append('"-17.0266":"Gln->pyro-Glu",')
        lines.append('"-18.0106":"Glu->pyro-Glu",')
        lines.append('"+42.0106":"Acetyl"')
        lines.append('},')
        lines.append('"fixed_modifications":{},')
        lines.append(f'"processes":{threads},')
        lines.append('"rename_to_usi": false,')
        lines.append('"fasta_file":null,')
        lines.append('"write_report": false,')
        lines.append('"rescoring_engine":{')
        lines.append('"percolator":{')
        lines.append('}},')
        lines.append('"psm_file_type" : "xtandem" ,')
        lines.append('"id_decoy_pattern": "_REVERSED"')
        lines.append('}}')

        f.writelines([line + "\n" for line in lines])


def main():
    outfile = snakemake.output[0]
    fragmentation_method = snakemake.params[0]
    fragment_tolerance = snakemake.params[1]
    spectrum_pattern = snakemake.params[2]
    threads = snakemake.params[3]
    
    createConfig(outfile, fragmentation_method, fragment_tolerance, spectrum_pattern, threads)


if __name__ == "__main__":
    main()
