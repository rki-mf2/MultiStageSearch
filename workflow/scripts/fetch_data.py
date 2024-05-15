import logging
import subprocess
import time
import pandas as pd
import numpy as np
from ete3 import NCBITaxa
from Bio import Entrez, SeqIO


def downloadCovidLineages(out_path, zip_name, APIkey):
    """Use the ncbi datasets query to download the metadata for all complete SARS-CoV-2 genomes.
    It downloads a zip file that is processed in formatCovidLineages().

    Args:
        out_path (string): Path where the zip file should be stored.
        zip_name (string): Name of the zip file.
        APIkey (string): The API key of your NCBI account (for the API).

    Raises:
        Exception: If something goes wrong while downloading.
    """
    with open(f"{out_path}/lineages_download.log", "w") as f:
        process = subprocess.Popen(
            [
                "datasets",
                "download",
                "virus",
                "genome",
                "taxon",
                f"2697049", # SARS-CoV-2
                "--complete-only",
                f"--include",
                "none",
                "--api-key",
                f"{APIkey}",
                "--filename",
                f"{out_path}/{zip_name}.zip"
            ],
            stdout=f,
            stderr=subprocess.STDOUT
        )
        process.wait()

        if process.returncode != 0:
            logging.info(f"downloadCovidLineages failed with code {process.returncode}!")
            raise Exception(f"downloadCovidLineages failed with code {process.returncode}!")
        else:
            logging.info("downloadCovidLineages finished successfully")
    

def formatCovidLineages(out_path, zip_name):
    """Formats the downloaded zip file so it can be further processed.

    Args:
        out_path (string): Path where the zip file should be stored.
        zip_name (string): Name of the zip file.

    Raises:
        Exception: If something goes wrong while formatting the zip (e.g. the file is corrupted).

    Returns:
        string: Path of the formatted tsv file.
    """
    tsv_file = f"{out_path}/{zip_name}.tsv"
    with open(tsv_file, "w") as f:
        process = subprocess.Popen(
        [
            "dataformat",
            "tsv",
            "virus-genome",
            "--fields",
            "accession,virus-pangolin",
            "--package",
            f"{out_path}/{zip_name}.zip",
        ],
        stdout=f,
        stderr=subprocess.STDOUT
        )
        process.wait()

        if process.returncode != 0:
            logging.info(f"formatCovidLineages failed with code {process.returncode}!")
            raise Exception(f"formatCovidLineages failed with code {process.returncode}!")
        else:
            logging.info("formatCovidLineages finished successfully")
            return tsv_file
    

def parseCovidTSV(tsv_file):
    """This function parses the formatted tsv and decides which lineages are used for the analysis.

    Args:
        tsv_file (string): Path of the formatted tsv file.

    Returns:
        dict: Mapping of the lineage and the corresponding genbank accessions.
        list: List of the genbank accessions used for the analysis.
    """
    lineages_df = pd.read_csv(tsv_file, sep="\t")
    lineages_df.rename(columns = {"Virus Pangolin Classification": "Lineage"}, inplace=True)
    lineages_df.dropna(axis=0, how="any", inplace=True)
    lineages_df[["Accession", "Lineage"]] = lineages_df[["Accession", "Lineage"]].astype(str)
    filtered_lineages_df = lineages_df[lineages_df["Lineage"].str.lower() != "unclassifiable"] # remove entries without lineage
    #filtered_lineages_df = filtered_lineages_df.sort_values(by="Length", ascending=False) # not usefull, since longer sequences could just mean more Ns, Example: OX448934.1

    big_lineages = filtered_lineages_df[filtered_lineages_df["Lineage"].str.count("\.") <= 1] # only use "major lineages" and "major sublineages" (with at most 1 dot, e.g. "B" and "B.1")
    big_unique_lineages = big_lineages["Lineage"].unique().tolist()

    genbank_accessions_lineage_mapping = {}
    genbank_accessions_lineage_mapping[2697049] = {}
    for lineage in big_unique_lineages:
        current_lineage_df = filtered_lineages_df.loc[filtered_lineages_df["Lineage"] == lineage]
        genbank_accessions_lineage_mapping[2697049][lineage] = current_lineage_df["Accession"].iloc[0] # take the first genbank accession for each lineage
    covid_genbank_accessions = list(genbank_accessions_lineage_mapping[2697049].values())
    
    return genbank_accessions_lineage_mapping, covid_genbank_accessions



def fetchCovidLineages(out_path, APIkey, max_sequence_length, sequence_length_diff):
    """Performs the fetching of SARS-CoV-2 genomes.

    Args:
        out_path (string): Path where the downloaded zip file containing the metadata should be stored.
        APIkey (string): The API key of your NCBI account (for the API).
        max_sequence_length (interger): the maximum lenght of the sequences, 
                                        all sequences above this threshold are ignored.
        sequence_length_diff (interger): the length difference threshold,
                                        if the lengths of two following sequences is to big, the search stops.
                                        The sequences are sorted by length (descending).

    Returns:
        list: List of all fetch genome records.
        dict: Mapping of the lineage and the corresponding genbank accessions.
    """
    zip_name = "SARS_CoV_2_complete_genomes"
    
    downloadCovidLineages(out_path, zip_name, APIkey)
    time.sleep(30) # for file system latency
    tsv_file = formatCovidLineages(out_path, zip_name)
    genbank_accessions_lineage_mapping, covid_genbank_accessions = parseCovidTSV(tsv_file)
    all_records = ncbiFetchGenomes(covid_genbank_accessions, max_sequence_length, sequence_length_diff, use_NCBI_taxonomy=False)
    return all_records, genbank_accessions_lineage_mapping


def getStrainNames(descendant_names, all_taxonomy_names):
    """get the strain names from the whole title

    Args:
        descendant_names (list): list of names from the NCBI taxonomy

    Returns:
        list: list of the strain names
    """

    strain_names = []
    for name in descendant_names:
        name = str(name)
        try:
            strain_names.append(name.split("strain ")[1])
        # if the name of the descendent has no "strain" in it, try to split after the name of the sample.
        # Do it for every alternative name found in the NCBI Taxonomy.
        except IndexError:
            for alternative_name in all_taxonomy_names:
                try:
                    alternative_name = str(alternative_name)
                    strain_names.append(name.split(f"{alternative_name} ")[1])
                except IndexError:
                    logging.info(f"Using full name for {name}, since no strain name could be found!")
                    strain_names.append(name)

    strain_names = [name.replace("(", "").replace(")", "") for name in strain_names]
    unique_strain_names = []
    for name in strain_names:
        if name not in unique_strain_names:
            unique_strain_names.append(name)
    return unique_strain_names


def getAllTaxonomyNames(taxid):
    """find all names of a species that can be used for querying the database

    Args:
        taxid (integer): taxon id of a species

    Returns:
        list: list of all names of a species
    """

    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    record = Entrez.read(handle)

    record_names = []
    record_names.append(record[0]["ScientificName"])
    for entry in record[0]["OtherNames"]["Synonym"]:
        record_names.append(entry)
    for entry in record[0]["OtherNames"]["EquivalentName"]:
        record_names.append(entry)
    return record_names


def ncbiFetchGenomes(genome_id_list, max_sequence_length, sequence_length_diff, mapping = None, use_NCBI_taxonomy=False):
    """ Fetches genomes with the genome ids in the NCBI nucleotide database.

    Args:
        genome_id_list (list): list of genome ids to query
        max_sequence_length (interger): the maximum lenght of the sequences, 
                                        all sequences above this threshold are ignored
        sequence_length_diff (interger): the length difference threshold,
                                        if the lengths of two following sequences is to big, the search stops.
                                        The sequences are sorted by length (descending)

    Returns:
        list: list of all fetch genome records
    """
    
    all_records = []
    used_records = []

    num_used_seqs = 0
    previous_seq_length = None
    for id in range(len(genome_id_list)):
        handle = Entrez.efetch(db="nucleotide", id=genome_id_list[id], rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        logging.info(f"{id+1}/max{len(genome_id_list)}; Sequence ID: {genome_id_list[id]}, Sequence Lenght: {len(record.seq)}")
        if (not use_NCBI_taxonomy) and previous_seq_length and (previous_seq_length >= len(record.seq) * sequence_length_diff):
            logging.info("Sequence length difference exceeds threshold, query ends here!")
            break
        previous_seq_length = len(record.seq)
        all_records.append(record)
        if mapping:
            for species_id in mapping.keys():
                if genome_id_list[id] in mapping[species_id].keys():
                    mapping[species_id][genome_id_list[id]] = record.id

    all_records = sorted(all_records, key=lambda x: len(x.seq), reverse=True)

    for record in all_records:
        if len(record.seq) > max_sequence_length:
            logging.info(f"Sequence with ID {genome_id_list[id]} is longer than allowed. It will be skipped!")                
            continue
        elif (num_used_seqs > 0) and (len(all_records[num_used_seqs-1].seq) >= (len(record.seq) * sequence_length_diff)):
            break
        else:
            used_records.append(record)
            num_used_seqs += 1
    logging.info(f"Only the {num_used_seqs} longest are used. All others are below the sequence length difference threshold!")

    if mapping:
        return used_records, mapping
    return used_records


def getTaxIdsToQuery(mapped_taxids, number_taxids, weight_diff):
    """Get taxon ids from the first search that are used for querying.

    Args:
        mapped_taxids (string): path to the taxid scores file
        number_taxids (integer): maximum number of taxon ids to consider
        weight_diff (integer): weight difference,
                                if the difference between two taxonids (sorted) is above the threshold,
                                the taxon id with the lower score and all below are not considered.

    Returns:
        list: list of taxon ids considered for querying
    """

    taxid_df = pd.read_csv(mapped_taxids, sep="\t", header=0, index_col=False)
    taxids = taxid_df["taxid"].to_list()
    relevant_taxids = taxids[:number_taxids]
    taxids_to_query = []

    for i in range(len(relevant_taxids)):
        if i == 0:
            taxids_to_query.append(relevant_taxids[i])
        else:
            previous_taxid = taxid_df.loc[taxid_df["taxid"] == relevant_taxids[i-1]]
            previous_taxid_score = previous_taxid["weight"].values[0]
            current_taxid = taxid_df.loc[taxid_df["taxid"] == relevant_taxids[i]]
            current_taxid_score = current_taxid["weight"].values[0]

            if previous_taxid_score <= (current_taxid_score * weight_diff):
                taxids_to_query.append(relevant_taxids[i])
            else:
                break

    logging.info(f"TaxIDs used: {taxids_to_query}")
    return taxids_to_query


def fetchData(taxids_to_query, max_number_accessions, max_sequence_length, sequence_length_diff, genome_query):
    """Retrieve genome records from the NCBI database.

    Args:
        taxids_to_query (list): taxon ids that are used for querying
        max_number_accessions (integer): maximum number of genomes to consider
        max_sequence_length (interger): the maximum lenght of the sequences, 
                                        all sequences above this threshold are ignored
        sequence_length_diff (interger): the length difference threshold,
                                        if the lengths of two following sequences is to big, the search stops.
                                        The sequences are sorted by length (descending)
        genome_query (string): either an empty string or "complete" --> only search for complete genomes or all genomes

    Returns:
        list: list of all genome records found
    """
    mapping = {}
    all_records = []
    genome_id_mapping = {}
    for taxid in taxids_to_query:
        mapping[taxid] = {}
        logging.info(f"Querying for TaxId: {taxid}")
        search_query = f'txid{taxid}[Organism]{genome_query}'
        search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_number_accessions, sort="Sequence Length")
        taxon_id_list = Entrez.read(search_handle)["IdList"]
        for id in taxon_id_list:
            mapping[taxid][id] = np.nan
        
        try:
            all_taxid_records, genome_taxid_mapping = ncbiFetchGenomes(taxon_id_list, max_sequence_length, sequence_length_diff, mapping, use_NCBI_taxonomy = False)
            all_records += all_taxid_records
            genome_id_mapping.update(genome_taxid_mapping)
        except TypeError:
            logging.info(f"No Sequences found for {taxid}!")
        
    return all_records, genome_id_mapping


def fetchDataNCBI(taxids_to_query, max_number_accessions, max_sequence_length, sequence_length_diff, genome_query, sqlite_db_path):
    """"Retrieve genome records from the NCBI database.
        In this function only strains are considered that are present in the NCBI taxonomy.
        This function is used when MSS is called within PepGM

    Args:
        taxids_to_query (list): taxon ids that are used for querying
        max_number_accessions (integer): maximum number of genomes to consider
        max_sequence_length (interger): the maximum lenght of the sequences, 
                                        all sequences above this threshold are ignored
        sequence_length_diff (interger): the length difference threshold,
                                        if the lengths of two following sequences is to big, the search stops.
                                        The sequences are sorted by length (descending)
        genome_query (string): either an empty string or "complete" --> only search for complete genomes or all genomes
        sqlite_db_path (string): path to the sqlite database file created by ete3. This can be required if the default path of ete3 is not accessible.

    Returns:
        list: list of all genome records found
        dict: dictionary that contains 
    """
    if sqlite_db_path != "None":
        ncbi = NCBITaxa(dbfile = sqlite_db_path)
    else:
        ncbi = NCBITaxa()

    logging.info("Using the NCBI Querying!")
    for taxid in taxids_to_query:
        logging.info(f"Querying for TaxId: {taxid}")
        all_taxonomy_names = getAllTaxonomyNames(taxid)
        descendants = ncbi.get_descendant_taxa(taxid)
        descendant_names = ncbi.translate_to_names(descendants)
        strain_names = getStrainNames(descendant_names, all_taxonomy_names)
        name_id_mapping = {key: (value1, value2, taxid) for key, value1, value2 in zip(descendants, descendant_names, strain_names)}
        genome_id_list = []

        for species_name in all_taxonomy_names:
            for strain_name in strain_names:
                search_query = f'{species_name} AND ((strain {strain_name}[All Fields]) or ({strain_name}[Title])){genome_query}'
                logging.info(f"Querying for Taxon Name: {search_query}")
                search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_number_accessions, sort="Sequence Length")
                taxon_ids = Entrez.read(search_handle)["IdList"]
                for taxon_id in taxon_ids:
                    genome_id_list.append(taxon_id)
        genome_id_list = list(set(genome_id_list))
        try:
            all_records = ncbiFetchGenomes(genome_id_list, max_sequence_length, sequence_length_diff, use_NCBI_taxonomy=True)
        except TypeError:
            logging.info(f"No Sequences found for {taxid}!")
    return all_records, name_id_mapping


def getStrainNameDF(all_records, mapping, NCBI):
    """Mapping the genomes to the strains and create a pandas DataFrame that contains unique 
       combinations of genome, species, strain and isolate.
       This function also creates a list containing all sequence records of the genomes in the Pandas DataFrame.

    Args:
        all_records (list): list of all genome records queried
        mapping (dictionary): Dictionary that maps the decendant taxon id, 
                                        descendant name and the taxon id of the ancestor.
        NCBI (boolean): True if the NCBI Taxonomy is used, False if not.

    Returns:
        pandas DataFrame: DataFrame containing unique 
                          combinations of genome, species, strain and isolate.
        list: containing all sequence records of the genomes in the Pandas DataFrame

    """
    df = pd.DataFrame(columns=["genbank_accession", "species", "strain", "isolate", "taxa", "HigherTaxa"])
    # SARS-CoV-2 handling
    if 2697049 in list(mapping.keys()):
        covid_taxon_id = 6000000 # start very high to avoid real taxon ids
        for record in all_records:
            accession = record.id
            species = record.features[0].qualifiers["organism"][0]
            for key, value in mapping[2697049].items(): # get the correct lineage for the sequence record
                if value == accession:
                    strain = key
                    break
            isolate = np.nan # could be filled, but omitted
            new_df = pd.DataFrame({"genbank_accession": accession, "species": species, "strain": strain, "isolate": isolate, "taxa": covid_taxon_id, "HigherTaxa": 2697049}, index=[0])
            df = pd.concat([df, new_df])
            covid_taxon_id += 1
        unique_df = df.drop_duplicates(subset=["genbank_accession", "species", "strain", "isolate"], keep="first").reset_index(drop=True)
    # handling of other species
    else:       
        # get df with strain names...
        for record in all_records:
            accession = record.id
            species = record.features[0].qualifiers["organism"][0]
            try:
                strain = record.features[0].qualifiers["strain"][0]
            except KeyError:
                strain = np.nan
            try:
                isolate = record.features[0].qualifiers["isolate"][0]
            except KeyError:
                isolate = np.nan
            new_df = pd.DataFrame({"genbank_accession": accession, "species": species, "strain": strain, "isolate": isolate}, index=[0])
            df = pd.concat([df, new_df])

        # get unique strains...
        unique_df = df.drop_duplicates(subset=["genbank_accession", "species", "strain", "isolate"], keep="first").reset_index(drop=True)
        # add taxon ids for strain if a mapping is provided and higher taxon id (ancestor)
        # when using the NCBI taxonomy
        if NCBI:
            for key in mapping:
                unique_df.loc[unique_df['species'].str.lower() == str(mapping[key][1]).lower(), ["taxa", "HigherTaxa"]] = (str(key), mapping[key][2])
                unique_df.loc[unique_df['strain'].str.lower() == str(mapping[key][1]).lower(), ["taxa", "HigherTaxa"]] = (str(key), mapping[key][2])
                unique_df.loc[(unique_df['isolate'].str.lower() == str(mapping[key][1]).lower()) & (unique_df['strain'] == np.nan), ["taxa", "HigherTaxa"]] = (str(key), mapping[key][2])

            # only use strains with a taxon id (for PepGM)
            number_genomes = unique_df.shape[0]
            unique_df = unique_df.dropna(subset="taxa").reset_index(drop=True)
            logging.info(f"Dropped {number_genomes - unique_df.shape[0]} genomes. No taxa was found for them!")
        # when not using the NCBI taxonomy
        else:
            number_genomes = unique_df.shape[0]
            unique_df = unique_df.dropna(subset=["strain", "isolate"], how="all").reset_index(drop=True) # should be specific for non NCBI Taxonomy
            logging.info(f"Dropped {number_genomes - unique_df.shape[0]} genomes. They had neither strain nor isolate information!")
            taxon_id = 5000000 # start very high to avoid real taxon ids
            for species_id in mapping:
                for genbank_id in mapping[species_id]:
                    unique_df.loc[unique_df["genbank_accession"] == mapping[species_id][genbank_id], ["taxa", "HigherTaxa"]] = (str(taxon_id), species_id)
                    taxon_id += 1

    strain_records = [record for record in all_records if (record.id in unique_df["genbank_accession"].to_list())]
    # filter for duplicate record ids in sequences. This can happend if the search finds the same sequences for different TaxonIDs (e.g. 10515 and 129951)
    record_ids = []
    unique_strain_records = []
    for record in strain_records:
        if record.id in record_ids:
            continue
        else:
            record_ids.append(record.id)
            unique_strain_records.append(record)

    return unique_df, unique_strain_records

def writeFiles(unique_df, unique_strain_records, activate_covid_mode, out_file_df, out_file_fasta, check_covid_mode_file):
    """Writes the pandas DataFrame and the sequence records to files for further downstream analysis.

    Args:
        unique_df (pandas DataFrame): DataFrame containing unique 
                          combinations of genome, species, strain and isolate.
        unique_strain_records (list): Containing all sequence records of the genomes in the Pandas DataFrame.
        activate_covid_mode (string): Path where a file will be stored that is used to decide whether the special SARS-CoV-2 Search approach is used.
        out_file_df (string): Path to where the csv file is created
        out_file_fasta (string): Path to where the fasta file is created
    """
    
    # filter problematic records
    record_to_filter = []
    for record in range(len(unique_strain_records)):
        try:
            unique_strain_records[record].seq[0]
        except:
            record_to_filter.append(record)

    record_to_filter.sort(reverse=True)
    for record in record_to_filter:
        logging.info(f"removing problematic sequence record {unique_strain_records[record]}")
        unique_strain_records.pop(record)
    # write files...
    unique_df.to_csv(out_file_df, index=False, sep=",")
    SeqIO.write(unique_strain_records, out_file_fasta, "fasta")
    with open(check_covid_mode_file, "w") as f:
        f.write(activate_covid_mode)


def main():
    # get all input, parameters and output paths
    mapped_taxids = snakemake.input[0]
    out_file_df = snakemake.output[0]
    out_file_fasta = snakemake.output[1]
    check_covid_mode_file = snakemake.output[2]
    APImail = snakemake.params[0]
    APIkey = snakemake.params[1]
    number_taxids = snakemake.params[2]
    weight_diff = snakemake.params[3]
    max_sequence_length = snakemake.params[4]
    sequence_length_diff = snakemake.params[5]
    max_number_accessions = snakemake.params[6]
    use_NCBI_Taxa = snakemake.params[7]
    only_use_complete_genomes = snakemake.params[8]
    sqlite_db_path = snakemake.params[9]
    out_path = snakemake.params[10]
    std_out = snakemake.log[1]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    # api info used for Entrez
    if APIkey and APImail:
        Entrez.email = APImail
        Entrez.api_key = APIkey

    # get input for query
    taxids_to_query = getTaxIdsToQuery(mapped_taxids, number_taxids, weight_diff)

    # either "complete genome" or "genome" in title
    if only_use_complete_genomes:
        genome_query = ' AND "complete "[Title]'
    else:
        genome_query = ""

    activate_covid_mode = "False"

    # which query to use
    if 2697049 in taxids_to_query: # SARS-CoV-2
        all_records, genbank_accessions_lineage_mapping = fetchCovidLineages(out_path, APIkey, max_sequence_length, sequence_length_diff)
        unique_df, unique_strain_records = getStrainNameDF(all_records, genbank_accessions_lineage_mapping, NCBI=False)
        activate_covid_mode = "True"
    elif use_NCBI_Taxa:
        all_records, mapping = fetchDataNCBI(taxids_to_query, max_number_accessions, max_sequence_length, sequence_length_diff, genome_query, sqlite_db_path)
        if len(all_records) == 0:
            logging.info("No Sequences found usign this query. Trying other query!")
            all_records, mapping = fetchData(taxids_to_query, max_number_accessions, max_sequence_length, sequence_length_diff, genome_query)
            unique_df, unique_strain_records = getStrainNameDF(all_records, mapping, NCBI=False)
        else:
            unique_df, unique_strain_records = getStrainNameDF(all_records, mapping, NCBI=True)
    else:
        all_records, mapping = fetchData(taxids_to_query, max_number_accessions, max_sequence_length, sequence_length_diff, genome_query)
        unique_df, unique_strain_records = getStrainNameDF(all_records, mapping, NCBI=False)
    
    # write output
    writeFiles(unique_df, unique_strain_records, activate_covid_mode, out_file_df, out_file_fasta, check_covid_mode_file)


if __name__ == "__main__":
    main()