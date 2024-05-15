import os
import glob 

class General:
    @staticmethod
    def get_input_all(wildcards):
        input_list = []
        input_list += expand(RESULT_DIR / "{sample}/Plots/create_plots.txt", sample=list(SAMPLES.index))       
        input_list += expand(RESULT_DIR / "{sample}/PepGMInput/MultiStageSearchResults.csv", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/PepGMInput/MultiStageSearchScores.csv", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/PepGMInput/MultiStageSearchStrainAccessions.csv", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_peptidome_heatmap.png", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/config.yaml", sample=list(SAMPLES.index))
        
        if COMPUTE_PHYLOGENY:
            input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_phylogeny.png", sample=list(SAMPLES.index))
            input_list += expand(RESULT_DIR / "{sample}/PepGMInput/phylogenetic_tree.nwk", sample=list(SAMPLES.index))
        
        if COMPUTE_SIMILARITY_MATRIX:
            input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_similarity_heatmap.png", sample=list(SAMPLES.index))

        input_list += expand(Plots.get_current_working_directory() / RESULT_DIR / "{sample}/ResultsReport.html", sample=list(SAMPLES.index))
       
        return input_list

class SearchDB:
    @staticmethod
    def get_input_MGF():
        if HOST_FILTERING:
            data_dict = {
                "mgf": RESULT_DIR / "{sample}/SpectraFilter/{sample}.mgf"
                }
        else:
            data_dict = {
                "mgf": MGF_FILE,
                }
        return data_dict
    

    @staticmethod
    def get_output_AddDecoysRef():
        ref_path = config["db_search"]["ref"]
        ref_path = ref_path.split(".")
        ref_name = ""
        num_name_segments = len(ref_path[:-1])
        for i in range(num_name_segments):
            if i == num_name_segments - 1:
                ref_name += ref_path[i]
            else:
                ref_name += ref_path[i] + "."

        ref_name += "_concatenated_target_decoy." + ref_path[-1]

        data_dict = {
            "ref_decoy_fasta": ref_name 
        }
        return data_dict

    
    @staticmethod
    def getSearchguiFDRsDBSearch():
        if config["db_search"]["scoring_engine"] == "PeptideShaker":
            data_dict = {
                "psm_fdr": config["db_search"]["peptide_shaker_params"]["psm_fdr"],
                "peptide_fdr": config["db_search"]["peptide_shaker_params"]["peptide_fdr"],
                "protein_fdr": config["db_search"]["peptide_shaker_params"]["protein_fdr"]
            } 
            return data_dict
        elif config["db_search"]["scoring_engine"] == "MS2Rescore":
            data_dict = {
                "psm_fdr": 100,
                "peptide_fdr": 100,
                "protein_fdr": 100
            } 
            return data_dict
        else:
            print('Please choose one of the following: "MS2Rescore", "PeptideShaker"')
            return None


    @staticmethod
    def getReports():
        data_dict = {}
        if config["db_search"]["scoring_engine"] == "PeptideShaker":
            data_dict["host_filtering_report"] = str(RESULT_DIR / "{sample}/SpectraFilter/host_Default_PSM_Report.txt")
            data_dict["first_search_report"] = str(RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt")
            data_dict["final_search_report"] = str(RESULT_DIR / "{sample}/FinalSearch/converted_proteomes_Default_PSM_Report.txt")
            data_dict["extra_search_report"] = str(RESULT_DIR / "{sample}/ExtraSearch/converted_proteomes_Default_PSM_Report.txt")
            data_dict["covid_search_report"] = str(RESULT_DIR / "{sample}/CovidMode/{lineage}/converted_proteomes_Default_PSM_Report.txt")
            data_dict["DB_suitability_host_filtering_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_host_Default_PSM_Report.txt")
            data_dict["DB_suitability_first_search_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_Default_PSM_Report.txt")
            data_dict["DB_suitability_final_search_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_converted_proteomes_Default_PSM_Report.txt")
            data_dict["DB_suitability_extra_search_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_converted_proteomes_Default_PSM_Report.txt")
        elif config["db_search"]["scoring_engine"] == "MS2Rescore":
            data_dict["host_filtering_report"] = str(RESULT_DIR / "{sample}/SpectraFilter/converted_ms2rescored.tsv")
            data_dict["first_search_report"] = str(RESULT_DIR / "{sample}/FirstSearch/converted_ms2rescored.tsv")
            data_dict["final_search_report"] = str(RESULT_DIR / "{sample}/FinalSearch/converted_ms2rescored.tsv")
            data_dict["extra_search_report"] = str(RESULT_DIR / "{sample}/ExtraSearch/converted_ms2rescored.tsv")
            data_dict["covid_search_report"] = str(RESULT_DIR / "{sample}/CovidMode/{lineage}/converted_ms2rescored.tsv")
            data_dict["DB_suitability_host_filtering_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_converted_ms2rescored.tsv")
            data_dict["DB_suitability_first_search_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_converted_ms2rescored.tsv")
            data_dict["DB_suitability_final_search_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_converted_ms2rescored.tsv")
            data_dict["DB_suitability_extra_search_report"] = str(RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_converted_ms2rescored.tsv")
        else:
            print('Please choose one of the following: "MS2Rescore", "PeptideShaker"')
            return None

        return data_dict
    


class PepGM:
    @staticmethod
    def get_input_CreatePepGMInput():
        input_list = []
        if not EXTRA_SEARCH:
            input_list.append(SearchDB.getReports()["final_search_report"])
            input_list.append(RESULT_DIR / "{sample}/FetchData/strain_accessions.csv")
            input_list.append(RESULT_DIR / "{sample}/taxids/ORF_scores.csv")       
        else:
            input_list.append(SearchDB.getReports()["extra_search_report"])
            input_list.append(RESULT_DIR / "{sample}/FetchData/strain_accessions.csv")
            input_list.append(RESULT_DIR / "{sample}/taxids/extra_search_ORF_scores.csv")  
        return input_list


class Plots:
    @staticmethod
    def get_input_PeptidomeHeatmap():
        if not EXTRA_SEARCH:
            data_dict = {
                "report": SearchDB.getReports()["final_search_report"],
                }
        else:
            data_dict = {
                "report": SearchDB.getReports()["extra_search_report"],
                }
        return data_dict

    @staticmethod
    def get_input_Plots():
        input_list = []
        input_list.append(RESULT_DIR / "{sample}/taxids/strain_name_counts.tsv")
        input_list.append(RESULT_DIR / "{sample}/taxids/taxid_scores.tsv")
        input_list.append(SearchDB.getReports()["first_search_report"])
        input_list.append(SearchDB.getReports()["final_search_report"])
        input_list.append(RESULT_DIR / "{sample}/taxids/ORF_scores.csv")
        if HOST_FILTERING:
            input_list.append(SearchDB.getReports()["host_filtering_report"])
        if EXTRA_SEARCH:
            input_list.append(RESULT_DIR / "{sample}/taxids/extra_search_strain_name_counts.tsv")
            input_list.append(SearchDB.getReports()["extra_search_report"])
            input_list.append(RESULT_DIR / "{sample}/taxids/extra_search_ORF_scores.csv")
        
        if COMPUTE_DB_SUITABILITY:
            if HOST_FILTERING:
                input_list.append(SearchDB.getReports()["DB_suitability_host_filtering_report"])
            input_list.append(SearchDB.getReports()["DB_suitability_first_search_report"])
            input_list.append(SearchDB.getReports()["DB_suitability_final_search_report"])
            if EXTRA_SEARCH:
                input_list.append(SearchDB.getReports()["DB_suitability_extra_search_report"])
        
        return input_list

    def get_current_working_directory():
        cwd = str(os.getcwd())
        return cwd


class CovidMode:
    @staticmethod
    def get_input_Covidgenome2Proteome(wildcards):
        checkpoint_output = checkpoints.splitLineages.get(**wildcards).output[0]
        subfolders = glob.glob(os.path.join(checkpoint_output, "*"))
        files = []
        for folder in subfolders:
            if os.path.isdir(folder):
                file = os.path.join(folder, os.path.basename(folder))
                files.append(file)
        lineages = [file.split("/")[-1] for file in files]
        return expand(RESULT_DIR / "{sample}/CovidMode/Lineages/{lineage}/{lineage}.fasta", sample=list(SAMPLES.index), lineage=lineages)


    @staticmethod
    def get_Covid_topScoring(wildcards):
        checkpoint_output = checkpoints.splitLineages.get(**wildcards).output[0]
        subfolders = glob.glob(os.path.join(checkpoint_output, "*"))
        files = []
        for folder in subfolders:
            if os.path.isdir(folder):
                file = os.path.join(folder, os.path.basename(folder))
                files.append(file)
        lineages = [file.split("/")[-1] for file in files]
        return expand(RESULT_DIR / "{sample}/CovidMode/{lineage}/taxids/ORF_scores.csv", sample=list(SAMPLES.index), lineage=lineages)


    def checkForCovidMode(wildcards):
        activation_file = checkpoints.fetchStrainGenomes.get(**wildcards).output[2]
        with open(activation_file, "r") as f:
            line = f.readline()
            if line == "True":
                return RESULT_DIR / "{sample}/CovidMode/top_scoring_lineages.fasta"
            else:
                return RESULT_DIR / "{sample}/FetchData/concat_strain_genomes.fasta"
