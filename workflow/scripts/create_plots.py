import logging
import pandas as pd
import matplotlib.pyplot as plt


def CreateStrainCountBarPlot(in_file, out_file):
    """Creates a bar plot of counts vs strains.

    Args:
        in_file (string): path to the strain mapping file
        out_file (string): path where the bar plot is stored
    """

    df = pd.read_csv(in_file, delimiter="\t")

    df["strain"] = df["strain"].astype(str)
    df["isolate"] = df["isolate"].astype(str)
    df["strain_isolate"] = df["strain"].replace("nan", "") + "_" + df["isolate"].replace("nan", "")
    df["strain_isolate"] = df["strain_isolate"].str.strip("_")
    df['strain_isolate'] = df["strain_isolate"].replace("", "NaN")
    sorted_df = df.sort_values(by="counts", ascending=False)

    counts = sorted_df["counts"]
    strain_isolate = sorted_df["strain_isolate"]

    fig, ax = plt.subplots()
    ax.bar(strain_isolate, counts)

    # Add labels and title
    ax.set_xlabel("Strain") 
    ax.set_ylabel("Counts")
    ax.set_title("Bar Plot of Counts vs Strain/Isolate")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateStrainCountConfidenceScoresBarPlot(in_file, out_file):
    """Creates a bar plot of confidence scores vs strains.

    Args:
        in_file (string): path to the strain mapping file
        out_file (string): path where the bar plot is stored
    """
        
    df = pd.read_csv(in_file, delimiter="\t")

    df["strain"] = df["strain"].astype(str)
    df["isolate"] = df["isolate"].astype(str)
    df["strain_isolate"] = df["strain"].replace("nan", "") + "_" + df["isolate"].replace("nan", "")
    df["strain_isolate"] = df["strain_isolate"].str.strip("_")
    df['strain_isolate'] = df["strain_isolate"].replace("", "NaN")
    sorted_df = df.sort_values(by="confidence_scoring", ascending=False)

    confidence_scoring = sorted_df["confidence_scoring"]
    strain_isolate = sorted_df["strain_isolate"]

    fig, ax = plt.subplots()
    ax.bar(strain_isolate, confidence_scoring)

    # Add labels and title
    ax.set_xlabel("Strain") 
    ax.set_ylabel("Confidence scoring")
    ax.set_title("Bar Plot of confidence scoring vs Strain/Isolate")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateStrainConfidenceBarPlot(in_file, out_file):
    """create bar plot of mean confidences vs strains

    Args:
        in_file (string): path to the strain mapping file
        out_file (string): path where the bar plot is stored
    """

    df = pd.read_csv(in_file, delimiter="\t")

    df["strain"] = df["strain"].astype(str)
    df["isolate"] = df["isolate"].astype(str)
    df["strain_isolate"] = df["strain"].replace("nan", "") + "_" + df["isolate"].replace("nan", "")
    df["strain_isolate"] = df["strain_isolate"].str.strip("_")
    df['strain_isolate'] = df["strain_isolate"].replace("", "NaN")
    sorted_df = df.sort_values(by="mean_confidence", ascending=False)

    mean_confidence = sorted_df["mean_confidence"]
    strain_isolate = sorted_df["strain_isolate"]

    fig, ax = plt.subplots()
    ax.bar(strain_isolate, mean_confidence)

    # Add labels and title
    ax.set_xlabel("Strain") 
    ax.set_ylabel("Mean Confidence")
    ax.set_title("Bar Plot of Mean Confidence vs Strain/Isolate")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateStrainProteomeLengthScoringBarPlot(in_file, out_file):
    """create bar plot of the proteome length scoring vs strains

    Args:
        in_file (string): path to the strain mapping file
        out_file (string): path where the bar plot is stored
    """

    df = pd.read_csv(in_file, delimiter="\t")

    df["strain"] = df["strain"].astype(str)
    df["isolate"] = df["isolate"].astype(str)
    df["strain_isolate"] = df["strain"].replace("nan", "") + "_" + df["isolate"].replace("nan", "")
    df["strain_isolate"] = df["strain_isolate"].str.strip("_")
    df['strain_isolate'] = df["strain_isolate"].replace("", "NaN")
    sorted_df = df.sort_values(by="length_scoring", ascending=False)

    length_scoring = sorted_df["length_scoring"]
    strain_isolate = sorted_df["strain_isolate"]

    fig, ax = plt.subplots()
    ax.bar(strain_isolate, length_scoring)

    # Add labels and title
    ax.set_xlabel("Strain") 
    ax.set_ylabel("Proteome Length Scoring")
    ax.set_title("Bar Plot of Proteome Length Scoring vs Strain/Isolate")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateTaxIdScoresBarPlot(in_file, out_file):
    """create bar plot of scores (weights) vs species taxon ids

    Args:
        in_file (string): path to the taxid scores file
        out_file (string): path where the bar plot is stored
    """

    df = pd.read_csv(in_file, sep="\t")

    sorted_df = df.sort_values(by="weight", ascending=False)

    taxids = sorted_df["taxid"].astype(str)[:30]
    weights = sorted_df["weight"][:30]

    plt.bar(taxids, weights)
    plt.xlabel("Tax ID")
    plt.ylabel("Weight")
    plt.title("Tax ID vs Weight")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()

def CreateORFScoresBarPlot(in_file, out_file):
    """create bar plot of scores (weights) vs strain taxon ids

    Args:
        in_file (string): path to the orf scores file
        out_file (string): path where the bar plot is stored
    """

    df = pd.read_csv(in_file, sep=",")
    sorted_df = df.sort_values(by="weight", ascending=False)

    strains = sorted_df["strain_isolate"].astype(str)
    weights = sorted_df["weight"]

    plt.bar(strains, weights)
    plt.xlabel("Strain")
    plt.ylabel("Weight")
    plt.title("Strain vs Weight")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateProportionsPieChart(first_search, final_search, out_file, host_filter_psm_report=None):
    """create pie chart to display the amount of matches from the first and the final search step

    Args:
        first_search (string): path to the PSM Report of the first search
        final_search (string): path to the PSM Report of the final search
        out_file (string): path where the pie chart is stored
    """
   
    df1 = pd.read_csv(first_search, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
    df1.columns = ["Peptide"]
    df1["Peptide"] = df1.Peptide.apply(lambda x: x.split(","))
    report_df_1 = df1.explode("Peptide", ignore_index=True)

    df2 = pd.read_csv(final_search, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
    df2.columns = ["ORF"]
    df2["ORF"] = df2.ORF.apply(lambda x: x.split(","))
    report_df_2 = df2.explode("ORF", ignore_index=True)

    if host_filter_psm_report:
        df3 = pd.read_csv(host_filter_psm_report, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
        df3.columns = ["Peptide"]
        df3["Peptide"] = df3.Peptide.apply(lambda x: x.split(","))
        report_df_3 = df3.explode("Peptide", ignore_index=True)

        num_entries = [len(report_df_3), len(report_df_1), len(report_df_2)]
        labels = ["Host Filtering", "First Search", "Final Search"]
    
    else:
        num_entries = [len(report_df_1), len(report_df_2)]
        labels = ["First Search", "Final Search"]

    def percent_and_amount(x):
        """Convert the amount into a percentage.

        Args:
            x (float): The percentage value.

        Returns:
            str: A string representing the percentage and the calculated amount.
        """

        amount = int(round(x/100.0 * sum(num_entries), 0))
        p = "{:.1f}%  ({:d})".format(x, amount)
        return p
    
    plt.pie(num_entries, labels=labels, autopct=percent_and_amount)
    plt.title("Number of PSMs")
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateConfidenceHistogram(in_file, out_file):
    """create a histogram of the confidences of the PSMs

    Args:
        in_file (string): path to the PSM Report
        out_file (string): path where the histogram is stored
    """

    df = pd.read_csv(in_file, sep="\t", on_bad_lines="skip")
    df.rename(columns={"Protein(s)": "Proteins"}, inplace=True)
    df["Proteins"] = df["Proteins"].astype(str)
    df["Proteins"] = df.Proteins.apply(lambda x: x.split(","))
    exploded_df = df.explode("Proteins", ignore_index=True)

    plt.hist(exploded_df["Confidence [%]"], bins=20)
    plt.xlabel("Confidence")
    plt.ylabel("Count")
    plt.title("Confidence Histogram")
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def computeMass(seq):

    amino_acid_masses = {
    "A": 71.08,
    "R": 156.19,
    "N": 114.10,
    "D": 115.09,
    "C": 103.14,
    "E": 129.12,
    "Q": 128.13,
    "G": 57.05,
    "H": 137.14,
    "I": 113.16,
    "L": 113.16,
    "K": 128.17,
    "M": 131.19,
    "F": 147.18,
    "P": 97.12,
    "S": 87.08,
    "T": 101.11,
    "W": 186.21,
    "Y": 163.18,
    "V": 99.13,
    }
    
    mass = 0
    for amino_acid in seq:
        mass += amino_acid_masses[amino_acid]
    return mass


def CreateSuitabilityBarPlot(search_steps, reports, out_file):
    """Creates a bar plot visualizing the suitabilities of the databases used in 
    the different search steps.

    Args:
        search_steps (list): A list containing the names of the different search steps.
        reports (list): A list containing the paths to the report of the different search steps.
        out_file (string): Path where the bar plot is stored.
    """

    db_suitabilities = []
    for report in reports:
        psm_report_df = pd.read_csv(report, sep="\t", index_col=None)

        psm_report_df["Protein(s)"] = psm_report_df["Protein(s)"].apply(lambda x: x.split(","))
        psm_report_df = psm_report_df.explode("Protein(s)", ignore_index=True)

        psm_report_df = psm_report_df.rename(columns={psm_report_df.columns[0]: "psm_id"})
        psm_report_df["computed_mass"] = psm_report_df["Sequence"].apply(computeMass)
        psm_report_df["computed_mass"] = round(psm_report_df["computed_mass"], 2)
        unique_masses = psm_report_df["computed_mass"].unique()

        for mass in unique_masses:
            current_mass_df = psm_report_df.loc[psm_report_df["computed_mass"] == mass]
            not_novor_df = current_mass_df.loc[~current_mass_df["Protein(s)"].str.startswith("NOVOR_PEPTIDE_")]
            novor_df = current_mass_df.loc[current_mass_df["Protein(s)"].str.startswith("NOVOR_PEPTIDE_")]
            if not_novor_df.shape[0] > 0:        
                psm_report_df.drop(novor_df.index, inplace=True)

        novor_hits = psm_report_df.loc[psm_report_df["Protein(s)"].str.startswith("NOVOR_PEPTIDE_")].shape[0]
        total_hits = psm_report_df.shape[0]
        db_suitability = round(((total_hits - novor_hits) / total_hits) * 100, 2)

        if db_suitability < 45.0:
            logging.info(f"The database suitability for the step: {search_steps[len(db_suitabilities)]} is below 45%. Please check the database and the quality of the sample!")

        db_suitabilities.append(db_suitability)
   
    colors = ['orange' if value < 45.0 else '#1f77b4' for value in db_suitabilities] # '#1f77b4' == standard blue of matplotlib
    barplot = plt.bar(search_steps, db_suitabilities, color=colors)

    plt.title('Database Suitabilities for the different Search Steps')
    plt.ylabel('Suitability')
    plt.xlabel('Search Steps')

    # print the values on top of the bars
    for i, search_step in enumerate(barplot):
        plt.text(search_step.get_x() + search_step.get_width()/2, db_suitabilities[i], round(db_suitabilities[i], 2), ha='center', va='bottom')
        
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def main():
    strain_mappings = snakemake.input[0]
    taxID_scores = snakemake.input[1]
    frist_search = snakemake.input[2]
    final_search = snakemake.input[3]
    ORF_scores = snakemake.input[4]

    sample_res_path = snakemake.params[0]
    sample_name = snakemake.params[1]
    host_filtering = snakemake.params[2]
    extra_search = snakemake.params[3]
    compute_db_suitability = snakemake.params[4]

    std_out = snakemake.log[1]
    logging.basicConfig(filename=std_out, level=logging.INFO)
    
    strain_counts_bar_plot = f"{sample_res_path}/{sample_name}_strain_counts_bar_plot.png"
    strain_conf_bar_plot = f"{sample_res_path}/{sample_name}_strain_conf_bar_plot.png"
    taxIdScores_bar_plot = f"{sample_res_path}/{sample_name}_taxIdScores_bar_plot.png"
    proportions_pie_chart = f"{sample_res_path}/{sample_name}_proportions_pie_chart.png"
    first_search_confidence_histogram = f"{sample_res_path}/{sample_name}_first_search_confidence_histogram.png"
    final_search_confidence_histogram = f"{sample_res_path}/{sample_name}_final_search_confidence_histogram.png"
    ORFScores_bar_plot = f"{sample_res_path}/{sample_name}_ORFScores_bar_plot.png"
    count_confidence_scores_bar_plot = f"{sample_res_path}/{sample_name}_count_confidence_scores_bar_plot.png"
    proteome_length_scoring_bar_plot = f"{sample_res_path}/{sample_name}_proteome_length_scoring_bar_plot.png"

    input_counter = 4
    if host_filtering:
        input_counter += 1
        host_filter_psm_report = snakemake.input[input_counter]
    else:
        host_filter_psm_report = None
    
    if extra_search:
        input_counter += 1
        extra_search_strain_mappings = snakemake.input[input_counter]
        input_counter += 1
        extra_search = snakemake.input[input_counter]
        input_counter += 1
        extra_search_ORF_scores = snakemake.input[input_counter]

        extra_search_strain_conf_bar_plot = f"{sample_res_path}/{sample_name}_extra_search_strain_conf_bar_plot.png"
        extra_search_confidence_histogram = f"{sample_res_path}/{sample_name}_extra_search_confidence_histogram.png"
        extra_search_ORFScores_bar_plot = f"{sample_res_path}/{sample_name}_extra_search_ORFScores_bar_plot.png"
        extra_search_count_confidence_scores_bar_plot = f"{sample_res_path}/{sample_name}_extra_search_count_confidence_scores_bar_plot.png"
        extra_search_strain_counts_bar_plot = f"{sample_res_path}/{sample_name}_extra_search_strain_counts_bar_plot.png"
        extra_proteome_length_scoring_bar_plot = f"{sample_res_path}/{sample_name}_extra_search_proteome_length_scoring_bar_plot.png"
    else:
        extra_search_strain_mappings = None
        extra_search = None
        extra_search_ORF_scores = None

    if compute_db_suitability:
        suitability_bar_plot = f"{sample_res_path}/{sample_name}_database_suitability_bar_plot.png"
        search_steps = []
        reports = []
        if host_filtering:
            input_counter += 1
            db_suitability_host_filtering_report = snakemake.input[input_counter]
            search_steps.append("Host Filtering")
            reports.append(db_suitability_host_filtering_report)

        input_counter += 1
        db_suitability_first_search_report = snakemake.input[input_counter]
        input_counter += 1
        db_suitability_final_search_report = snakemake.input[input_counter]

        search_steps += ["First Search", "Final Search"]
        reports += [db_suitability_first_search_report, db_suitability_final_search_report]

        if extra_search:
            input_counter += 1
            db_suitability_extra_search_report = snakemake.input[input_counter]
            search_steps.append("Top-Scoring Search")
            reports.append(db_suitability_extra_search_report)
        
        CreateSuitabilityBarPlot(search_steps, reports, suitability_bar_plot)


    CreateStrainCountBarPlot(strain_mappings, strain_counts_bar_plot)
    CreateStrainConfidenceBarPlot(strain_mappings, strain_conf_bar_plot)
    CreateTaxIdScoresBarPlot(taxID_scores, taxIdScores_bar_plot)
    CreateConfidenceHistogram(frist_search, first_search_confidence_histogram)
    CreateConfidenceHistogram(final_search, final_search_confidence_histogram)
    # new scorings
    CreateORFScoresBarPlot(ORF_scores, ORFScores_bar_plot)
    CreateStrainCountConfidenceScoresBarPlot(strain_mappings, count_confidence_scores_bar_plot)
    CreateStrainProteomeLengthScoringBarPlot(strain_mappings, proteome_length_scoring_bar_plot)
    # Proportions Pie Chart
    CreateProportionsPieChart(frist_search, final_search, proportions_pie_chart, host_filter_psm_report)

    if extra_search:
        CreateStrainConfidenceBarPlot(extra_search_strain_mappings, extra_search_strain_conf_bar_plot)
        CreateConfidenceHistogram(extra_search, extra_search_confidence_histogram)
        CreateORFScoresBarPlot(extra_search_ORF_scores, extra_search_ORFScores_bar_plot)
        CreateStrainCountConfidenceScoresBarPlot(extra_search_strain_mappings, extra_search_count_confidence_scores_bar_plot)
        CreateStrainCountBarPlot(extra_search_strain_mappings, extra_search_strain_counts_bar_plot)
        CreateStrainProteomeLengthScoringBarPlot(extra_search_strain_mappings, extra_proteome_length_scoring_bar_plot)


if __name__ == "__main__":
    main()