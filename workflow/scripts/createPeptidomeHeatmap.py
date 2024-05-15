import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def getPsmDF(csv):
    """Preprocesses the report of the last search step. 
    The protein list is splitted and each PSM corresponds to one row.

    Args:
        csv (string): Path to the report of the last search step.

    Returns:
        pandas DataFrame: DataFrame where each row represents one PSM.
    """
    psm_report_df = pd.read_csv(csv, sep="\t", index_col=None)

    psm_report_df["Protein(s)"] = psm_report_df["Protein(s)"].astype(str)
    psm_report_df["Protein(s)"] = psm_report_df["Protein(s)"].apply(lambda x: x.split(","))
    psm_report_df = psm_report_df.explode("Protein(s)", ignore_index=True)
    psm_report_df = psm_report_df.rename(columns={"Protein(s)": "Accessions"})
    psm_report_df["Accessions"] = psm_report_df["Accessions"].str.split(".").str[0]

    psm_report_df = psm_report_df.rename(columns={psm_report_df.columns[0]: "psm_id"})
    psm_report_df = psm_report_df[["Accessions", "psm_id"]]

    return psm_report_df

def getAccessionsToCompare(df):
    """Gets all accession in the pandas DataFrame.

    Args:
        df (pandas DataFrame): preprocessed DataFrame from getPsmDF().

    Returns:
        list: List of all accession in the report.
    """
    accessions_to_compare = df["Accessions"].drop_duplicates().to_list()
    return accessions_to_compare


def computeSimilarityMatrix(df, accessions_to_compare):
    """Computes the pairwise similarity matrix of the identified peptidomes.

    Args:
        df (pandas DataFrame): DataFrame containing 
        accessions_to_compare (_type_): _description_

    Returns:
        numpy array: Matrix containing the pairwise similarities of the identified peptidomes.
    """
    similarity_matrix = np.zeros((len(accessions_to_compare), len(accessions_to_compare)))
    for accession1 in accessions_to_compare:
        list1  = df[df["Accessions"] == accession1]["psm_id"].to_list()
        for accession2 in accessions_to_compare[accessions_to_compare.index(accession1) + 1:]:
            list2  = df[df["Accessions"] == accession2]["psm_id"].to_list()
            count = len(set(list1) & set(list2)) # number of identical PSMs
            similarity = count / max(len(list1), len(list2))
            similarity_matrix[accessions_to_compare.index(accession1), accessions_to_compare.index(accession2)] = similarity
            similarity_matrix[accessions_to_compare.index(accession2), accessions_to_compare.index(accession1)] = similarity
            logging.info(f"Similarity of identified peptidomes of {accession1} and {accession2} = {similarity}")
    np.fill_diagonal(similarity_matrix, 1)

    return similarity_matrix


def createHeatmap(similarity_matrix, accessions_to_compare, output):
    """Creates a heatmap visualizing the pairwise similarities of the idenified peptidomes. 

    Args:
        similarity_matrix (numpy array): Matrix containing the pairwise similarities of the identified peptidomes.
        accessions_to_compare (list): List of accessions in the report.
        output (string): Path where the heatmap is stored as a .png.
    """

    plt.figure(figsize=(len(similarity_matrix), round((len(similarity_matrix) / 5) * 4))) # 5:4 ratio
    sns.heatmap(similarity_matrix, annot=True, cmap='coolwarm', fmt=".2f")

    plt.title('Pairwise Similarities Heatmap')
    plt.xlabel('Accessions')
    plt.ylabel('Accessions')

    tick_labels = accessions_to_compare 
    plt.xticks(np.arange(len(tick_labels)) + 0.5, tick_labels, rotation=45)
    plt.yticks(np.arange(len(tick_labels)) + 0.5, tick_labels, rotation=0)

    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    psm_report = snakemake.input[0]
    heatmap_png = snakemake.output[0]
    std_out = snakemake.log[1]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    psm_report_df = getPsmDF(psm_report)
    accessions_to_compare = getAccessionsToCompare(psm_report_df)
    similarity_matrix = computeSimilarityMatrix(psm_report_df, accessions_to_compare)
    createHeatmap(similarity_matrix, accessions_to_compare, heatmap_png)


if __name__ == '__main__':
    main()