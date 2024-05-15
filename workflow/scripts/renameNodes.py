import logging
import pandas as pd
from ete3 import Tree


def renameNodes(newick_file, strain_accessions_df):
    """Renames the nodes in the newick file with the taxon ids.

    Args:
        newick_file (string): path to the newick file
        strain_accessions_df (pandas DataFrame): pandas DataFrame containing the strain accessions
        renamed_newick (string): path to the newick file with the renamed nodes
    
    Returns:
        ete3.Tree: phylogenetic tree with the renamed nodes
    """

    tree = Tree(newick_file)

    for node in tree.traverse():
        if node.is_leaf():
            taxon_id = strain_accessions_df.loc[strain_accessions_df["genbank_accession"] == node.name, "taxa"].values[0]
            node.name = taxon_id
    
    return tree


def writeNewick(phylo_tree, out_file):
    """Writes the newick file.

    Args:
        phylo_tree (ete3.Tree): phylogenetic tree
        out_file (string): path to the newick file with the renamed nodes
    """

    phylo_tree.write(outfile=out_file)


def main():
    newick_file = snakemake.input[0]
    strain_accessions = snakemake.input[1]
    renamed_newick = snakemake.output[0]
    std_out = snakemake.log[0]

    logging.basicConfig(filename=std_out, level=logging.INFO)

    strain_accessions_df = pd.read_csv(strain_accessions, sep=",")
    tree = renameNodes(newick_file, strain_accessions_df)
    writeNewick(tree, renamed_newick)


if __name__ == "__main__":
    main()