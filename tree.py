import click
from Bio import Phylo


@click.group()
def cli():
    pass


def reroot_tree(input_newick: str, output_newick: str) -> None:
    tree = Phylo.read(input_newick, "newick")
    tree.root_at_midpoint()
    Phylo.write(tree, output_newick, "newick")

@cli.command()
@click.argument('input_newick', type=click.Path(exists=True))
@click.argument('output_newick', type=click.Path(exists=False))
def reroot(input_newick: str, output_newick: str) -> None:
    """
    Reroots the input newick-formatted tree at midpoint
    """
    reroot_tree(input_newick, output_newick)


if __name__ == "__main__":
    cli()
