
import click
import numpy as np
from Bio import SeqIO


@click.group()
def cli():
    pass


def print_stats(input_fasta: str) -> None:
    with open(input_fasta) as handle:
        count = 0
        avg = 0
        min = np.Inf
        max = np.NINF
        for rec in SeqIO.parse(handle, "fasta"):
            l = len(rec)
            count += 1
            avg += l
            
            if l > max:
                max = l
            if l < min:
                min = l

        avg /= count

        print(f"Num sequences: {count}")
        print(f"Min length: {min}")
        print(f"Max length: {max}")
        print(f"Avg length: {avg}")
            


@cli.command()
@click.argument('input_fasta', type=click.Path(exists=True))
def stats(input_fasta: str) -> None:
    """
    Prints short statistics about the file
    """
    print_stats(input_fasta)


if __name__ == "__main__":
    cli()
