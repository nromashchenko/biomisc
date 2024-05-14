import click
import gzip
import numpy as np
from Bio import SeqIO


@click.group()
def cli():
    pass


def print_stats(input_fasta: str) -> None:
    with gzip.open(input_fasta, "rt") if input_fasta.endswith('.gz') else open(input_fasta) as handle:
        count = 0
        avg = 0
        min = np.Inf
        min_name = None
        max = np.NINF
        max_name = None

        for rec in SeqIO.parse(handle, "fasta"):
            l = len(rec)
            count += 1
            avg += l
            
            if l > max:
                max = l
                max_name = rec.id
            if l < min:
                min = l
                min_name = rec.id

        avg /= count

        print(f"Num sequences: {count}")
        print(f"Min length: {min} {min_name}")
        print(f"Max length: {max} {max_name}")
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
