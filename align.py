
import click
from Bio import AlignIO
from Bio import SeqIO


@click.group()
def cli():
    pass

# https://www.biostars.org/p/326044/
def column(alignment):
    '''Generator for getting list of characters in each column in the alignment'''
    for i in range(alignment.get_alignment_length()):
        c = []

        for record in alignment:
            c.append(record.seq[i])

        yield c

def is_gap(base):
    return base.lower() in ['-', 'n']

def is_invariant(col):
    filtered = [val for val in col if not is_gap(val)]
    return len(set(filtered)) < 2


def print_stats(input_fasta: str) -> None:
    with open(input_fasta) as handle:
        alignment = AlignIO.read(handle, "fasta")
        print("Input file:", input_fasta)
        print(f"Alignment size: {alignment.get_alignment_length()}bp x {len(alignment)}")

        g = sum([sum(is_gap(base) for base in x.seq) for x in alignment])
        l = sum([len(x.seq) for x in alignment])
        
        print(f"Gaps: {(g / l * 100):.2f}%")
        print(f"Mean alignment length (no gaps): {(l - g) / len(alignment)}")

        invariant = [col for col in column(alignment) if is_invariant(col)]
        print(f"Invariant sites: {(len(invariant) / alignment.get_alignment_length() * 100):.2f}%")


def cut_alignment(input_fasta: str, output_file: str, begin: int, end: int) -> None:
    with open(input_fasta) as handle:
        alignment = AlignIO.read(handle, "fasta")

        length = alignment.get_alignment_length()
        new_length = end - begin

        if new_length < 0:
            raise RuntimeError("End position must be higher than the begin position")

        if begin + new_length > length:
            raise RuntimeError("Alignment is too short")

        new_records = []
        for record in alignment:
            new_records.append(record[begin:end])

        new_alignment = AlignIO.MultipleSeqAlignment(new_records)
        AlignIO.write(new_alignment, output_file, "fasta")


@cli.command()
@click.argument('input_fasta', type=click.Path(exists=True))
def stats(input_fasta: str) -> None:
    """
    Prints short statistics about the alignment
    """
    print_stats(input_fasta)


@cli.command()
@click.argument('input_fasta', type=click.Path(exists=True))
@click.option('-o', '--output-file', type=click.Path(exists=False), required=True)
@click.option('-b', '--begin', type=int, required=True)
@click.option('-e', '--end', type=int, required=True)
def cut(input_fasta: str, output_file: str, begin: int, end: int) -> None:
    """
    Takes the slice [begin:end] of the input alignment and saves it to a file
    """
    cut_alignment(input_fasta, output_file, begin, end)


@cli.command()
@click.argument('input_fasta', type=click.Path(exists=True))
@click.argument('output_fasta')
def ungap(input_fasta: str, output_fasta: str) -> None:
    sequences = [SeqIO.SeqRecord(record.seq.ungap("-"), id=record.id, name=record.name, description=record.description)
                    for record in SeqIO.parse(input_fasta, "fasta")]
    SeqIO.write(sequences, output_fasta, "fasta")


if __name__ == "__main__":
    cli()
