
import os
import sys
import subprocess
import tempfile
import csv
from dataclasses import dataclass

import util.subtypes as st


@dataclass
class BlastRow:
    qseqid: str
    qlen: int
    sseqid: str
    sgi: str
    slen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    length: int
    pident: float
    nident: float
    btop: int
    stitle: str
    sstrand: str


def blast(subtype, input_file, output_file):
    db_file = st.alignment_file(subtype)

    subprocess.call(
        ["blastn",
         "-query", input_file,
         "-db", db_file,
         "-num_alignments", "1",
         "-reward", "1",
         "-penalty", "-1",
         "-gapopen", "2",
         "-gapextend", "1",
         "-out", output_file,
         "-outfmt", "6 qseqid qlen sseqid sgi slen qstart qend sstart send evalue bitscore length pident nident btop stitle sstrand"],
        shell=False)


def init_blast_row(row):
    it = iter(row)
    return BlastRow(
        qseqid=next(it),
        qlen=int(next(it)),
        sseqid=next(it),
        sgi=next(it),
        slen=int(next(it)),
        qstart=int(next(it)),
        qend=int(next(it)),
        sstart=int(next(it)),
        send=int(next(it)),
        evalue=float(next(it)),
        bitscore=float(next(it)),
        length=int(next(it)),
        pident=float(next(it)),
        nident=float(next(it)),
        btop=next(it),
        stitle=next(it),
        sstrand=next(it),
    )


def iterate_values_from_tsv(file_path):
    with open(file_path, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        for row in reader:
            yield row


def iterate_blast_rows_from_tsv(file_path):
    previous_key = None
    values = []

    for row in iterate_values_from_tsv(file_path):
        key = row[0]
        typed = init_blast_row(row)

        if key != previous_key and previous_key is not None:
            yield values
            values = []

        values.append(typed)
        previous_key = key

    if values:
        yield values


def blast_interate(subtype, input_file):
    with tempfile.NamedTemporaryFile() as output_file:
        blast(subtype, input_file, output_file.name)
        for seq in iterate_blast_rows_from_tsv(output_file.name):
            yield seq


def is_sorted(lst):
    last = None
    for x in lst:
        if last is not None and x < last:
            return False
        last = x
    return True


def check_scramble(blast_rows):
    # HIV 5' region can easily map to its 3' region because they are identical.
    # Such a maping would not constitute a scramble, so we ignore the 5' region for this check.
    ignored_5_prime = [x for x in blast_rows if x.sstart > 622 and x.send > 622]

    if not ignored_5_prime:
        # No alignment.
        # It should be an error normally, yet not a scramble error.
        return None

    all_same = len(set(x.sstrand for x in ignored_5_prime)) == 1
    if not all_same:
        # Some parts of the sequence were aligned
        # in forward direction (plus)
        # and some in reverse (minus).
        # This indicates an internal inversion.
        return "mix"

    ignored_5_prime.sort(key=lambda x: x.qstart)
    direction = ignored_5_prime[0].sstrand
    if direction == "plus" and is_sorted(x.sstart for x in ignored_5_prime):
        return None
    elif direction == "minus" and is_sorted(x.send for x in reversed(ignored_5_prime)):
        return None
    else:
        return direction + "Scramble"


def check_nonhiv(blast_rows):
    aligned_length = sum(abs(x.qend - x.qstart) + 1 for x in blast_rows)
    total_length = blast_rows[0].qlen if blast_rows else 1
    ratio = aligned_length / total_length

    if ratio < 0.8:
        return "CHIMERA"
    else:
        return None
