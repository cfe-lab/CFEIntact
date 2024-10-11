import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Iterator, List, Dict
import contextlib
import importlib.resources as resources


@contextlib.contextmanager
def reference_dir() -> Iterator[Path]:
    """
    A context manager handling subtype sequences paths packaged with CFEIntact.

    The complexity of the function arises from the need to maintain compatibility with
    multiple python versions due to changes in APIs of the `importlib.resources` package.

    It first tries to fetch the resource using `resources.files` function introduced in
    Python 3.9. If it fails, it falls back on `resources.path`.
    It further ensures that the obtained resource is returned
    as a Path instance regardless of it being a string, Path, or contextlib context-manager instance.

    Note: `resources.path` is set to be deprecated in future Python versions, hence the
    intended primary method is using `resources.files`.

    Yields:
        Path: A path-like object pointing to the 'subtype_alignments' resource within 'cfeintact'.
    """

    try:
        ret = resources.as_file(resources.files('cfeintact').joinpath('subtype_alignments'))
    except AttributeError:
        ret = resources.path('cfeintact', 'subtype_alignments')

    if isinstance(ret, str):
        yield Path(ret)
    elif isinstance(ret, Path):
        yield ret
    else:
        with ret as path:
            yield path


# This variable holds a mapping between subtype names and their sequences.
# Prevents re-reading the same sequence from disk multiple time.
SEQUENCE_CACHE: Dict[str, List[SeqRecord]] = {}


def HXB2() -> SeqRecord:
    """
    Return the sequence of HXB2, the standard HIV reference sequence
    """
    return subtype_sequences("HXB2")[0]


def subtypes() -> List[str]:
    """
List all currently available HIV subtypes
"""
    with reference_dir() as REFERENCE_DIR:
        return [f.replace(".fasta", "") for f in os.listdir(REFERENCE_DIR) if f.endswith(".fasta")]


@contextlib.contextmanager
def alignment_file(subtype: str) -> Iterator[Path]:
    with reference_dir() as REFERENCE_DIR:
        yield Path(os.path.join(REFERENCE_DIR, subtype + ".fasta"))


def subtype_sequences_iter(subtype: str) -> Iterator[SeqRecord]:
    with alignment_file(subtype) as path:
        yield from SeqIO.parse(path, "fasta")


def subtype_sequences(subtype: str) -> List[SeqRecord]:
    """
Return a list of sequences associated with an HIV subtype.

Args:
    subtype: folder in which to put temporary files.
    """

    if subtype not in SEQUENCE_CACHE:
        SEQUENCE_CACHE[subtype] = list(subtype_sequences_iter(subtype))

    return SEQUENCE_CACHE[subtype]
