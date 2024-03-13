import dataclasses
from dataclasses import dataclass
from Bio import Seq, SeqRecord

import cfeintact.coordinates as coords
import cfeintact.wrappers as wrappers
from cfeintact.mapped_orf import MappedORF
from cfeintact.original_orf import OriginalORF
from cfeintact.reference_index import ReferenceIndex
from cfeintact.find_orf import find_orf

@dataclass
class AlignedSequence:
    this: Seq
    reference: Seq
    orfs: dict[str, MappedORF]     = dataclasses.field(default=None)
    alignment: (str, str)          = dataclasses.field(default=None)
    coordinates_mapping: list[int] = dataclasses.field(default=None)


    def get_alignment(self):
        if not self.alignment:
            self.alignment = wrappers.mafft([self.reference, self.this])

        return self.alignment


    def aligned_reference(self):
        return self.get_alignment()[0]


    def aligned_this(self):
        return self.get_alignment()[1]


    def get_coordinates_mapping(self):
        if not self.coordinates_mapping:
            self.coordinates_mapping = coords.map_positions(self.aligned_reference(), self.aligned_this())

        return self.coordinates_mapping


    def map_index(self, index):
        if isinstance(index, ReferenceIndex):
            index = index.value

        if not isinstance(index, int):
            raise TypeError(f"Expected integer as index", index)

        mapping = self.get_coordinates_mapping()
        if index < len(mapping):
            return mapping[index]
        else:
            return mapping[-1]


    def index(self, index):
        if isinstance(index, ReferenceIndex):
            index = self.map_index(index)

        return self.this[index]


    def slice(self, first, last):
        if isinstance(first, ReferenceIndex):
            first = self.map_index(first)
        if isinstance(last, ReferenceIndex):
            last = self.map_index(last)

        newthis = self.this[first:(last + 1)]
        newreference = self.reference[self.map_index(first):(self.map_index(last) + 1)]
        # TODO: calculate new "coordinates_mapping" and new "alignment" from these indexes.
        return AlignedSequence(this=newthis, reference=newreference)


    def reverse(self):
        newthis = SeqRecord.SeqRecord(Seq.reverse_complement(self.this.seq),
                                      id = self.this.id + "[REVERSE_COMPLEMENT]",
                                      name = self.this.name
                                      )

        return AlignedSequence(this=newthis, reference=self.reference)


    def alignment_score(self):
        return sum([a==b for a, b in zip(self.get_alignment()[0].seq, self.get_alignment()[1].seq)])


    def get_orf(self, expected_orf: OriginalORF):
        if self.orfs is None:
            self.orfs = {}

        if expected_orf.name not in self.orfs:
            self.orfs[expected_orf.name] = find_orf(self, expected_orf)

        return self.orfs[expected_orf.name]
