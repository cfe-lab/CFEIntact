from dataclasses import dataclass

@dataclass
class ReferenceIndex:
    value: int

    def mapto(self, aligned):
        return aligned.map_index(self.value)
