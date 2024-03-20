from dataclasses import dataclass
from typing import Any


@dataclass
class ReferenceIndex:
    value: int

    def mapto(self, aligned: Any) -> int:
        ret: int = aligned.map_index(self.value)
        return ret
