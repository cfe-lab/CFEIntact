from Bio.Seq import translate, Seq, MutableSeq
from typing import Union


def translate_to_aminoacids(seq: Union[str, Seq, MutableSeq],
                            frame: int = 0,
                            to_stop: bool = False
                            ) -> str:
    for_translation = seq[frame:]
    for_translation += 'N' * ({0: 0, 1: 2, 2: 1}[len(for_translation) % 3])
    ret: str = translate(for_translation, to_stop=to_stop)
    return ret
