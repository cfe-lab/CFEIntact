from Bio import Seq


def translate_to_aminoacids(seq: str, frame: int = 0, to_stop: bool = False) -> str:
    for_translation = seq[frame:]
    for_translation += 'N' * ({0: 0, 1: 2, 2: 1}[len(for_translation) % 3])
    ret: str = Seq.translate(for_translation, to_stop=to_stop)
    return ret
