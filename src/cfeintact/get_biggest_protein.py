from typing import List, Tuple, Optional

def get_biggest_protein(has_start_codon: bool, aminoacids: str) -> Optional[Tuple[str, int, int]]:
    def skip_to_startcodon(x):
        index = x.find("M")
        if index >= 0:
            return x[index:]
        else:
            return ""

    candidates: List[Tuple[str, int, int]] = []
    current_protein = ""

    if has_start_codon:
        skipping_to_startcodon = True
        taking_until_stopcodon = False
    else:
        skipping_to_startcodon = False
        taking_until_stopcodon = True

    for i, aa in enumerate(aminoacids):
        if skipping_to_startcodon:
            if aa == "M":
                skipping_to_startcodon = False
                taking_until_stopcodon = True
            else:
                continue

        if taking_until_stopcodon:
            if aa == "*":
                if current_protein:
                    start = i - len(current_protein)
                    end = i - 1
                    candidates.append((current_protein, start, end))

                current_protein = ""
                if has_start_codon:
                    skipping_to_startcodon = True
                    taking_until_stopcodon = False
            else:
                current_protein += aa

    if current_protein:
        i += 1
        start = i - len(current_protein)
        end = i - 1
        candidates.append((current_protein, start, end))

    if candidates:
        return max(candidates, key=lambda p: len(p[0]))
    else:
        return None
