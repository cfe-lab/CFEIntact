
def get_biggest_protein(has_start_codon, aminoacids):
    def skip_to_startcodon(x):
        index = x.find("M")
        if index >= 0:
            return x[index:]
        else:
            return ""

    parts = aminoacids.split("*")
    subparts = [skip_to_startcodon(x)
                for x in parts] if has_start_codon else parts
    longest = max(subparts, key=len)
    return longest
