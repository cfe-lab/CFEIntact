
def map_positions(reference, query):
    """
    Returns a mapping that for each position in the reference sequence,
    returns the correspoding position in the query sequence.

    Args:
        reference: refererence sequence (aligned to query)
        query: query sequence (aligned to reference)
    """

    assert len(reference) == len(query)

    reference_pos = 0
    query_pos = 0
    reference_len = len(reference) - reference.count("-")
    mapping = [None] * reference_len
    last_pos = 0

    for i in range(len(reference)):
        if query_pos >= reference_len:
            break
        if query_pos == last_pos:
            last_pos += 1
            mapping[query_pos] = reference_pos
        if reference[i] != "-":
            query_pos += 1
        else:
            if query[i] != "-":
                reference_pos += 1
                mapping[query_pos] = reference_pos
                continue

        if query[i] != "-":
            reference_pos += 1

    return mapping
