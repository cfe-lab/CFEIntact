from typing import Iterable, TypeVar

T = TypeVar("T")


def get_indel_impact(reference: Iterable[T], query: Iterable[T]) -> float:
    shift = 0
    impacted = 0
    counter = 0
    good = 0

    for (x, y) in zip(reference, query):
        if x == "-" and y != "-":
            shift += 1
        if x != "-" and y == "-":
            shift -= 1
        if shift % 3 != 0:
            counter += 1
            good = 0
        else:
            good += 1
            if good > 5:
                impacted += counter if counter > 30 else 0
                counter = 0

    impacted += counter if counter > 30 else 0
    return impacted
