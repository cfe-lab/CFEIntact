from typing import Iterable, TypeVar

T = TypeVar("T")


def get_indel_impact(reference: Iterable[T], query: Iterable[T]) -> int:
    shift = 0
    impacted = 0

    for (x, y) in zip(reference, query):
        if x == "-" and y != "-":
            shift += 1
        if x != "-" and y == "-":
            shift -= 1

        if shift % 3 != 0:
            impacted += 1

    return impacted
