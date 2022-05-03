from itertools import combinations
from typing import List, Tuple

import numpy as np

WEIGHT = [0.6, 0.6, 0.6, 1.0, 1.0]
POSSIBLE = tuple("YGHXW")
GREEN = tuple("YGH")
RED = tuple("XW")


def recommendations(gene_pool: List[str]) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    pool = np.array(tuple(map(tuple, gene_pool)))
    n = min(len(pool) + 1, 9)

    viable_combinations = []
    viable_outcomes = []
    for i in range(n - 1):
        index_combinations = np.array(tuple(combinations(range(len(pool)), i + 1)))
        possible = pool[index_combinations]
        weights = (possible[:, :, :, None] == POSSIBLE).sum(axis=1) * WEIGHT
        outcomes = (weights.max(axis=2)[:, :, None] == weights)[:, :, :2]
        viable = (
            ((outcomes[:, :, [0]] ^ outcomes[:, :, [1]]) & outcomes[:, :, :2]).sum(
                axis=1
            )
            >= 3
        ).all(axis=1)
        viable_combinations += [possible[viable]]
        viable_outcomes += [outcomes[viable]]

    return viable_combinations, viable_outcomes
