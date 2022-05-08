from itertools import combinations
from typing import List, Tuple

import numpy as np

WEIGHT = [0.6, 0.6, 0.6, 1.0, 1.0]
POSSIBLE = tuple("YGHXW")
GREEN = tuple("YGH")
RED = tuple("XW")


def recommendations(
    gene_pool: List[str], viability_method="yg"
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    pool = np.array(tuple(map(tuple, gene_pool)))
    pool = pool[(pool[:, :, None] == RED).any(axis=2).sum(axis=1) <= 2]
    n = min(len(pool) + 1, 9)
    viable_combinations = []
    viable_outcomes = []
    for i in range(n - 1):
        index_combinations = np.array(tuple(combinations(range(len(pool)), i + 1)))
        possible = pool[index_combinations]
        weights = (possible[:, :, :, None] == POSSIBLE).sum(axis=1) * WEIGHT
        outcomes = weights.max(axis=2)[:, :, None] == weights
        if viability_method == "yg":
            viable = outcomes[:, :, :2].any(axis=2).sum(axis=1) >= 6
        elif viability_method == "all_green":
            viable = outcomes[:, :, :3].any(axis=2).sum(axis=1) >= 6
        elif viability_method == "yg-1":
            viable = outcomes[:, :, :3].any(axis=2).sum(axis=1) >= 5
        else:
            raise NotImplementedError()
        print(f"Combination Size: {i + 1}    Viable: {viable.sum()}")
        viable_combinations += [possible[viable]]
        viable_outcomes += [outcomes[viable]]

    return viable_combinations, viable_outcomes
