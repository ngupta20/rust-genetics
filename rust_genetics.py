from itertools import combinations

import numpy as np

POSSIBLE = tuple("YGHXW")
GREEN = tuple("YGH")
RED = tuple("XW")
WEIGHT = [0.6, 0.6, 0.6, 1.0, 1.0]


def recommend(file: str = "genes.csv"):
    with open(file, "r") as file:
        gene_pool = np.array(tuple(map(tuple, file.read().split("\n"))))
    cs = [
        gene_pool[np.array(tuple(combinations(range(len(gene_pool)), i)))]
        for i in range(1, 9)
    ]
    ws = [(cs[i][:, :, :, None] == POSSIBLE).sum(axis=1) * WEIGHT for i in range(8)]
    vs = [
        cs[i][
            (ws[i].max(axis=2)[:, :, None] == ws[i])[:, :, :3].any(axis=2).all(axis=1)
        ]
        for i in range(8)
    ]
    return sum([[["".join(z) for z in y] for y in x] for x in vs if x.shape[0] > 0], [])
