import numpy as np
from block_kutin import synthesis, depth, upper_bound_depth

topology = "ladder_3"
block_size = 3
n = 15  # must be a multiple of block_size

operator = np.eye(n, dtype=np.byte)
for _ in range(100):
    i, j = np.random.choice(list(range(n)), size=2, replace=False)
    operator[i] ^= operator[j]

circuit = synthesis(operator, block_size, topology)
print(
    "cnot depth:",
    depth(circuit),
    "| upper bound:",
    upper_bound_depth(n, topology, block_size),
)
