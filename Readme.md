# Shallow CNOT circuit synthesis for constrained architecture

This repository contains an implementation of the code presented in [this work](https://arxiv.org/abs/2303.07302).

## Installing

Start by cloning the repository:

```bash
git clone https://github.com/smartiel/block_kutin.git
```

Everything is handled by pip:

```bash
pip install .
```

You might have [to install Rust](https://www.rust-lang.org/tools/install).


## The library

### Enumerating block operators

The synthesis algorithm relies on an exhaustive pre-enumeration of shallow circuit simplementing a set of block operators.
This enumeration procedure needs to be run once per given block architecture.

It is written in Rust and binded in python is a nice wrapper.
The following code runs the exhaustive enumeration and stores the resulting database on disk for later usage:

```python
from block_kutin import Database

database = Database("ladder_2")
# ladder_2 layout has a block size of k=2
# Operators in the database have shape (2k, k) for step 1, (2k, 2k) for step 2, (k, k) for step 3

import numpy as np

# This is a valid Step 1 operator (it is full rank)
array = np.array([[1, 0], [0, 1], [0, 1], [1, 0]])
# circuit is a list of pairs of integers (control, target)
circuit = database.get_circuit(array, step=1)
print(circuit)


# This is a valid Step 2 operator (it is full rank and its bottom right block is 0)
array = np.array([[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
# circuit is a list of pairs of integers (control, target)
circuit = database.get_circuit(array, step=2)
print(circuit)


# This is a valid Step 3 operator (it is full rank )
array = np.array([[1, 0], [1, 1]])
# circuit is a list of pairs of integers (control, target)
circuit = database.get_circuit(array, step=3)
print(circuit)

```

### Full synthesis algorithm

The synthesis algorithm can be called as follows:
```python
from block_kutin import synthesis
import numpy as np
from block_kutin import synthesis, depth, upper_bound_depth

topology = "ladder_3"
block_size = 3
n = 15  # must be a multiple of block_size

# generating a random operator
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
```

In this code:
- `synthesis` performs the synthesis and returns a list of pairs (control, target)
- `depth` is a helper function that computes the depth of the resulting circuit if needed
- `upper_bound_depth` computes the depth upper bound as derived in the paper