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
