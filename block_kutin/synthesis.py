import numpy as np
import networkx as nx
import copy
from block_kutin import Database

# Display


def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:" + fmt + "}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:" + str(col_maxes[i]) + fmt + "}").format(y), end="  ")
        print("")


# Preliminaries circuit and graph functions


def apply_circuit(A, circ):
    """Apply a circuit to an operator A.

    Inputs:
        A: the np.array of the input linear reversible
            operator
        circ: the circuit to apply, given as a list of
                tuples (control, target) for each CNOT.

    Output:
        the output operator after execution of the circuit.
    """

    for gate in circ:
        A[gate[1], :] ^= A[gate[0], :]

    return A


def shift(circ, p):
    """Shifts a n-qubit circuit such that it is applied on the qubits p+k, k=1..n.
    Input:
        circ: the circuit to be shifted
        p: the number of qubits to shift
    Output:
        the circuit shifted.
    """

    new_circ = []
    for gate in circ:
        new_circ.append(([x + p for x in gate]))
    return new_circ


def depth(circ):
    """Computes the depth of the circuit.
    Input:
        circ: the circuit
    Output:
        the depth of the circuit
    """

    slices = []
    flag = False
    j = 0

    for gate in circ:
        flag = False
        for i in range(len(slices) - 1, -1, -1):
            for qubit in slices[i]:
                for qubit_gate in gate:
                    if qubit_gate == qubit:
                        flag = True
                        break
                if flag:
                    break

            if flag:
                j = i + 1
                break
            if i == 0:
                j = 0

        if j >= len(slices):
            slices.append(copy.copy(gate))
        else:
            for qubit in gate:
                slices[j].append(qubit)

    return len(slices)


def is_hardware_compliant(circ, archi, p, n):
    """Check if the circuit is hardware compliant."""
    G = create_connectivity_graph(n, archi, p)
    for gate in circ:
        assert [
            gate[1],
            gate[0],
        ] in G.edges, f"The circuit is not hardware compliant:{[gate[1], gate[0]]}"


def check_circuit(A, circ, archi, p):
    """Checks if a given CNOT circuit implements a given
    operator and respects the hardware constraints.

    Inputs:
        A: the np.array of the linear reversible
           operator that we want to implement.
        circ: the CNOT circuit to check
        archi: the name of the connectivity
        p: the block size (only useful when archi="all_to_all")

    Raises an error if something went wrong.
    """
    B = np.eye(np.size(A, 0), dtype=bool)
    G = create_connectivity_graph(np.size(A, 0), archi, p)

    for gate in circ:
        assert [
            gate[1],
            gate[0],
        ] in G.edges, f"The circuit is not hardware compliant:{[gate[1], gate[0]]}"
        B[gate[1], :] ^= B[gate[0], :]

    assert (A == B).all(), "The circuit does not implement the desired operator."


def create_connectivity_graph(n, archi, p):
    G = nx.Graph()
    assert n % p == 0
    number_of_blocks = n // p

    if archi == "lnn":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)

    elif archi == "line_2":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)
            G.add_edge(i * p + 1, i * p + 2)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)

    elif archi == "ladder_2":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)

            G.add_edge(i * p, i * p + 2)
            G.add_edge(i * p + 1, i * p + 3)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)

    elif archi == "ladder_2_diagonals":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)

            G.add_edge(i * p, i * p + 2)
            G.add_edge(i * p + 1, i * p + 3)
            G.add_edge(i * p, i * p + 3)
            G.add_edge(i * p + 1, i * p + 2)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)

    elif archi == "ladder_3":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)
            G.add_edge(i * p + 1, i * p + 2)

            G.add_edge(i * p, i * p + 3)
            G.add_edge(i * p + 1, i * p + 4)
            G.add_edge(i * p + 2, i * p + 5)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)
        G.add_edge((number_of_blocks - 1) * p + 1, (number_of_blocks - 1) * p + 2)

    elif archi == "ladder_3_diagonals":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)
            G.add_edge(i * p + 1, i * p + 2)

            G.add_edge(i * p, i * p + 3)
            G.add_edge(i * p, i * p + 4)
            G.add_edge(i * p + 1, i * p + 3)
            G.add_edge(i * p + 1, i * p + 4)
            G.add_edge(i * p + 1, i * p + 5)
            G.add_edge(i * p + 2, i * p + 4)
            G.add_edge(i * p + 2, i * p + 5)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)
        G.add_edge((number_of_blocks - 1) * p + 1, (number_of_blocks - 1) * p + 2)

    elif archi == "ladder_4":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)
            G.add_edge(i * p + 1, i * p + 2)
            G.add_edge(i * p + 2, i * p + 3)

            G.add_edge(i * p, i * p + 4)
            G.add_edge(i * p + 1, i * p + 5)
            G.add_edge(i * p + 2, i * p + 6)
            G.add_edge(i * p + 3, i * p + 7)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)
        G.add_edge((number_of_blocks - 1) * p + 1, (number_of_blocks - 1) * p + 2)
        G.add_edge((number_of_blocks - 1) * p + 2, (number_of_blocks - 1) * p + 3)

    elif archi == "ladder_4_diagonals":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)
            G.add_edge(i * p + 1, i * p + 2)
            G.add_edge(i * p + 2, i * p + 3)

            G.add_edge(i * p, i * p + 4)
            G.add_edge(i * p, i * p + 5)
            G.add_edge(i * p + 1, i * p + 4)
            G.add_edge(i * p + 1, i * p + 5)
            G.add_edge(i * p + 1, i * p + 6)
            G.add_edge(i * p + 2, i * p + 5)
            G.add_edge(i * p + 2, i * p + 6)
            G.add_edge(i * p + 2, i * p + 7)
            G.add_edge(i * p + 3, i * p + 6)
            G.add_edge(i * p + 3, i * p + 7)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)
        G.add_edge((number_of_blocks - 1) * p + 1, (number_of_blocks - 1) * p + 2)
        G.add_edge((number_of_blocks - 1) * p + 2, (number_of_blocks - 1) * p + 3)

    elif archi == "grid":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)
            G.add_edge(i * p, i * p + 2)
            G.add_edge(i * p + 1, i * p + 3)
            G.add_edge(i * p + 2, i * p + 3)

            G.add_edge(i * p + 2, i * p + 4)
            G.add_edge(i * p + 3, i * p + 5)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)
        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 2)
        G.add_edge((number_of_blocks - 1) * p + 1, (number_of_blocks - 1) * p + 3)
        G.add_edge((number_of_blocks - 1) * p + 2, (number_of_blocks - 1) * p + 3)

    elif archi == "grid_diagonals":
        for i in range(number_of_blocks - 1):
            G.add_edge(i * p, i * p + 1)
            G.add_edge(i * p, i * p + 2)
            G.add_edge(i * p, i * p + 3)
            G.add_edge(i * p + 1, i * p + 2)
            G.add_edge(i * p + 1, i * p + 3)
            G.add_edge(i * p + 2, i * p + 3)

            G.add_edge(i * p + 2, i * p + 4)
            G.add_edge(i * p + 2, i * p + 5)
            G.add_edge(i * p + 3, i * p + 4)
            G.add_edge(i * p + 3, i * p + 5)

        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 1)
        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 2)
        G.add_edge((number_of_blocks - 1) * p, (number_of_blocks - 1) * p + 3)
        G.add_edge((number_of_blocks - 1) * p + 1, (number_of_blocks - 1) * p + 2)
        G.add_edge((number_of_blocks - 1) * p + 1, (number_of_blocks - 1) * p + 3)
        G.add_edge((number_of_blocks - 1) * p + 2, (number_of_blocks - 1) * p + 3)

    elif archi == "all_to_all" or archi == "all_to_all_3" or archi == "all_to_all_4":
        for i in range(number_of_blocks - 1):
            for j in range(p):
                for k in range(p):
                    if j != k:
                        G.add_edge(i * p + j, i * p + k)
                    G.add_edge(i * p + j, (i + 1) * p + k)

        for j in range(p):
            for k in range(p):
                if j != k:
                    G.add_edge(
                        (number_of_blocks - 1) * p + j, (number_of_blocks - 1) * p + k
                    )

    return G


def upper_bound_depth(n, archi, p):
    if archi == "all_to_all":
        return int((2 + 3 / p) * n + 2 * p + 6)
    elif archi == "lnn":
        return int(5 * n)
    elif archi == "line_2":
        return int(5 * n)
    elif archi == "ladder_2":
        return int(4 * n + 3)
    elif archi == "ladder_2_diagonals":
        return int(7 * n / 2 + 3)
    elif archi == "ladder_3":
        return int(11 * n / 3 + 8)
    elif archi == "ladder_3_diagonals":
        return int(3 * n + 8)
    elif archi == "ladder_4":
        return int(13 * n / 4 + 10)
    elif archi == "ladder_4_diagonals":
        return int(5 * n / 2 + 10)
    elif archi == "grid":
        return int(4 * n + 8)
    elif archi == "grid_diagonals":
        return int(15 * n / 4 + 6)
    elif archi == "all_to_all_3":
        return int(7 * n / 3 + 6)
    elif archi == "all_to_all_4":
        return int(7 * n / 4 + 6)


# Linear algebra functions


def random_operator(n, k):
    A = np.eye(n, dtype=bool)
    for i in range(k):
        control = np.random.randint(n)
        target = np.random.randint(n)
        while target == control:
            target = np.random.randint(n)

        A[target, :] ^= A[control, :]

    return A


def row_echelon(M):
    m = np.size(M, 0)
    n = np.size(M, 1)
    A = np.copy(M)
    h = 0
    k = 0
    while h < m and k < n:
        pivot = -1
        for i in range(h, m):
            if A[i, k]:
                pivot = i
                break
        if pivot == -1:
            k += 1
        else:
            A[[h, pivot], :] = A[[pivot, h], :]
            for i in range(h + 1, m):
                if A[i, k]:
                    for j in range(n):
                        A[i, j] = A[i, j] ^ A[h, j]
            h += 1
            k += 1
    return A


def rank(M):
    A_reduced = row_echelon(M)
    rank = 0
    for i in range(np.size(M, 0)):
        if (A_reduced[i, :] != np.zeros(np.size(M, 1))).any():
            rank += 1
    return rank


def inverse(M):
    n = np.size(M, 0)
    assert np.size(M, 1) == n
    B = np.copy(M)
    A = np.eye(n, dtype=bool)
    for i in range(n):
        if not B[i, i]:
            j = i + 1
            while not B[j, i]:
                j += 1
            for k in range(n):
                B[i, k] = B[i, k] ^ B[j, k]
                A[i, k] = A[i, k] ^ A[j, k]
        for j in range(n):
            if B[j, i] and i != j:
                for k in range(n):
                    B[j, k] = B[i, k] ^ B[j, k]
                    A[j, k] = A[i, k] ^ A[j, k]
    assert (B == np.eye(n)).all()
    return A


def GF2_prod(A, B):
    n = np.size(A, 0)
    m = np.size(B, 1)
    C = np.zeros((n, m), dtype=bool)
    for i in range(n):
        for j in range(m):
            for k in range(np.size(A, 1)):
                C[i, j] ^= A[i, k] * B[k, j]
    return C


def UPL(A):
    """The UPL decomposition.
    Input:
        A: the boolean np.array of the linear reversible
           operator to synthesize.
    Returns:
        U: the upper triangular matrix U with
        columns  permuted that has to be reduces to
        an upper triangular operator

        labels: the labels to sort
    """

    B = np.copy(A)
    n = np.size(A, 0)
    U = np.eye(n, dtype=bool)

    for i in range(n - 1, -1, -1):
        pivot = n - 1
        if np.count_nonzero(A[:, i]):
            while not A[pivot, i]:
                pivot -= 1

            for j in range(n):
                if j != pivot and A[j, i]:
                    A[j, :] ^= A[pivot, :]
                    U[:, pivot] ^= U[:, j]
            for j in range(i - 1, -1, -1):
                if A[pivot, j]:
                    A[:, j] ^= A[:, i]

    labels = -np.ones(n, dtype=int)
    for i in range(n):
        j = 0
        while not A[i, j]:
            j += 1

        labels[i] = n - j - 1

    inv_labels = np.zeros(len(labels), dtype=int)
    for i in range(len(labels)):
        inv_labels[labels[i]] = i

    U = U[:, inv_labels]

    return U, labels


# The synthesis functions


def zero_block(V, step):
    p = np.size(V, 1)
    circ = []

    if step == 2:
        p = p // 2
        V = V[:, :p]

        for i in range(p):
            circ.append([i, i + p])
            circ.append([i + p, i])

        B = V[:p, :] ^ V[p:, :]
        V[:p, :] = V[p:, :]
        V[p:, :] = B
        circ += standard_zero_block(V[:p, :], V[p:, :])
    else:
        j = 0
        rows = [i for i in range(p)]
        r = rank(V[:p, :])
        while r < p:
            for row in rows:
                V[j, :] ^= V[row + p, :]
                new_r = rank(V[:p, :])
                if new_r > r:
                    r = new_r
                    circ.append([row + p, j])
                    rows = list(set(rows) - set([row]))
                    break
                else:
                    V[j, :] ^= V[row + p, :]
            j += 1

        circ += standard_zero_block(V[:p, :], V[p:, :])

    return circ


def standard_zero_block(A, B):
    p = np.size(A, 0)
    circ = []

    invA = inverse(A)
    C = GF2_prod(B, invA)

    while np.count_nonzero(C) > 0:
        G = nx.Graph()
        G.add_nodes_from([i for i in range(p)], bipartite=0)
        G.add_nodes_from([i for i in range(p, 2 * p)], bipartite=1)
        for i in range(p):
            for j in range(p):
                if C[j, i]:
                    G.add_edges_from([(i, j + p)])

        matching = nx.bipartite.maximum_matching(G, top_nodes=[i for i in range(p)])

        treated_vertices = {}
        for v1, v2 in matching.items():
            if v1 not in treated_vertices:
                if v1 < v2:
                    circ.append([v1, v2])
                    treated_vertices[v1] = 1
                    treated_vertices[v2] = 1
                    C[v2 - p, v1] ^= True
                else:
                    circ.append([v2, v1])
                    treated_vertices[v1] = 1
                    treated_vertices[v2] = 1
                    C[v1 - p, v2] ^= True

    return circ


def sort_labels_block_local(U, p, labels, rows, step, database, circ):
    """Sort two blocks of labels."""

    sorted_labels = np.sort(labels[rows])

    if step == 1:
        columns_to_consider = sorted_labels[:p]
    elif step == 2:
        columns_to_consider = sorted_labels

    V = U[rows, :][:, columns_to_consider]

    if database == []:
        # all_to_all case
        local_circ = zero_block(V, step)
    else:
        local_circ = database.get_circuit(V, step=step)

    for gate in local_circ:
        circ.append(([rows[gate[0]], rows[gate[1]]]))
        U[rows[gate[1]], :] ^= U[rows[gate[0]], :]

    labels[rows] = sorted_labels


def sort_labels_block(U, labels, p, step, database):
    """An intermediate function to sort the labels by block."""

    n = np.size(U, 0)
    assert (
        n % p == 0
    ), "The size of the operator has to be a multiple of the block size."
    number_of_blocks = n // p
    circ = []

    shift = 0
    for i in range(n // p):
        for j in range((number_of_blocks - shift) // 2):
            block1 = [
                i for i in range(2 * p * j + p * shift, 2 * p * j + p * shift + p)
            ]
            block2 = [
                i for i in range(2 * p * j + p * shift + p, 2 * p * (j + 1) + p * shift)
            ]
            sort_labels_block_local(U, p, labels, block1 + block2, step, database, circ)

        shift = (shift + 1) % 2

    return circ


def direct_synthesis(A, archi, database):
    if archi == "all_to_all" or archi == "lnn":
        circ = []
        for i in range(np.size(A, 0)):
            if not A[i, i]:
                for j in range(i + 1, np.size(A, 0)):
                    if A[j, i]:
                        circ.append([j, i])
                        A[i, :] ^= A[j, :]
                        break
            for j in range(np.size(A, 0)):
                if A[j, i] and i != j:
                    circ.append([i, j])
                    A[j, :] ^= A[i, :]
    else:
        circ = database.get_circuit(A, step=3)

    return circ


def synthesis(A, p, archi):
    """The function to perform the synthesis.

    Inputs:
        A: the boolean np.array of the linear reversible
           operator to synthesize.
        p: the block size. The size of A has to be a multiple of p!
        archi: the name of the architecture of the hardware.
               archi is a string!
               You can choose between:
                    - lnn
                    - line_2
                    - ladder_2
                    - ladder_2_diagonals
                    - ladder_3
                    - ladder_3_diagonals
                    - ladder_4
                    - ladder_4_diagonals
                    - grid
                    - grid_diagonals
                    - all_to_all_3
                    - all_to_all_4
                    - all_to_all
    Output:
        circ: the circuit hardware compliant implementing A.
              circ is given as a list of tuples of integers (control, target)
              for each CNOT.
    """

    n = np.size(A, 0)
    assert (
        n % p == 0
    ), "The size of the operator has to be a multiple of the block size."

    if archi == "all_to_all" or archi == "lnn":
        database = []
    else:
        database = Database(archi)

    U, labels = UPL(np.copy(A))
    circ = sort_labels_block(U, labels, p, 1, database)
    d1 = depth(circ)
    # print("Depth step 1 : ", d1)
    NW = apply_circuit(np.copy(A), circ)

    labels = np.array([np.size(A, 0) - i - 1 for i in range(np.size(A, 0))])

    circ2 = sort_labels_block(NW, labels, p, 2, database)
    d2 = depth(circ2)
    # print("Depth step 2 : ", d2)
    circ += circ2

    for i in range(n // p):
        circ += shift(
            direct_synthesis(
                NW[p * i : p * (i + 1), p * i : p * (i + 1)], archi, database
            ),
            p * i,
        )

    # print("Depth step 3 : ", depth(circ)-d1-d2)
    return circ[::-1]


# Tests
def test_synthesis(n, p, archi, niter):
    print("Upper bound   : ", upper_bound_depth(n, archi, p))
    for i in range(niter):
        A = random_operator(n, n * n)
        circ = synthesis(np.copy(A), p, archi)
        check_circuit(A, circ, archi, p)
        print("Depth circuit : ", depth(circ))
        assert depth(circ) <= upper_bound_depth(
            n, archi, p
        ), "The depth does not match the upper bound."
