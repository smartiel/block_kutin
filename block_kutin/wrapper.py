"""
This module provides a wrapper around the thrift database generation tools
"""
import os
from .block_kutin import generate_database, load_database, get_circuit_from_matrix


class Database:
    """
    A class managing database generation, loading, and queries for
    the block Kutin algorithm.


    Arguments:
        topology_name (str): the name of the topology (see the paper or the src/topology.rs file for more details)
        database_folder (optional, str): the directory in which to store the database files. Default is `~/.database`
    """

    def __init__(self, topology_name, database_folder="~/.database"):
        self.topology_name = topology_name
        self.folder = os.path.expanduser(database_folder)
        self.databases = [self._load_db(1), self._load_db(2)]

    def _load_db(self, step):
        """
        Loads the database.
        If the database file does not exist, generates it and loads it
        """
        fname = os.path.join(self.folder, f"{self.topology_name}.step_{step}.db")
        if os.path.exists(fname):
            return load_database(fname)
        if not os.path.exists(self.folder):
            os.mkdir(self.folder)
        generate_database(self.topology_name, fname, step)
        return load_database(fname)

    def get_circuit(self, matrix, step=1):
        """
        Queries the database to find the shallowest implementation of a given operator.

        The operator is specified by its F_2 matrix:

        - for step 1 operators: the matrix is expected to be full rank and to have shape (2k, k)
          where k is the block size
        - for step 2 operators: the matrix is expected to be full rank, to have shape (2k, 2k)
          where k is the block size, and to have only 0s in its bottom right k x k block.

        Argument:
            matrix (np.ndarray): a numpy array representing the operator
            step (optional, int): the step (either 1 or 2, see paper for more details)

        Returns:
            list: a list of pairs (control, target) describing the CNOT implementation
        """
        return get_circuit_from_matrix(
            matrix.T.flatten(), self.databases[step == 2], step
        )
