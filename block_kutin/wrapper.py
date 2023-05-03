"""
This module provides a wrapper around the thrift database generation tools
"""
import os
from .block_kutin import generate_database, load_database, get_circuit_from_matrix


class Database:
    """
    A class managing database generation, loading, and queries for
    the block Kutin algorithm.
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
        return get_circuit_from_matrix(
            matrix.T.flatten(), self.databases[step == 2], step
        )
