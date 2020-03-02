#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to read the output files from the massif valgrind tool.
"""
import numpy as np

from src.io_utils import general_io as gen_io
from src.parsing import general_parsing as gen_parse
from src.system import type_checking as type_check


class Massif_File(gen_io.DataFileStorage):
    """
    A container to store massif memory profiling data in.


    Inputs/Attributes:
        * xyz_data <numpy.array> => The parsed xyz data from an xyz file.
        * cols <numpy.array> => The parsed column data from the xyz file.
        * timesteps <numpy.array> => The parsed timesteps from the xyz file.
    """
    def __init__(self, filepath):
        self.data, self.heap_tree = {}, {}
        super().__init__(filepath)

    def parse(self):
        """
        Will do the parsing of the massif file.
        """
        ltxt = self.file_txt.split("snapshot=")
        self.get_metadata(ltxt[0])
        for snapshot in ltxt[1:]:
            self.get_vals_from_snapshot(snapshot)

        for key in ('mem_heap_B', 'mem_heap_extra_B', 'mem_stacks_B',):
            self.data[key] = np.array(self.data[key]) / 1e6

        self.peak_mem = max(self.data['mem_heap_B'])

    def get_metadata(self, section):
        """
        Will pull whatever metadata there is from the section given.

        Inputs:
            * section <str> => The text from which to pull the metadata
        Outputs:
            Everything is saved as an attribute
        """
        for line in section.split('\n'):
            splitter = line.split(':')
            if len(splitter) == 2:
                setattr(self, splitter[0].strip(),
                        type_check.eval_type(splitter[1].strip()))

    def get_vals_from_snapshot(self, snapshot):
        """
        Will get the values that have been saved in each snapshot

        Inputs:
            * snapshot <str> => The snapshot text
        Outputs:
            Everything is saved as an attribute.
        """
        ltxt = snapshot.split("\n")
        detailed_heap = False
        # Parse the bits of data at the top
        for line in ltxt:
            if '=' in line:
                splitter = line.split("=")
                if len(splitter) == 2:
                    key, val = splitter
                    val = type_check.eval_type(val)
                    self.data.setdefault(key, []).append(val)

                    if key == "heap_tree" and val == "detailed":
                        detailed_heap = True

        # Save the heap_tree
        if detailed_heap:
            heap_tree = ltxt[7:]
            snapshot_num = len(self.data['time'])
            self.heap_tree[snapshot_num] = heap_tree
