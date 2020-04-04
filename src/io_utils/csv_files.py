#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to write and read CSV files.
"""
import pandas as pd

from src.io_utils import general_io as gen_io


class Write_CSV(object):
    """
    Will write a CSV file or list of CSVs files.

    This doesn't inherit from Write_File as there it would have to overwrite
    all the methods and there wouldn't be any point.

    Inputs:
        * Data_Class <class> => The class containing all the data to be written
        * filepath <str>     => The path to the file to be written.
    """
    def __init__(self, Data_Class, filepath):
        filepath, _ = gen_io.remove_file_extension(filepath)

        csv_data = Data_Class.get_csv_data()

        # The CSV data may be stored in many formats
        # Write lists of csvs
        if type(csv_data) == list:
            if len(csv_data) == 1:
                self.write_csv(csv_data[0], f"{filepath}.csv")
            else:
                for i, df in csv_data:
                    new_filepath = f"{filepath}_{i}.csv"
                    self.write_csv(df, new_filepath)

        # Write dicts of csvs
        elif type(csv_data) == dict:
            for key in csv_data:
                new_filepath = f"{filepath}_{key}.csv"
                self.write_csv(csv_data[key], new_filepath)

        # Write plain DataFrames
        else:
            self.write_csv(csv_data, f"{filepath}.csv")

    def write_csv(self, df, filepath):
        """
        Will write the csv after checking that it is a dataframe
        """
        type_ = type(df)
        if type_ == pd.DataFrame:
            df.to_csv(filepath, index=False)
        else:
            raise IOError(f"Can't csv data to '{filepath}'."
                          + " It is of unknown type {type_}.")
