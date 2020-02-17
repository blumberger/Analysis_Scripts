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
        if type(Data_Class.df_data) == dict:
            for key in Data_Class.df_data:
                new_filepath = f"{filepath}_{key}.csv"
                if type(Data_Class.df_data[key]) == pd.DataFrame:
                    Data_Class.df_data[key].to_csv(new_filepath, index=False)
                else:
                    print(f"Can't save data {key}")
