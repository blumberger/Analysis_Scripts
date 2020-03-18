#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to parse the parameter files.
"""
import re

from src.io_utils import general_io as gen_io

from src.system import type_checking as type_check

from src.parsing import general_parsing as gen_parse


class Params(gen_io.DataFileStorage):
    """
    Will parse the parameter files.

    A parameter file takes the form:
        SECTION:
        param1 = var1 var2 ... varN
        .
        .
        .

    One can specify any number (greater than 1) of sections and params and vars.
    """
    def parse(self):
        """
        Will parse the file
        """
        self.data = {}
        self.file_txt = '\n'.join((i for i in self.file_txt.split("\n") if i))

        self.get_sects()

        for sect in self.sect_txt:
            self.data[sect] = {}
            for line in self.sect_txt[sect]:
                splitter = line.split("=")
                vals = '='.join(splitter[1:]).strip().split()
                if len(vals) == 1:
                    self.data[sect][splitter[0].strip()] = type_check.eval_type(vals[0])
                else:
                    self.data[sect][splitter[0].strip()] = [type_check.eval_type(i)
                                                            for i in vals]

    def get_sects(self):
        """
        Will get the sections and store the raw txt in self.data

        Outputs:
            <dict> self.sect_txt => Keys are sections, vals are txt split by \n
        """
        self.sect_txt = {}
        all_sects = re.findall("[a-zA-Z_ ]+:", self.file_txt)

        # Define a couple of functions
        tidy_name = lambda s: s.rstrip(':').lower()
        tidy_sect = lambda s: s.strip("\n").split("\n")

        # Inital Section
        tmp = all_sects[0].join(self.file_txt.split(all_sects[0])[1:])
        prev_name = tidy_name(all_sects[0])

        # Middle Sections
        for sect_i in range(len(all_sects)-1):
            sect = all_sects[sect_i + 1]

            sect_name = tidy_name(sect)
            splitter = tmp.split(sect)
            prev_sect = splitter[0]

            self.sect_txt[prev_name] = tidy_sect(prev_sect)
            prev_name = sect_name
            tmp = sect.join(splitter[1:])

        # Last Section
        sect_name = tidy_name(all_sects[-1])
        self.sect_txt[prev_name] = tidy_sect(tmp)
