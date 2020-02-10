"""
A module to store derived types (classes) used in the parsing of input files.
"""

class Variable(object):
   """
   A class to hold the variable objects parsed from an input file.

   Everytime a variable is declared this class should be used to store it. It
   can be used to store not just the variable data but also any metadata declared
   in the input file etc...

   Inputs:
      * variable_name <str> => The name of the variable that is to be stored.
      * varaible_data <*>   => The data which should be stored under the variable_name

   Important Attributes:
      * name <str> => The name of the variable.
      * data <*>   => The data declared.
      * metadata <dict> => Any extra data that is relevant to the variable.
   """
   def __init__(self, var_name, var_data, metadata={}):
      self.name = var_name
      self.data = var_data
      self.metadata = metadata
      

   # Overload the type convertors
   def __str__(self):
       """return str(data)"""
       return str(self.data)
   def __float__(self):
       """return float(data)"""
       return float(self.data)
   def __int__(self):
       """return int(data)"""
       return int(self.data)

   # Overload mathematical operators
   def __add__(self, val):
       """Ammend data attribute and return self"""
       self.data += val
       return self
   def __sub__(self, val):
       """Ammend data attribute and return self"""
       self.data -= val
       return self
   def __mul__(self, val):
       """Ammend data attribute and return self"""
       self.data *= val
       return self
   def __truediv__(self, val):
       """Ammend data attribute and return self"""
       self.data /= val
       return self
   def __floordiv__(self, val):
       """Ammend data attribute and return self"""
       self.data //= val
       return self
   def __pow__(self, val):
       """Ammend data attribute and return self"""
       self.data **= val
       return self
