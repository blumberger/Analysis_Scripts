"""
A module to store derived types (classes) used in the parsing of input files.
"""

class Variable(object):
    """
    A class to hold the variables parsed from an input file.

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
        return f"Variable '{self.name}':" + "\n" + str(self.data)
    def __float__(self):
        """return float(data)"""
        return float(self.data)
    def __int__(self):
        """return int(data)"""
        return int(self.data)

    # Overload the indexing Functions
    def __getitem__(self, key):
        """         Indexing affects metadata        """
        if key in self.metadata:
            return self.metadata[key]
        return False

    def __setitem__(self, key, value):
        """         Indexing affects metadata        """
        self.metadata[key] = value
        try:  # set data's metadata too
            self.data.metadata[key] = value
        except AttributeError:
            pass

    def set_data_var(self):
        """
        A function to find what sort of data variable the data holds
        """
        for var_name in dir(self.data):
            attr = getattr(self.data, var_name)
            if 'data' in var_name and not callable(attr):
                if 'meta' in var_name or '_' == var_name[0]:
                    continue
                elif isinstance(attr, (int, float)):
                    continue
                elif type(attr) == str and len(attr) < 10:
                    continue
                else:
                    self.data.data = attr

    # Overload appending
    def append(self, val):
        """
        Check the data is a list, if not don't append
        """
        if isinstance(self.data, (int, str, float, tuple)):
            raise TypeError(f"Cannot append to variable '{self.name}' which is of type {type(self.data)}.")

        elif type(self.data) == dict:
            self.data.setdefault('appended_vals', []).append(self.data)

        else:
            self.data.append(val)
            return self


    def __len__(self):
        """
        Return a length of the variable
        """
        if isinstance(self.data, (list, dict, tuple)):
            return len(self.data)
        else:
            return 1

    # Overload mathematical operators
    def __add__(self, val):
        """Ammend data attribute and return self"""
        self.data += val
        return self
    def __radd__(self, val):
        """Ammend data attribute and return self"""
        self.data += val
        return self
    def __sub__(self, val):
        """Ammend data attribute and return self"""
        self.data -= val
        return self
    def __rsub__(self, val):
        """Ammend data attribute and return self"""
        self.data -= val
        return self
    def __mul__(self, val):
        """Ammend data attribute and return self"""
        self.data *= val
        return self
    def __rmul__(self, val):
        """Ammend data attribute and return self"""
        self.data *= val
        return self
    def __truediv__(self, val):
        """Ammend data attribute and return self"""
        self.data /= val
        return self
    def __rtruediv__(self, val):
        """Ammend data attribute and return self"""
        self.data = val / self.data
        return self
    def __floordiv__(self, val):
        """Ammend data attribute and return self"""
        self.data //= val
        return self
    def __rfloordiv__(self, val):
        """Ammend data attribute and return self"""
        self.data = val // self.data
        return self
    def __pow__(self, val):
        """Ammend data attribute and return self"""
        self.data **= val
        return self
    def __rpow__(self, val):
        """Ammend data attribute and return self"""
        self.data **= val
        return self
