"""
A module to store derived types (classes) used in the parsing of input files.
"""
import numpy as np
from collections import Counter

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
    name = "Input Variable"
    def __init__(self, var_name, var_data, metadata={}):
       self.name = var_name

       self.is_dict = False
       if type(var_data) == dict:
            self.is_dict = False
            self.data = Vars(var_data)
       else:
            self.data = var_data
       self.metadata = metadata


    # Overload the type convertors
    def __str__(self):
        """return str(data)"""
        return f"Variable '{self.name}':" + "\n" + str(self.data)
    def __repr__(self):
        return str(self.data)
    def __float__(self):
        """return float(data)"""
        return float(self.data)
    def __int__(self):
        """return int(data)"""
        return int(self.data)

    def __contains__(self, val):
        """Check if the value is in the metadata."""
        return val in self.metadata

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

    def metadata_update(self, dict_):
        self.metadata.update(dict_)
        self.data.metadata.update(dict_)

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

    def get_xyz_data(self):
        return self.data.get_xyz_data()

    def get_xyz_cols(self):
        return self.data.get_xyz_cols()

    def get_xyz_timesteps(self):
        return self.data.get_xyz_timesteps()

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
        """
        Will add the data stored in this variable by a value.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = self.data[key] + val[val_keys[0]]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : self.data + val[val_keys[0]]})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = self.data[key] + val[key]
                else:
                    raise SystemError("No scheme for adding data.\n\nTrying to multiply"
                                + f"{self.data} + {val}")

            else:
                raise SystemError("No scheme for adding data.\n\nTrying to multiply"
                                + f"{self.data} + {val}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = self.data[key] + val
            else:
                self.data = self.data + val

        return self

    def __radd__(self, val):
        """
        Will add the data stored in this variable by a value on the right.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = val[val_keys[0]] + self.data[key]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : val[val_keys[0]] + self.data})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = val[key] + self.data[key]
                else:
                    raise SystemError("No scheme for adding data.\n\nTrying to multiply"
                                + f"{val} + {self.data}")

            else:
                raise SystemError("No scheme for adding data.\n\nTrying to multiply"
                                + f"{val} + {self.data}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = val + self.data[key]
            else:
                self.data = val + self.data

        return self

    def __sub__(self, val):
        """
        Will subtract the data stored in this variable by a value.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = self.data[key] - val[val_keys[0]]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : self.data - val[val_keys[0]]})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = self.data[key] - val[key]
                else:
                    raise SystemError("No scheme for subtracting data.\n\nTrying to multiply"
                                + f"{self.data} - {val}")

            else:
                raise SystemError("No scheme for subtracting data.\n\nTrying to multiply"
                                + f"{self.data} - {val}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = self.data[key] - val
            else:
                self.data = self.data - val

        return self

    def __rsub__(self, val):
        """
        Will subtract the data stored in this variable by a value on the right.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = val[val_keys[0]] - self.data[key]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : val[val_keys[0]] - self.data})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = val[key] - self.data[key]
                else:
                    raise SystemError("No scheme for subtracting data.\n\nTrying to multiply"
                                + f"{val} - {self.data}")

            else:
                raise SystemError("No scheme for subtracting data.\n\nTrying to multiply"
                                + f"{val} - {self.data}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = val - self.data[key]
            else:
                self.data = val - self.data

        return self

    def __mul__(self, val):
        """
        Will multiply the data stored in this variable by a value.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = self.data[key] * val[val_keys[0]]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : self.data * val[val_keys[0]]})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = self.data[key] * val[key]
                else:
                    raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{self.data} * {val}")

            else:
                raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{self.data} * {val}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = self.data[key] * val
            else:
                self.data = self.data * val

        return self

    def __rmul__(self, val):
        """
        Will multiply the data stored in this variable by a value on the right.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = val[val_keys[0]] * self.data[key]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : val[val_keys[0]] * self.data})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = val[key] * self.data[key]
                else:
                    raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{val} * {self.data}")

            else:
                raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{val} * {self.data}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = val * self.data[key]
            else:
                self.data = val * self.data

        return self

    def __truediv__(self, val):
        """
        Will divide the data stored in this variable by a value.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = self.data[key] / val[val_keys[0]]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : self.data / val[val_keys[0]]})

            # If there are multiple keys check we can divide the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = self.data[key] / val[key]
                else:
                    raise SystemError("No scheme for dividing data.\n\nTrying to divide"
                                + f"{self.data} / {val}")

            else:
                raise SystemError("No scheme for dividing data.\n\nTrying to divide"
                                + f"{self.data} / {val}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = self.data[key] / val
            else:
                self.data = self.data / val

        return self

    def __rtruediv__(self, val):
        """
        Will divide the data stored in this variable by a value on the right.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = val[val_keys[0]] / self.data[key]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : val[val_keys[0]] / self.data})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = val[key] / self.data[key]
                else:
                    raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{val} / {self.data}")

            else:
                raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{val} / {self.data}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = val / self.data[key]
            else:
                self.data = val / self.data

        return self

    def __floordiv__(self, val):
        """
        Will divide the data stored in this variable by a value and truncate to an integer.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = self.data[key] // val[val_keys[0]]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : self.data // val[val_keys[0]]})


            # If there are multiple keys check we can divide the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = self.data[key] // val[key]
                else:
                    raise SystemError("No scheme for dividing data.\n\nTrying to divide"
                                + f"{self.data} // {val}")

            else:
                raise SystemError("No scheme for dividing data.\n\nTrying to divide"
                                + f"{self.data} // {val}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = self.data[key] // val
            else:
                self.data = self.data // val

        return self

    def __rfloordiv__(self, val):
        """
        Will divide the data stored in this variable by a value on the right and truncate to an integer.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = val[val_keys[0]] // self.data[key]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : val[val_keys[0]] // self.data})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = val[key] // self.data[key]
                else:
                    raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{val} // {self.data}")

            else:
                raise SystemError("No scheme for multiplying data.\n\nTrying to multiply"
                                + f"{val} // {self.data}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = val // self.data[key]
            else:
                self.data = val // self.data

        return self

    def __pow__(self, val):
        """
        Will exponentiate the data stored in this variable by a value.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = self.data[key] ** val[val_keys[0]]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : self.data ** val[val_keys[0]]})

            # If there are multiple keys check we can exponentiate the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = self.data[key] ** val[key]
                else:
                    raise SystemError("No scheme for exponentiating data.\n\nTrying to exponentiate"
                                + f"{self.data} ** {val}")

            else:
                raise SystemError("No scheme for exponentiating data.\n\nTrying to exponentiate"
                                + f"{self.data} ** {val}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = self.data[key] ** val
            else:
                self.data = self.data ** val

        return self

    def __rpow__(self, val):
        """
        Will exponentiate the data stored in this variable by a value on the right.

        Because this data structure has become a bit complicated the mathematical
        operations need to be handled in a special way.

        The general scheme is outlined below:

            * If the val to operate on this data is a dict (Vars) obj then check the keys.
                * If there is only one key, use that as the data.
                * If there are many keys and if the keys are the same as in this structure
                  then use the operator on like keys.

            * If the val isn't a dict then try operating on all data in this structure.
        """
        # Check if the variable is a dictionary or not.
        self.is_dict = type(self.data) == Vars

        # If the val is a dict
        if type(val) == Vars:
            val_keys = list(val.keys())

            # If there is only one key then treat the val like 1 variable
            if len(val_keys) == 1:
                if self.is_dict:
                    for key in self.data:
                        self.data[key] = val[val_keys[0]] ** self.data[key]
                else:
                    self.is_dict = True
                    self.data = Vars({val_keys[0] : val[val_keys[0]] ** self.data})

            # If there are multiple keys check we can multiply the data
            elif self.is_dict:
                self_keys = list(self.data.keys())
                if self_keys == val_keys:
                    for key in self.data:
                        self.data[key] = val[key] ** self.data[key]
                else:
                    raise SystemError("No scheme for exponentiating data.\n\nTrying to raise"
                                + f"{val} to the power of {self.data}")

            else:
                raise SystemError("No scheme for exponentiating data.\n\nTrying to raise"
                                + f"{val} to the power of {self.data}")

        # If the val isn't a dict
        else:
            if self.is_dict:
                for key in self.data:
                    self.data[key] = val ** self.data[key]
            else:
                self.data = val ** self.data

        return self

    def __eq__(self, val):
        return self.data == val

    def __ge__(self, val):
        return self.data > val



class Vars(dict):
    """
    A subclass of the dictionary object.

    This acts exactly the same as a dictionary but has a few more bells and
    whistles to make the data easier to call forward.
    """
    name = "Vars Dict"
    metadata = {}

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __repr__(self):
        dictrepr = dict.__repr__(self)
        return '%s(%s)' % (type(self).__name__, dictrepr)

    def update(self, *args, **kwargs):
        for k, v in dict(*args, **kwargs).items():
            self[k] = v


    ###############################################################################
    ##                XYZ Data Handling Functions                                ##
    ###############################################################################

    ###########################################
    ### XYZ columns
    def __get_ats_per_mol_from_lammps_dump__(self, df):
        """
        Will get the number of atoms per molecule from the lammps dump file.

        Inputs:
            * df <DataFrame> => The lammps dump data.

        N.B. Only works for systems with 1 mol type.
        """
        if 'mol' not in df.columns:
            raise SystemExit("I can't find which molecule each atom belongs to...\n"
                             +"The Lammps dump file needs a mol column.\n\n")

        unique_vals = np.unique(df['mol'])
        unique_steps = df['timestep'].unique()
        return len(df[(df['mol'] == unique_vals[0]) & (df['timestep'] == unique_steps[0])])

    def __get_xyz_cols_from_lammps_dump__(self):
        """
        Will get the element type columns from a csv file.
        """
        if self.number_each_atom is False:
            raise SystemExit("I need to know how many of each type of atom there are in each mol."
                           + "\n\nSet this by using the command:\n\n\t`set system <data_name>"
                           + " to <mol_type>`\n\nin the input file.")

        df = self.__get_lammps_csv_data__()

        ats_per_mol = self.__get_ats_per_mol_from_lammps_dump__(df)

        # Error check
        if 'type' not in df.columns:
            raise SystemError("\n\n\nI can't calculate the xyz cols of the mols."
                              + " I can't find the atoms types in the data."
                              + "I don't what atom types are in the mol")

        # Error check -if there is a mol with equal nums of atom of x and y type.
        num_elm_in_mol = [self.number_each_atom[i] for i in self.number_each_atom]
        if len(set(num_elm_in_mol)) != len(num_elm_in_mol):
            raise SystemError("\n\n\nCan't find the xyz columns."
                            + " I don't know what types the atoms are and can't"
                            + " work them out as the molecule has 2 or more "
                            + "elements with the same number of atoms.\n\n")

        # Compare how many atoms of each type are in each molecule and how many there should be.
        num_at_types = Counter(df.loc[:ats_per_mol-1, 'type'])
        at_types = {i: self.number_each_atom[i] for i in self.number_each_atom}
        cvt_type = {}
        for i in num_at_types:
            for elm in at_types:
                if at_types[elm] == num_at_types[i]:
                    # print(f"Atom type '{i}' is element {elm}.")
                    cvt_type[str(i)] = elm
                    break
            else:
                if sum(at_types.values()) != sum(num_at_types.values()):
                    msg = "Could not determine what element each atom type is."
                    msg = "\n\nExpected num atoms = " + f"{sum(at_types.values())}"
                    msg = "\n\nTotal num atoms = " + f"{sum(num_at_types.values())}"
                    raise SystemError("Could not determine what element each atom type is.")
                else:
                    for i in num_at_types:
                        for elm in at_types:
                            if at_types[elm] == num_at_types[i]:
                                cvt_type[str(i)] = elm 
                                break
                    at_types.pop(elm)
                    num_at_types.pop(i)

                    if len(at_types) == 1:
                        print(num_at_types, at_types)
                    print(cvt_type)
                    raise SystemExit
            at_types.pop(elm)

        # Create the cols array
        cols = [df['type'][df['timestep'] == i].to_numpy().astype(str)
                for i in df['timestep'].unique()]
        cols = np.array(cols)

        self.natom = len(cols[0])
        for i in np.unique(cols):
            cols[cols == i] = cvt_type[i]

        return np.array(cols)

    def __get_xyz_cols_from_xyz__(self):
        """Will return the xyz columns."""
        return self['xyz'].cols

    def __get_xyz_cols_from_lammps_inp__(self):
        at_types = self['lammps_input'].csv_data['atoms']['at_type'].unique()
        try:
            new_at_types = {at_type: self['lammps_input'].metadata['atom_types'][at_type]['details']['abbreviation'] 
                         for at_type in at_types}
        except KeyError:
            new_at_types =at_types
        at_types = new_at_types


        cols = self['lammps_input'].csv_data['atoms']['at_type'].to_numpy().astype(str)
        for i in at_types:
            cols[cols == str(i)] = at_types[i]

        cols = np.reshape(cols, (len(cols), 1))
        return cols

    def __get_xyz_cols_from_lammps_data__(self):
        at_types = self['lammps_data'].csv_data['atoms']['at_type'].unique()
        try:
            new_at_types = {at_type: self['lammps_data'].metadata['atom_types'][at_type]['details']['abbreviation'] 
                         for at_type in at_types}
        except KeyError:
             new_at_types =at_types

        at_types = new_at_types


        cols = self['lammps_data'].csv_data['atoms']['at_type'].to_numpy().astype(str)
        for i in at_types:
            cols[cols == str(i)] = at_types[i]

        cols = np.reshape(cols, (len(cols), 1))
        return cols


    def get_xyz_cols(self, number_each_atom=False):
        """
        Will try to get the atom element types from the number of each type in a
        molecule.
        """
        self.number_each_atom = number_each_atom
        if number_each_atom is False and 'number_each_atom' in self.metadata:
            self.number_each_atom = self.metadata['number_each_atom']

        XYZ_DATA_KEYS = {'xyz': self.__get_xyz_cols_from_xyz__,
                         'lammps_dump': self.__get_xyz_cols_from_lammps_dump__,
                         'lammps_input': self.__get_xyz_cols_from_lammps_inp__,
                         'lammps_data': self.__get_xyz_cols_from_lammps_data__,
                        }
        xyz_data = [XYZ_DATA_KEYS[i]() for i in self if i in XYZ_DATA_KEYS]

        return np.array(xyz_data)

    ###########################################
    ### Timesteps
    def get_xyz_timesteps(self):
        """
        Will grab the timestep of each step.
        """
        XYZ_DATA_KEYS = {'xyz': self.__get_xyz_timesteps_from_xyz__,
                         'lammps_dump': self.__get_xyz_timesteps_from_lammps_dump__,
                         'lammps_input': lambda: [0],
                         'lammps_data': lambda: [0],
                        }
        xyz_data = [XYZ_DATA_KEYS[i]() for i in self if i in XYZ_DATA_KEYS]
        return np.array(xyz_data)

    def __get_xyz_timesteps_from_lammps_inp__(self):
        #if 'time' in self['lammps_input'].csv_data['atoms']
        return [0]
    

    def __get_xyz_timesteps_from_xyz__(self):
        """Will get the xyz timesteps from xyz file container."""
        return self['xyz'].timesteps

    def __get_xyz_timesteps_from_lammps_dump__(self):
        """Will return the xyz timesteps from lammps dump file container."""
        return self['lammps_dump'].csv_data['timestep'].unique()


    ###########################################
    ### XYZ Data
    def get_xyz_data(self):
        """
        Will search the data dictionary and if any keys have xyz data the will return those.

        The data types with xyz data in are given in XYZ_DATA_KEYS, this contains references
        to the functions that grab the xyz data for that data type.
        """
        XYZ_DATA_KEYS = {'xyz': self.__get_xyz_data_from_xyz__,
                         'lammps_dump': self.__get_xyz_data_from_lammps_dump__,
                         'lammps_input': self.__get_xyz_data_from_lammps_inp__,
                         'lammps_data': self.__get_xyz_data_from_lammps_data__,
                        }
        xyz_data = [XYZ_DATA_KEYS[i]() for i in dict.keys(self) if i in XYZ_DATA_KEYS]

        return np.array(xyz_data)

    def __get_xyz_data_from_xyz__(self):
        """Will return the xyz data an XYZ file type."""
        return self['xyz'].xyz_data

    def __get_xyz_data_from_lammps_inp__(self):
        xyz = self['lammps_input'].csv_data['atoms'][['x', 'y', 'z']]
        xyz = xyz.to_numpy()
        return xyz

    def __get_xyz_data_from_lammps_data__(self):
        xyz = self['lammps_data'].csv_data['atoms'][['x', 'y', 'z']]
        xyz = xyz.to_numpy()
        return xyz

    def __get_xyz_data_from_lammps_dump__(self):
        """Will return the xyz data from a lammps dump file type."""
        xyz = ('x', 'y', 'z',)
        csv_data = self.__get_lammps_csv_data__()

        # Return the wrapped data
        if all(j in csv_data for j in xyz):

            unique_steps = csv_data['timestep'].unique()
            return np.array([csv_data[['x', 'y', 'z']][csv_data['timestep'] == i].to_numpy()
                             for i in unique_steps])
        else:
            raise SystemError("\n\nNo x, y, z data in the lammps dump file.\n\n")

    def __get_lammps_csv_data__(self):
        """Will return the correct csv data from a Lammps file obj."""
        wrap = ""
        if 'coordinate_wrapping' in self.metadata:
            if self.metadata['coordinate_wrapping'] == 'unwrapped':
                wrap = "unwrap"

            elif self.metadata['coordinate_wrapping'] == 'remove_split_mols':
                self['lammps_dump'].remove_split_mols()
                wrap = "rem_split"

            elif self.metadata['coordinate_wrapping'] == "wrapped":
                wrap = "wrap"

            elif self.metadata['coordinate_wrapping'] == "unwrap_split_mols":
                self['lammps_dump'].unwrap_split_mols()
                wrap = "unwrap_split_mols"

            elif self.metadata['coordinate_wrapping'] == "add_images":
                self['lammps_dump'].unwrap_add_images()
                wrap = "add_img"

            else:
                raise SystemExit("I don't understand the coordinate_wrapping instruction." +
                                 "\n\nChoose from:" + 
                                 "\n\t* unwrapped" + 
                                 "\n\t* remove_split_mols" + 
                                 "\n\t* wrapped" + 
                                 "\n\t* unwrap_split_mols" +
                                 "\n\t* add_img"
                                 )

        # Set which data to use
        if wrap == "":  csv_data = self['lammps_dump'].csv_data
        elif wrap == "unwrap": csv_data = self['lammps_dump'].unwrapped_csv
        elif wrap == "rem_split": csv_data = self['lammps_dump'].rm_split_mols
        elif wrap == "wrap": csv_data = self['lammps_dump'].wrapped_csv
        elif wrap == "unwrap_split_mols": csv_data = self['lammps_dump'].unwrapped_csv
        elif wrap == "add_img": csv_data = self['lammps_dump'].add_img_csv
        else: raise SystemExit("Can't find correct csv data from lammps log. Please point me to where it is.")

        csv_data.index = range(len(csv_data))
        return csv_data







    def splice_xyz_data(self, xmin=False, xmax=False, ymin=False, ymax=False, zmin=False, zmax=False):
        """
        Will splice the data and make changes inplace.
        """
        XYZ_DATA_KEYS = {'xyz': self.__splice_xyz_data_in_xyz__,
                         'lammps_dump': self.__splice_xyz_data_in_lammps_dump__}

        for i in XYZ_DATA_KEYS:
            if i in self:
                XYZ_DATA_KEYS[i](xmin, xmax, ymin, ymax, zmin, zmax)


    def __splice_xyz_data_in_lammps_dump__(self, xmin=False, xmax=False, ymin=False, ymax=False, zmin=False, zmax=False):
        """Splice xyz data from lammps dump file."""
        csv_data = self['lammps_dump'].csv_data

        print(csv_data)
        raise SystemExit


    def __splice_xyz_data_in_xyz__(self, xmin=False, xmax=False, ymin=False, ymax=False, zmin=False, zmax=False):
        """Splice xyz data from a xyz file."""
        xyz_data = __get_xyz_data_from_xyz__()

        mask = np.array([np.ones(np.shape(i)[:-1], dtype=bool) for i in xyz_data])

        if xmin: mask = mask & (xyz_data[:, :, 0] > xmin)
        if xmax: mask = mask & (xyz_data[:, :, 0] < xmax)
        if ymin: mask = mask & (xyz_data[:, :, 1] > ymin)
        if ymax: mask = mask & (xyz_data[:, :, 1] < ymax)
        if zmin: mask = mask & (xyz_data[:, :, 2] > zmin)
        if zmax: mask = mask & (xyz_data[:, :, 2] < zmax)

        return np.array([xyz_data[i][mask[i]] for i in range(len(xyz_data))])
