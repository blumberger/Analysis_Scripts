import numpy as np
import os

from src.io_utils import general_io as gen_io
from src.input_file import input_file_types as inp_types


class Numpy(object):
	"""
	Will write numpy files.
	"""
	def __init__(self, Data_Class, filepath):

		if type(Data_Class) == inp_types.Vars:
			for data_key in Data_Class:
				Data = Data_Class[data_key]
				if hasattr(Data, "get_numpy_arrays"):
					numpy_arrs = Data.get_numpy_arrays()
				else:
					raise SystemError(f"Can't write numpy arrays for '{Data}'.")

				self._write_files_(numpy_arrs, filepath)

		else:
			raise SystemError(f"Can't currently write numpy files from {Data_Class}")

	def _write_files_(self, arrs, filepath):
		"""
		Will loop over numpy arrays and write them
		"""
		# If there's a dot in the filepath prepend the filepath
		if '.' in filepath:
			filename, ext = gen_io.remove_file_extension(filepath)
			arrs = {f"{filename}_{key}.{ext}": arrs[key] for key in arrs}
			print(arrs.keys())
		
		# Else just assume it's a folder
		else:
			if not os.path.isdir(filepath):
				os.makedirs(filepath)

			arrs = {f"{filepath}/{key}.npy": arrs[key] for key in arrs}

		# Now write the files.
		for key in arrs:
			filepath = gen_io.create_unique_filepath(key)
			np.save(key, arrs[key])