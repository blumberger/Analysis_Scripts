import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN

from src.io_utils import general_io as gen_io
from src.system import type_checking as type_check
from src.calc import molecule_utils as mol_utils


class ESP_File(gen_io.DataFileStorage):
    """
    Will read a gaussian ESP file that gives charge information.

    This inherits from gen_io.DataFileStorage see this for more info.

    Inputs:
        * filepath <str> => The path to the file to be loaded.
    """
    _write_types = ('xyz',)
    name = "EPS File"
    _defaults = {}

    def _parse_(self):
        """
        Will parse the file text and store the data.
        """
        self._ltxt = self.file_txt.split("\n")
        self._double_float = lambda str_: float(str_.replace("D", "e"))

        self.__get_units__()
        self.__get_charges_and_crds__()
        self.__get_dipole_moment__()
        self.__get_quadrupole_moment__()
        self.__get_charge_density__()



    ######################################
    ##             Writing              ##
    ######################################
    def get_numpy_arrays(self):
        """
        Will return a dictionary with the named numpy arrays that need writing.

        The dict will have the name of the array as a key and the array as a value.
        """
        return {'grid_points': self.grid_points,
                'charge_density': self.charge_density,
                'atomic_crds': self.xyz_data,
                'partial_charges': self.partial_charges,
               }

    ######################################
    ##           Visualising            ##
    ######################################
    def plot_crds(self, ax=False, show=True):
        """
        Will plot the partial charges.

        Inputs:
            * ax <plt.axis> => An optional axis, if none is provided then one will be created.
                                This must have 3D projection. The created axis will have
                                orthographic projection.
            * show <bool> => Whether to show the plot or just hold it in RAM.
        """
        if ax is False:
            f = plt.figure()
            ax = f.add_subplot(111, projection="3d", proj_type="ortho")


        colors = ['k', 'b', 'g', 'r', 'y', '#aaaaaa', '#fababa', '#ff00ff', '#00ffff', '#f37f29']
        for elm in set(self.cols):
            mask = self.cols == elm
            if elm not in mol_utils.PT_abbrv:
                raise SystemError(f"Don't recognise element type '{elm}'")

            if 'plot_color' not in mol_utils.PT_abbrv[elm]:
                raise SystemError(f"I don't know what color to plot element '{elm}'. "
                                  + "Please set the value 'plot_color' in the file: "
                                  + "'src/data/periodic_table.json'.\n\nYou can use matplotlib colors or hex.")

            color = mol_utils.PT_abbrv[elm]['plot_color']
            size = 6 * (mol_utils.PT_abbrv[elm]['atomic_weight']) ** 0.2
            crds = self.xyz_data[mask]
            ax.plot(crds[:,0], crds[:,1], crds[:,2], '.',
                    ls="none", ms=size, color=color, alpha=0.5)

        # Reshape axes
        xlim = ax.get_xlim(); ylim = ax.get_ylim(); zlim = ax.get_zlim()
        xdiff = np.diff(xlim); ydiff = np.diff(ylim); zdiff = np.diff(zlim)
        max_diff = max([xdiff, ydiff, zdiff])
        x_ext = abs(max_diff - xdiff) / 2.
        y_ext = abs(max_diff - ydiff) / 2.
        z_ext = abs(max_diff - zdiff) / 2.
        ax.set_xlim([xlim[0] - x_ext, xlim[1] + x_ext])
        ax.set_ylim([ylim[0] - y_ext, ylim[1] + y_ext])
        ax.set_zlim([zlim[0] - z_ext, zlim[1] + z_ext])

        # Make pretty
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z");

        if show:
            plt.show()

        return ax

    def plot_partial_charges(self, ax=False, show=True, sym_n_norm=False, dp=3):
        """
        Will plot the partial charges and crds on a 3D matplotlib figure.

        Inputs:
            * ax <plt.axis> => An optional axis, if none is provided then one will be created.
                                This must have 3D projection. The created axis will have
                                orthographic projection.
            * show <bool> => Whether to show the plot or just hold it in RAM.
        """
        #ax = self.plot_crds(ax, show=False)

        if ax is False:
            f = plt.figure()
            ax = f.add_subplot(111, projection="3d", proj_type="ortho")

        if sym_n_norm is True:
            groups = {0: [22, 25, 4, 7], 1: [23, 24, 5, 6], 2: [32, 33, 15, 14], 3: [31, 34, 16, 13], 4: [21, 26, 8, 3], 5: [20, 9, 2, 27], 6: [30, 35, 17, 12], 7: [19, 10, 28, 1], 8: [29, 11], 9: [18, 0]}
            for i, at_group in enumerate(groups):
                at_inds = groups[at_group]
                crds = self.xyz_data[at_inds]
                self.partial_charges[at_inds] = float(f"{np.mean(self.partial_charges[at_inds]):{dp}f}")
            self.partial_charges /= sum(self.partial_charges)

        print(self.partial_charges)
        for i, (c, (x, y, z)) in enumerate(zip(self.partial_charges, self.xyz_data)):
            color = 'b' if c >= 0 else 'r'
            ax.scatter([x], [y], [z], s=(abs(c)**2)*5000, color=color, alpha=0.2)
            ax.text(x-0.3, y, z, f"{c:.{dp}f}", fontsize=17, ha="right", va="top")  # f"{c:.2f}", fontsize=17)

        # Reshape axes
        xlim = ax.get_xlim(); ylim = ax.get_ylim(); zlim = ax.get_zlim()
        xdiff = np.diff(xlim); ydiff = np.diff(ylim); zdiff = np.diff(zlim)
        max_diff = max([xdiff, ydiff, zdiff])
        x_ext = abs(max_diff - xdiff) / 2.
        y_ext = abs(max_diff - ydiff) / 2.
        z_ext = abs(max_diff - zdiff) / 2.
        ax.set_xlim([xlim[0] - x_ext, xlim[1] + x_ext])
        ax.set_ylim([ylim[0] - y_ext, ylim[1] + y_ext])
        ax.set_zlim([zlim[0] - z_ext, zlim[1] + z_ext])

        # Make pretty
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z");

        if show:
            plt.show()


    def plot_charge_density(self, ax=False, show=True):
        """
        Will plot the charge density as an isosurface.

        Inputs:
            * ax <plt.axis> => An optional axis, if none is provided then one will be created.
                                This must have 3D projection. The created axis will have
                                orthographic projection.
            * show <bool> => Whether to show the plot or just hold it in RAM.
        """
        print("Can't visualise charge density with this software as mayavi won't install via pip!")

        # if ax is False:
        #     f = plt.figure()
        #     ax = f.add_subplot(111, projection="3d", proj_type="ortho")


        # if show:
        #     plt.show()

    ######################################
    ##             Parsing              ##
    ######################################

    def __get_units__(self):
        """Will parse the units from the first line and store them in the metadata."""
        line = [i.strip() for i in self._ltxt[0].split(' - ')]
        if len(line) == 2:
            self.metadata['gauss_file_type'], self.metadata['units'] = line

        else:
            raise SystemError("Parsing error in line 1.\n\n\tLine = '%s'" % self._ltxt[0]
                              + "\n\tError: Line split by ' - ' is not length 2.")

    def __get_charges_and_crds__(self):
        """
        Will parse the coordinate, charge and metadata info from the esp file.

        This will first get some metadata then call the function __get_atom_data__ which
        will parse the xyz data as well as the partial charge data.
        """
        for line in self._ltxt:
            if 'CHARGE' in line:
                line = self._ltxt[1]
                for i in line.split(' - '):
                    splitter = i.split('=')
                    if len(splitter) == 2:
                        self.metadata[splitter[0].strip().lower()] = type_check.eval_type(splitter[1].strip())

                self.__get_atom_data__()
                break

    def __get_atom_data__(self):
        """
        Will get the charge data and the coordinates from the xyz section.
        """
        for line_num, line in enumerate(self._ltxt):
            if 'ATOMIC COORDINATES AND ESP CHARGES.' in line:
                break

        # Get the number of atoms
        natoms = re.findall("#ATOMS *= *\d+", line)
        if len(natoms) == 1: natoms = natoms[0]
        else:
            raise SystemError("Parsing Error: Can't find number of atoms.")
        self.metadata['number_atoms'] = int(re.findall("\d+", natoms)[0])

        # Parse the lines
        cols, xyz, partial_charges = [], [], []
        for line_num1, line in enumerate(self._ltxt[line_num + 1: line_num + self.metadata['number_atoms'] + 1]):
            splitter = line.split()
            if len(splitter) != 5:
                args = (line_num + 2 + line_num1, len(splitter), line)
                raise SystemError("Parsing Error: Can't parse the coordinate and partial charge data."
                                  + "\n\n\tBad Line = %i\n\n\tNum Cols = %i\n\n\tBad Line = '%s'" % args)

            cols.append(splitter[0])
            xyz.append([self._double_float(i) for i in splitter[1:4]])
            partial_charges.append(self._double_float(splitter[4]))

        self.xyz_data = np.array(xyz)
        self.cols = np.array(cols)
        self.partial_charges = np.array(partial_charges)

    def __get_dipole_moment__(self):
        """Will parse the lines containing dipole moment info and store it in the metadata."""
        for line_num, line in enumerate(self._ltxt):
            if 'DIPOLE MOMENT' in line:
                break

        line = self._ltxt[line_num + 1]
        vals = re.findall("[a-zA-Z]+ *= *[-D0-9.]+", line)
        for i in vals:
            splitter = i.split('=')
            if len(splitter) == 2:
                self.metadata[f"dipole_{splitter[0].lower()}"] = self._double_float(splitter[1])

    def __get_quadrupole_moment__(self):
        """Will parse the lines containing dipole moment info and store it in the metadata."""
        for line_num, line in enumerate(self._ltxt):
            if 'QUADRUPOLE MOMENT' in line:
                break

        line = '  '.join(self._ltxt[line_num+1 : line_num+3])

        vals = re.findall("[a-zA-Z]+ *= *[-+D0-9.]+", line)
        vals = {i.split('=')[0].strip().lower(): self._double_float(i.split('=')[1]) for i in vals}

        self.metadata['quadrupole_moment'] = np.array([[vals['xx'], vals['xy'], vals['xz']],
                                                       [vals['xy'], vals['yy'], vals['yz']],
                                                       [vals['xz'], vals['yz'], vals['zz']]])

    def __get_charge_density__(self):
        """
        Will loop over the large block in the esp file and get the charge density data.

        In each row of the block there are 4 values:  x, y, z, val
            * x is x coord
            * y is y coord
            * z is z coord
            * val is the density at that point
        """
        for line_num, line in enumerate(self._ltxt):
            if 'ESP VALUES AND GRID' in line:
                break

        # Get the number of grid points
        npoints = re.findall("#POINTS *= *\d+", line)
        if len(npoints) == 1: npoints = npoints[0]
        else:
            raise SystemError("Parsing Error: Can't find number of grid points.")
        self.metadata['number_grid_points'] = int(re.findall("\d+", npoints)[0])

        xyz, vals = [], []
        for line_num1, line in enumerate(self._ltxt[line_num+1: line_num+1+self.metadata['number_grid_points']]):
            splitter = line.split()
            if len(splitter) != 4:
                args = (line_num + 2 + line_num1, len(splitter), line)
                raise SystemError("Parsing Error: Can't parse the charge density data."
                                  + "\n\n\tBad Line = %i\n\n\tNum Cols = %i\n\n\tBad Line = '%s'" % args)

            xyz.append([self._double_float(i) for i in splitter[:3]])
            vals.append(self._double_float(splitter[3]))

        self.grid_points = np.array(xyz)
        self.charge_density = np.array(vals)
