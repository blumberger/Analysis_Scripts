"""
A script from Sam to add restraints to molecules
"""

natom = 36
at_per_mol = 36
mol2_file = "pos-init.mol2"

# put number of atoms per mol and index atoms to restrain for the single molecule (zero-based).
atoms_list_of_list_dict = { str(at_per_mol) : [[23, 22, 21, 20, 19],
                                          [24, 25, 26, 27, 10],
                                          [28, 9,  8,  7,  6],
                                          [5,  4,  3,  2,  1]]
                          }


nmolecules = natom // at_per_mol
datadict = "./"
mol2 = datadict + mol2_file


# coding: utf-8

# In[1]:

import numpy as np
import sys

masses = {
    'x': 0.000,
    'h': 1.008,
    'he': 4.003,
    'li': 7.016,
    'be': 9.012,
    'b': 11.009,
    'c': 12.000,
    'n': 14.003,
    'o': 15.995,
    'f': 18.998,
    'ne': 19.992,
    'na': 22.990,
    'mg': 23.985,
    'al': 26.981,
    'si': 27.977,
    'p': 30.974,
    's': 31.972,
    'cl': 34.969,
    'ar': 39.962,
    'k': 38.964,
    'ca': 39.963,
    'sc': 44.956,
    'ti': 47.948,
    'v': 50.944,
    'cr': 51.941,
    'mn': 54.938,
    'fe': 55.935,
    'co': 58.933,
    'ni': 57.935,
    'cu': 62.930,
    'zn': 63.929,
    'ga': 68.926,
    'ge': 73.921,
    'as': 74.922,
    'se': 79.917,
    'br': 78.918,
    'kr': 83.912,
    'rb': 84.912,
    'sr': 87.906,
    'y': 88.906,
    'zr': 89.905,
    'nb': 92.906,
    'mo': 97.905,
    'tc': 98.906,
    'ru': 101.904,
    'rh': 102.906,
    'pd': 107.904,
    'ag': 106.905,
    'cd': 113.903,
    'in': 114.904,
    'sn': 119.902,
    'sb': 120.904,
    'te': 129.906,
    'i': 126.904,
    'xe': 131.904
}

def print_(natoms, comment, extras, names, coordinates):
    #coordinates must be in array form
    print "%s" % natoms
    print "%s" % comment

    for i in range(natoms):
        extra = "   ".join(["%14.10f" % float(x) for x in extras[i]]) if extras else ""
        print "%3s   %14.10f   %14.10f   %14.10f   %s" % (names[i], coordinates[i, 0],
            coordinates[i, 1], coordinates[i, 2], extra)

def read_xyz(filename):
    '''reads xyz file and returns Geometry object'''

    #there might be multiple steps in the xyz files and this is storing all of them
    out = []

    with open(filename, "r") as f:
        line = f.readline()
        while line != "":
            natoms = int(line)
            comment = f.readline().rstrip()

            names = []
            coords = []
            extras = []

            for i in range(natoms):
                line = f.readline()
                data = line.split()
                name, x, y, z = data[0:4]
                extra = data[4:]

                names.append(name.capitalize())
                coords.append( [float(x), float(y), float(z)] )
                if extra:
                    extras.append(extra)
            out.append([names, np.array(coords), comment, extras,natoms])

            line = f.readline()

    return out

def calc_com_for_single_list(residue, single_list):
    """residue: is a list of list for a iven residues with the coordinates
       singles_list: is a list with the coordinates to calculate the COM for.
                if empty then it calculates the com for the whole mol"""

    #start center of mass calc
    numx = 0.0
    numy = 0.0
    numz = 0.0
    mass_tot = 0.0
    index_atoms_taken = []
    #loop over atoms of a residue
    if single_list == []:
        for index in range(len(residue[0])):#loop over the len of a single list
            mass = residue[index][0]
            numx += residue[index][1]*mass
            numy += residue[index][2]*mass
            numz += residue[index][3]*mass
            mass_tot += mass
            index_taken = residue[index][4]
            index_atoms_taken.append(index_taken)
    else:
        for index in single_list:
            mass = residue[index][0]
            numx += residue[index][1]*mass
            numy += residue[index][2]*mass
            numz += residue[index][3]*mass
            mass_tot += mass
            index_taken = residue[index][4]
            index_atoms_taken.append(index_taken)
    coms = np.array([numx/mass_tot, numy/mass_tot, numz/mass_tot])
    return coms, index_atoms_taken

def computeCOM(names, coordinates, masses):
    '''
    Returns the center of mass of the geometry.
    '''
    com = np.array([])
    if (len(com) == 3):
        return self.com
    else:
        sums = np.zeros(3)
        totMass = 0.0
        for i, name in enumerate(names):
            sums[:] += masses[name.lower()]*coordinates[i, :]
            totMass += masses[name.lower()]

        sums[:] /= totMass

        com = sums
        return sums

def merge_xyz_together(filename, filename_to_merge, filename_out):
    out = read_xyz(filename)
    names, coordinates, comment, extras,natoms = out[0]
    center = computeCOM(names, coordinates, masses )
    print_(natoms, comment, extras, names, np.array(coordinates))

    out2 = read_xyz(filename_to_merge)
    names1, coordinates1, comment1, extras1,natoms1 = out2[0]
    center = computeCOM(names1, coordinates1, masses )
    print_(natoms1, comment1, extras1, names1, np.array(coordinates1))

    #merge two xz together
    merged_names = names + names1
    merged_natoms = natoms1 + natoms
    merged_coordinates =  np.concatenate((np.array(coordinates), np.array(coordinates1)), axis=0)
    comment = comment + comment1
    extras = extras + extras1

    print_(merged_natoms, comment, extras, merged_names, np.array(merged_coordinates))
    write_out(merged_natoms, comment, extras, merged_names, np.array(merged_coordinates), filename_out)

#write out coordinates for each frames in separate files
#it could be useful if you want selectr fewer frames from an MD
def write_out(natoms, comment, extras, names, coordinates, filename_out):
    out_file = open(filename_out, "w")
    out_file.write("%s \n" % natoms)
    out_file.write("%s \n" % comment)

    for i in range(natoms):
        extra = "   ".join(["%14.10f" % float(x) for x in extras[i]]) if extras else ""
        out_file.write("%3s   %14.10f   %14.10f   %14.10f   %s \n" % (names[i], coordinates[i, 0],
           coordinates[i, 1], coordinates[i, 2], extra))
    out_file.close()


def _costruct_chain(filename, outfile, shift_mol, nmol, second_phase=False):
    """Constrcut a chain of mol with a given padding (shift_mol)
       NOTE: you need to start the original mol from the middle and with predetrmined orientation"""

    #this would read all the snaphshot
    out = read_xyz(filename)
    #take just the first snap
    names, coords, comment, extras,natoms =  out[0]
    coms = computeCOM(names, coords, masses)

    new_names = []
    new_coord = []
    new_natoms = 0

    #new_coord = list(coords)
    #print new_coord

    for mol in range(nmol):
        new_names.extend(names)
        new_natoms += natoms
        # shift the various atoms along shift_mol and append
        shift_mol_use = np.array(shift_mol)*(mol)
        # shift second phase if present wrt the first chain
        if second_phase:
            shift_mol_use += np.array(shift_mol)*(3)

        for atom in coords:
            shifted_coord = np.sum([atom, np.array(shift_mol_use)], axis=0)
            new_coord.append(shifted_coord)

    #print len(new_names)
    print "N atoms:", new_natoms
    #print len(new_coord)
    write_out(new_natoms, comment, extras, new_names, np.array(new_coord), outfile)
    print "DONE"


# this function calculates the coms and store the coordinates of the atoms in different lists per molecules
def com_calc(filename, natoms, n_BN_atoms):
    """NB: this function works only for systems with one kind of molecules with same number of atoms per molecule
      disctonary masses"""
    dict_masses = {
        "H" : 1.00794,
        "C" : 12.0107,
        "CP" : 12.0107,
        "N"  : 14.0067,
        "O"  : 15.999,
        "S" : 32.065,
        "F" : 18.9984
    }

    coms = []
    mol_specific_list = []
    f = open(filename, "r")

    nmol = len(f.readlines()[2+n_BN_atoms:])/natoms

    print "NMOL IN 3 PENTACENE LAYERS", nmol
    f.seek(0)
    #loop over the whole number of molecules that are all the same
    for itera in range(0,nmol):
        numx = 0.0
        numy = 0.0
        numz = 0.0
        mass_tot = 0.0
        #iterate over different mols and store the
        list_specific_mol = []
        num_iter = 0
        for line in f.readlines()[2+n_BN_atoms+(natoms*itera):2+natoms+n_BN_atoms+ ((natoms)*itera)]:

            #index atoms 1-based
            num_iter += 1
            index_atom = n_BN_atoms+(natoms*itera) + num_iter

            line_list = line.split()
            try:
                mass = dict_masses[line_list[0]]
            except:
                msg = "Atom mass not found"
                sys.exit(msg)
            line_list[0] = mass
            line_list = [float(i) for i in line_list]
            numx += line_list[1]*mass
            numy += line_list[2]*mass
            numz += line_list[3]*mass
            mass_tot += mass

            #append the atom index at the end of the list
            line_list.append(index_atom)

            #append MODiFIED line_list
            list_specific_mol.append(line_list)

        #center of mass append
        coms.append([numx/mass_tot, numy/mass_tot, numz/mass_tot])
        #store specific mol list
        mol_specific_list.append(list_specific_mol)


        # reinitialize the pointer otherwise it does not read following lines
        f.seek(0)
    #print "COMS is after applying com_calc:", coms
    return coms, mol_specific_list, nmol


#create three fakes layers with coms
def create_coms_external_file(filename, coms, layers_decomposition, total_n_equal_mol):
    print layers_decomposition
    print total_n_equal_mol
    coms_out = open(filename, "w")
    coms_out.write("%s \n title \n" %total_n_equal_mol)
    for num, point in enumerate(coms):
        if num < layers_decomposition[0]:
            coms_out.write("C \t %.10f \t %.10f \t %.10f \n" %(point[0], point[1], point[2]))
        elif layers_decomposition[0]<= num < layers_decomposition[0]+layers_decomposition[1]:
            coms_out.write("N \t %.10f \t %.10f \t %.10f \n" %(point[0], point[1], point[2]))
        else:
            coms_out.write("S \t %.10f \t %.10f \t %.10f \n" %(point[0], point[1], point[2]))

    coms_out.close()

#write out geo_center of the absolute restraint to check it
def create_res_coms_external_file(filename, coms,):

    coms_out = open(filename, "w")
    coms_out.write("%s \n title \n" %len(coms))
    for num, point in enumerate(coms):
        coms_out.write("C \t %.10f \t %.10f \t %.10f \n" %(point[0], point[1], point[2]))

def write_colvar_collective_files(collective_file, colvar_file, atoms_list_of_list_dict, list_mol_to_use, force_constant, target_absolute):
    out = open(collective_file, "w")
    out2 = open(colvar_file, "w")

    coms_restraints_list = []
    count = 0
    #loop over residues list _of_list which is the list of the residues
    for res in  list_mol_to_use:

        atoms_list_of_list = atoms_list_of_list_dict[str(len(res))]
        print atoms_list_of_list
        #loop over list with atom index for coms
        for single_list in atoms_list_of_list:
           coms_list_taken, index_atoms_taken = calc_com_for_single_list(res, single_list)
           count += 1

           index_atoms_taken = [str(i) for i in index_atoms_taken]
           list_ = " ".join(index_atoms_taken)

           coms_restraints_list.append(coms_list_taken)

           #convert to bhor the absolute point (coms)
           coms_list_taken = [str(float(i)*1.88973) for i in coms_list_taken]
           list_com = " ".join(coms_list_taken)


           print "COUNT:", count
           print "TARGET ABSOLUTE RESTR:", target_absolute
           print "ATOMS: ", index_atoms_taken
           print "COMS: ", coms_list_taken #bhor


           #here write COLLECTIVE.include
           out.write("""
   &COLLECTIVE
           COLVAR %s
           TARGET [angstrom] %s
           &RESTRAINT
                   K %s
           &END RESTRAINT
           INTERMOLECULAR
   &END COLLECTIVE\n\n
   """ %(count,str(target_absolute),str(force_constant)))

           #here write COLVAR.include
           out2.write("""
   &COLVAR
      &DISTANCE
           &POINT
                ATOMS %s
                TYPE GEO_CENTER
           &END POINT
           &POINT
               XYZ %s
               TYPE FIX_POINT
           &END POINT
           ATOMS 1 2
      &END DISTANCE
   &END COLVAR\n
   """ %(list_, list_com))

    out.close()
    out2.close()

    #return the COM of the restraints
    return coms_restraints_list


def mol_specific_list_from_mol2(mol2):
    """Construct a list of list [res 1: [ [coords], [coords]], res 2: [ [coords], [coords]] ... ]"""
    #read mol2
    mol2_list = []
    with open(mol2, "r") as file_:
        copy = False
        f = file_.readlines()
        # get the residues number from 3rd line
        number_of_resides =  f[2].split()[1]
        for line in f:
            splitter = line.split()
            if "@<TRIPOS>ATOM" in line:
                copy = True
                #continue
            elif "@<TRIPOS>BOND" in line:
                copy = False
                break
            elif copy:
                # take what we want: index mass coord res
                index = int(splitter[0])
                try:
                    mass = masses[splitter[1].lower()]
                except:
                    msg = "Atom mass not found"
                    sys.exit(msg)

                cx =  float(splitter[2])
                cy =  float(splitter[3])
                cz =  float(splitter[4])
                res = int(splitter[6])

                mol2_list.append([mass, cx, cy,cz, index, res])

    # get the correct mol specific list format
    mol_specific_list = []
    for res in range(int(number_of_resides)):
        res = res+1
        temp= []
        for list_ in  mol2_list:
            if list_[-1]  == res:
                temp.append(list_[:-1])
        mol_specific_list.append(temp)
    #print file_.read()
    print "DONE mol_specific_list"
    return mol_specific_list, mol2_list, number_of_resides


# add solvent around

def _get_os_dist( pos_solvent, coords):
    list = []
    for pos_atoms in coords:
        dist = np.linalg.norm(np.array(pos_atoms) - np.array(pos_solvent))
        list.append(dist)
    return min(list)

def _add_solvent(filename_out, _grid, _realgrid, _closest_dist):
    #self._get_grid()
    #self._get_carbon_pos(organic_crystal)

    #this would read all the snaphshot
    out = read_xyz(filename_out)
    #take just the first snap
    names, coords, comment, extras,natoms =  out[0]

    #list_of_list = [list(array(range(self._grid[i])) - int(self._grid[i] / 2)) for i in range(3)]
    list_of_list = [list(np.array(range(_grid[i])) - 1) for i in range(3)]
    print "list_of_list", list_of_list
    grid = [[i, j, k] for i in list_of_list[0]
            for j in list_of_list[1]
            for k in list_of_list[2]]
    print "grid", len(grid)
    solvent_crystal = []
    for pos in grid:
        realpos = np.array(pos) * np.array(_realgrid)
        #print realpos
        dist = _get_os_dist(realpos, coords)
        if dist >= _closest_dist:
            result = "%s  %f  %f  %f \n" % (_kind_solvent, realpos[0], realpos[1], realpos[2])
            solvent_crystal.append(result)
    return solvent_crystal


def _get_grid(_sizebox, _density):
    volume = np.prod(_sizebox)
    pre_natoms = _density * volume
    cubic_roots = int(np.power(pre_natoms, 1.0 / 3.0))
    print "CUBICCCC", cubic_roots
    _grid = [cubic_roots, cubic_roots, cubic_roots]
    _realgrid = np.array(_sizebox) / np.array(_grid)
    print "self._realgrid", _realgrid
    return _grid, _realgrid

def write_solv(solvent_file):
    with open(solvent_file, "w") as f:
        f.write("%s \n" %len(solvent_crystal))
        f.write("\n")
        for line in solvent_crystal:
            f.write(line)


# In[2]:
###################################################### MAIN ################################################################
# construct COLVAR files
# NOTE: you need the mol2 created by mercury for the whole system

#name output coms of the xyz
coms_out = datadict + "coms_slice.xyz"
#name output of the restraints position generated the xyz
rescom_out = datadict + "coms_of_restraints_slice.xyz"
#output file names
collective_file = datadict + "COLLECTIVE.include"
colvar_file = datadict + "COLVAR.include"

# force constant for the restraints
force_constant = 0.005  #force constant in au
# distance between the restraints and the target absolute position
target_absolute = 0.0 #target value in Ang

#list of molecules to apply restraint to. From start to end
# range starts from 0 (and ends to -1)
mol_to_take = range(0, nmolecules)


mol_specific_list, mol2_list, number_of_resides = mol_specific_list_from_mol2(mol2)


##############################################   START THE PROGRAM #########################################################

list_mol_to_use = [mol_specific_list[i] for i in mol_to_take]

coms_restraints_list = write_colvar_collective_files(collective_file, colvar_file, atoms_list_of_list_dict, list_mol_to_use, force_constant, target_absolute)
#print coms_restraints_list

create_res_coms_external_file(rescom_out, coms_restraints_list)
print "DONE writing coms restraints"


# In[4]:

# write out coms of all the molecules

tot_coms = []
for res in mol_specific_list:
    mass_tot = 0.0
    xm = 0.0
    ym = 0.0
    zm = 0.0
    for elem in res:
        mass_tot += elem[0]
        xm += elem[1]*elem[0]
        ym += elem[2]*elem[0]
        zm += elem[3]*elem[0]

    tot_coms.append([xm/mass_tot, ym/mass_tot, zm/mass_tot])


with open(coms_out, "w") as f:
    f.write("%s \n" %number_of_resides)
    f.write("\n")

    for coms in tot_coms:
        coms = map(str,coms)
        f.write("C %s %s %s \n" %(coms[0], coms[1], coms[2]))

print "DONE"


# In[ ]:



