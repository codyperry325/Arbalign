import ase
from ase.io import read,write
import numpy as np

from scipy.optimize import linear_sum_assignment


reflects = []
reflects.append(np.array([-1,0,0,0,-1,0,0,0,-1]).reshape(3,3))
reflects.append(np.array([1,0,0,0,-1,0,0,0,-1]).reshape(3,3))
reflects.append(np.array([-1,0,0,0,1,0,0,0,-1]).reshape(3,3))
reflects.append(np.array([-1,0,0,0,-1,0,0,0,1]).reshape(3,3))
reflects.append(np.array([-1,0,0,0,1,0,0,0,1]).reshape(3,3))
reflects.append(np.array([1,0,0,0,-1,0,0,0,1]).reshape(3,3))
reflects.append(np.array([1,0,0,0,1,0,0,0,-1]).reshape(3,3))
reflects.append(np.eye(3,3))

swaps = []
swaps.append([0,1,2])
swaps.append([0,2,1])
swaps.append([2,1,0])
swaps.append([1,0,2])
swaps.append([2,0,1])
swaps.append([1,2,0])



def rmsd(V, W):
    """
    V - set of coordinates
    W - set of coordinates
    Returns root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

def kabsch(A, B):
    """
    Kabsch Algorithm as implemented by Jimmy Charnley Kromann
    Calculate RMSD between two XYZ files
    by: Jimmy Charnley Kromann <jimmy@charnley.dk> and
    Lars Andersen Bratholm <larsbratholm@gmail.com>
    project: https://github.com/charnley/rmsd
    license: https://github.com/charnley/rmsd/blob/master/LICENSE
    A - set of coordinates
    B - set of coordinates
    Performs the kabsch algorithm to calculate the RMSD between A and B
    Returns an RMSD
    """
    A_new = np.array(A)
    A_new = A_new - sum(A_new) / len(A_new)
    A = A_new
    B_new = np.array(B)
    B_new = B_new - sum(B_new) / len(B_new)
    B = B_new

    # Compute covariance matrix
    C = np.dot(np.transpose(A), B)

    # Compute singular value decomposition (SVD)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Compute rotation matrix
    U = np.dot(V, W)

    # Rotate A
    A = np.dot(A, U)

    # print(U)
    err =  rmsd(A, B)
    return err , U


def get_distance(c1, c2):
    x = (c1[0] - c2[0]) ** 2
    y = (c1[1] - c2[1]) ** 2
    z = (c1[2] - c2[2]) ** 2
    D = (x + y + z) ** 0.5
    return D

def align_axes(atoms,m_i):
    newcoords = []
    coords = atoms.get_positions()
    for c in coords:
        newcoords.append(np.matmul(m_i,c))
    atoms.set_positions(newcoords)

def sort_atoms(atoms):
    syms = atoms.get_chemical_symbols()
    uniq_elements = sorted(list(set(syms)))
    sorted_ids = []
    for e in uniq_elements:
        element_list = []
        for i in range(len(syms)):
            if syms[i] == e:
                element_list.append(i)
        sorted_ids.append(element_list)
    return sorted_ids

def rotops(theta, ID):
    if ID == 'Rx':
        return np.array([1,0,0,0,np.cos(theta),-np.sin(theta),0,np.sin(theta),np.cos(theta)]).reshape(3,3)
    elif ID == 'Ry':
        return np.array([np.cos(theta),0,np.sin(theta),0,1,0,-np.sin(theta),0,np.cos(theta)]).reshape(3,3)
    elif ID == 'Rz':
        return np.array([np.cos(theta),-np.sin(theta),0,np.sin(theta),np.cos(theta),0,0,0,1]).reshape(3,3)

def rotate_opt(monB,op,theta_deg):

    bcoords = monB.get_positions()
    bsyms = monB.get_chemical_symbols()


    rotated = []
    theta = np.pi * theta_deg / 180
    R_i = rotops(theta, op)
    for c in range(len(bcoords)):
        rot = np.matmul(R_i, bcoords[c])
        rotated.append(rot)

    rotated_B = ase.Atoms(positions=rotated,symbols=bsyms)
    return rotated_B

class CostMatrix():
    def __init__(self,A,B):
        self.monA = A
        self.monB = B

        self.sorted_A = self.monA[self.monA.numbers.argsort()]
        self.sorted_B = self.monB[self.monB.numbers.argsort()]

    def reorder(self):
        all_swapped_b = []
        all_b_sym = []
        acoords = self.monA.get_positions()
        bcoords = self.monB.get_positions()
        sorted_A_ids = sort_atoms(self.monA)
        sorted_B_ids = sort_atoms(self.monB)
        for i in range(len(sorted_B_ids)):
            assert len(sorted_A_ids[i]) == len(sorted_B_ids[i]), "Natoms in test and base structures not equal."
            self.A_coords = [acoords[x] for x in sorted_A_ids[i]]
            self.B_coords = [bcoords[x] for x in sorted_B_ids[i]]
            B_syms = [self.sorted_B.get_chemical_symbols()[x] for x in sorted_B_ids[i]]
            self.distance()
            B_reordered = [self.B_coords[x] for x in self.LSA[1]]
            for b in range(len(B_reordered)):
                all_swapped_b.append(B_reordered[b])
                all_b_sym.append(B_syms[b])
        Bfinal = ase.Atoms(positions=all_swapped_b,symbols=all_b_sym)
        return Bfinal

    def distance(self):
        D_ij = []
        for a in self.A_coords:
            for b in self.B_coords:
                D_ij.append(get_distance(a,b))
        D_ij = np.array(D_ij)
        D_ij = D_ij.reshape(len(self.A_coords),len(self.B_coords))
        self.LSA = linear_sum_assignment(D_ij)

class Structure():
    def __init__(self,lat,ss,packing,ucname='unitcell.cif'):
        self.lattice = lat
        self.ss = read(ss)
        self.packing = packing
        self.name = ucname

        self.packed = np.zeros(3)
        self.packed[0] = self.lattice[0] * self.packing[0]
        self.packed[1] = self.lattice[1] * self.packing[1]
        self.packed[2] = self.lattice[2] * self.packing[2]

    def supercell_xyz(self):
            a = self.packed[0]
            b = self.packed[1]
            c = self.packed[2]
            alpha = self.lattice[3] * np.pi / 180
            beta = self.lattice[4] * np.pi / 180
            gamma = self.lattice[5] * np.pi / 180
            n = (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)
            M = np.array([a, 0, 0, b * np.cos(gamma), b * np.sin(gamma), 0, c * np.cos(beta), c * n, c * np.sqrt(np.sin(beta) ** 2 - n ** 2)])
            M = M.reshape(3, 3)
            return (M)

    def get_distance(self,c1,c2):
        x = (c1[0] - c2[0])**2
        y = (c1[1] - c2[1])**2
        z = (c1[2] - c2[2])**2
        D = (x + y + z)**0.5
        return D

    def make_unitcell(self):
        allcoords = self.ss.get_positions()
        allatoms = self.ss.get_chemical_symbols()

        #Supercell sometimes necessary for pdim replacement!
        M = self.supercell_xyz()
        UC_xyz = []

        with open(self.name,'w') as f:
            f.write('data_Unitcell\n\n')
            f.write('_symmetry_space_group_name_H-M "P1"\n')
            f.write('_symmetry_Int_Tables_number  1\n')
            f.write('_symmetry_cell_setting TRICLINIC\n')
            f.write('_cell_length_a {:f}\n'.format(self.packed[0]))
            f.write('_cell_length_b {:f}\n'.format(self.packed[1]))
            f.write('_cell_length_c {:f}\n'.format(self.packed[2]))
            f.write('_cell_angle_alpha {:f}\n'.format(self.lattice[3]))
            f.write('_cell_angle_beta {:f}\n'.format(self.lattice[4]))
            f.write('_cell_angle_gamma {:f}\n\n'.format(self.lattice[5]))
            f.write('loop_\n')
            f.write('\t_symmetry_equiv_pos_site_id\n')
            f.write('\t_symmetry_equiv_pos_as_xyz\n')
            f.write('\t(X,Y,Z)\n\n')
            f.write('loop_\n')
            f.write('\t_atom_site_label\n')
            f.write('\t_atom_site_fract_x\n')
            f.write('\t_atom_site_fract_y\n')
            f.write('\t_atom_site_fract_z\n')

            for j in range(len(allcoords)):
                skip = False
                abc = np.matmul(allcoords[j],np.linalg.inv(M))
                for i in range(len(abc)):
                    k = (abc[i])
                    if k < 0:
                        t = k
                        while t < 0:
                            t += 1
                        loc = np.where(abc==k)[0][0]
                        abc[loc] = t
                    if k > 1:
                        t = k
                        while t > 1:
                            t -= 1
                        loc = np.where(abc==k)[0][0]
                        abc[loc] = t
                xyz = np.matmul(abc,M)
                UC_xyz.append(xyz)

                for i,c1 in enumerate(UC_xyz):
                    if i != j:
                        D = self.get_distance(xyz, c1)
                        if D < 0.6:
                            skip = True
                            break

                if not skip:
                    atom = allatoms[j]
                    [a,b,c]= abc
                    f.write('{:s}\t{:f}\t{:f}\t{:f}\n'.format(atom,a,b,c))

def align_product(monA,monB):
    noH = True

    m_a = monA.get_moments_of_inertia(vectors=True)[1]
    m_b = monB.get_moments_of_inertia(vectors=True)[1]
    align_axes(monA,m_a)
    align_axes(monB,m_b) #May be better with m_a alignment?

    sorted_A = monA[monA.numbers.argsort()]
    sorted_A_ids = sort_atoms(sorted_A)

    #get sorted base
    sorted_c = []
    sorted_s = []
    for i in range(len(sorted_A_ids)):
        A_coords = [sorted_A.get_positions()[x] for x in sorted_A_ids[i]]
        A_syms = [sorted_A.get_chemical_symbols()[x] for x in sorted_A_ids[i]]
        for x in range(len(A_coords)):
            sorted_s.append(A_syms[x])
            sorted_c.append(A_coords[x])
    monA = ase.Atoms(positions=sorted_c, symbols=sorted_s)

    write('aligned_B.xyz',monB)
    # write('aligned_A.xyz',monA)

    allstructures = []
    count = 0
    for r in reflects:
        for s in swaps:
            base = read('aligned_B.xyz')
            newcoords = []
            for coord in base.get_positions():
                line = np.zeros(3)
                for index,val in enumerate(s):
                    line[index] = coord[val]
                line = np.matmul(r,line)
                newcoords.append(line)
            base.set_positions(newcoords)
            allstructures.append(base)
            count += 1

    """
    Goal here is to loop through each symmetry operated dimer in allstructures and evaluate the distance between typed atoms.
    These distances will be fed into scipy.optimize.linear_sum_assignment (Hungarian algorithm) to find the optimal matches.
    Each distance is a 'cost matrix' - returns the idea arrangement. Would need to arrange EACH, calc RMSD, and choose the 
    lowest RMSD among each operation and ordering. 
    
    Methyl is problematic. ordering is correct but inverted for 9MA. FIX - NO H!!
    Toggle H on and off - code this!! 
    """

    count = 0
    rmse_list = []
    rot_list = []
    ordered_strucs = []
    #print(uniq_elements)
    for A in range(len(allstructures)):
        sorted_B = allstructures[A][allstructures[A].numbers.argsort()]

        cost = CostMatrix(sorted_A,sorted_B)
        monB = cost.reorder()

        Acoords = monA.get_positions()
        Bcoords = monB.get_positions()
        Asyms = monA.get_chemical_symbols()
        Bsyms = monB.get_chemical_symbols()
        [err,rot] = kabsch(Acoords,Bcoords)
        rotated = []
        for c in Bcoords:
            rotated.append(np.matmul(rot,c))

        monB = ase.Atoms(positions=Bcoords,symbols=Bsyms)

        if noH == True:
            hloc = sorted(list(set(monB.get_chemical_symbols()))).index('H')
            atomIDs = sort_atoms(monB)
            atomIDs.pop(hloc)
            A_noH = []
            B_noH = []
            for i in range(len(atomIDs)):
                A_noH.append(monA.get_positions()[atomIDs[i]])
                B_noH.append(monB.get_positions()[atomIDs[i]])
            [err, rot] = kabsch(A_noH[0],B_noH[0])
        else:
            [err, rot] = kabsch(monA.get_positions(), rotated)

        rmse_list.append(err)
        rot_list.append(rot)
        ordered_strucs.append(monB)
        count += 1

    minstruc = rmse_list.index(min(rmse_list))
    minrot = rot_list[minstruc]
    mincoords = ordered_strucs[minstruc].get_positions()
    rotated = []
    for c in mincoords:
        rotated.append(np.matmul(minrot,c))
    rotated_min = ase.Atoms(positions=rotated, symbols=ordered_strucs[minstruc].get_chemical_symbols())

    theta_range = np.arange(-20,30,2)
    allrotops = ['Rx','Ry','Rz']
    for op in allrotops:
        rmses = []
        for theta in theta_range:
            rotated_B = rotate_opt(rotated_min,op,theta)
            cost = CostMatrix(monA, rotated_B)
            monB = cost.reorder()
            newBcoords = monB.get_positions()
            [err, rot] = kabsch(monA.get_positions(), newBcoords)

            rmses.append(err)

        theta_min = theta_range[rmses.index(min(rmses))]
        rotated_B = rotate_opt(rotated_min,op,theta_min).get_positions()
        rotated_min.set_positions(rotated_B)

        cost = CostMatrix(monA, rotated_min)
        monB = cost.reorder()
        newBcoords = monB.get_positions()

        [err, rot] = kabsch(monA.get_positions(), newBcoords)
    test_A = monA.get_positions()
    test_B = monB.get_positions()
    [err, rot] = kabsch(test_A, test_B)
    rotated = []
    for c in test_B:
        rotated.append(np.matmul(rot, c))
    rotated_min = ase.Atoms(positions=rotated, symbols=ordered_strucs[minstruc].get_chemical_symbols())

    align_axes(rotated_min, m_a.T)
    return rotated_min


