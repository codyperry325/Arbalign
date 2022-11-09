import numpy as np
import os
import ase
from ase.io import read,write
import shutil


class ReadData():
    def __init__(self,f):
        self.file = f

    def cif(self):
        coorlabel = False
        coords = []
        with open(self.file,'r') as fh:
            lines = fh.readlines()
            for line in enumerate(lines):
                if '_cell_length_a' in line[1]:
                    a = float(line[1].split()[1])
                elif '_cell_length_b' in line[1]:
                    b = float(line[1].split()[1])
                elif '_cell_length_c' in line[1]:
                    c = float(line[1].split()[1])
                elif '_cell_angle_alpha' in line[1]:
                    alpha = float(line[1].split()[1])
                elif '_cell_angle_beta' in line[1]:
                    beta = float(line[1].split()[1])
                elif '_cell_angle_gamma' in line[1]:
                    gamma = float(line[1].split()[1])
                elif '_atom_site_fract_z' in line[1]:
                    coorlabel = True
                    continue
                if coorlabel == True:
                    line = line[1].split()
                    line[1] = float(line[1])
                    line[2] = float(line[2])
                    line[3] = float(line[3])
                    coords.append(line)
        lattparams = [a,b,c,alpha,beta,gamma]
        #print(coords)
        return coords,lattparams

    def xyz(self):
        #print(self.file)
        coords = []
        with open(self.file,'r') as fh:
            lines = fh.readlines()
            coorRead = False
            for line in enumerate(lines):
                if 'Monomer' in line[1]:
                    coorRead = True
                    continue
                if coorRead == True:
                    coords.append(line[1])
        return coords

class Setup():
    def __init__(self,c,z,job='monomer',dim=[2,2,2],fname='smartalign_inps'):
        self.data = c
        self.Z = z
        self.jobtype = job
        self.ss_dim = dim
        self.fname = fname

        self.target = os.path.join(os.getcwd(),self.fname)


    def supercell(self):
        print("HERE")
        path = os.getcwd() + os.path.sep
        xyzpath = path + 'allxyz' + os.path.sep
        os.makedirs(xyzpath,exist_ok=True)
        coords = self.data[0]
        Trans = np.zeros(3)
        count = 0
        SScoords = []
        SScoords_xyz = []

        #3x3x3 Supercell
        #ixjxk
        for i in range(self.ss_dim[0]):
            for j in range(self.ss_dim[1]):
                for k in range(self.ss_dim[2]):
                    Trans[0] = i
                    Trans[1] = j
                    Trans[2] = k
                    #print(Trans)
                    for line in coords:
                        img = (line[1:4] + Trans)
                        imgline = [line[0],img[0],img[1],img[2]]
                        SScoords.append(imgline)
                    count +=1
        xyzCell = self.supercell_xyz()
        for line in SScoords:
            a = line[1]
            b = line[2]
            c = line[3]
            coord = [a,b,c]
            tmp = (np.matmul(coord,xyzCell))
            xyzline = [line[0],tmp[0],tmp[1],tmp[2]]
            SScoords_xyz.append(xyzline)
        with open('SS.xyz','w') as f:
            f.write('{:s}\n'.format(str(len(SScoords_xyz))))
            f.write('supercell\n')
            for line in SScoords_xyz:
                for item in line:
                    f.write('{:s}\t'.format(str(item)))
                f.write('\n')

        #get monomers
        natoms = len(coords) / int(self.Z)
        moncount = 0
        atomcount = 0
        monomer = []
        centers = []
        for i in range(len(SScoords_xyz)):
            atomcount += 1
            start_ind = int(moncount * natoms)
            index = start_ind + atomcount - 1
            monomer.append(SScoords_xyz[index])
            if np.mod(i + 1, natoms) == 0:
                syms = []
                moncoords = []
                fname = os.path.join(xyzpath,str(moncount + 1).zfill(3) + '.xyz')

                for line in monomer:
                    syms.append(line[0])
                    moncoords.append(line[1:])
                mon_ase = ase.Atoms(symbols=syms,positions=moncoords)
                write(fname,mon_ase)

                centers.append(self.COM(monomer))
                moncount += 1
                atomcount = 0
                monomer = []
        all_dcom = self.distances(centers)
        close = all_dcom[0]
        uniq = all_dcom[1]

        seen = set()
        if self.jobtype == 'dimer':
            os.makedirs(self.target, exist_ok=True)
            for ID in close:
                mcount = 0
                id1 = ID[0]
                id2 = ID[1]
                if id1 in seen or id2 in seen:
                    #print("Duplicate: ", id1, id2)
                    continue
                else:
                    seen.add(id1)
                    seen.add(id2)
                dfile = os.path.join(self.target,str(id1).zfill(3) + '-' + str(id2).zfill(3) + '.xyz')
                for entry in sorted(os.listdir(xyzpath)):
                    monnum = int(entry.strip('.xyz'))
                    if monnum == id1:
                        mon1 = read(os.path.join(xyzpath,entry))
                        a1 = mon1.get_chemical_symbols()
                        c1 = mon1.get_positions()
                        mcount += 1
                    if monnum == id2:
                        mon2 = read(os.path.join(xyzpath, entry))
                        a2 = mon2.get_chemical_symbols()
                        c2 = mon2.get_positions()
                        mcount += 1
                    if mcount == 2:
                        break

                dimer_syms = np.append(a1,a2,axis=0)
                dimer_coords = np.append(c1,c2,axis=0)
                reacted = ase.Atoms(symbols=dimer_syms,positions=dimer_coords)
                write(dfile,reacted)

        elif self.jobtype == 'monomer':
            shutil.copytree(xyzpath,self.target)
            print('MONOMER COUNT:', moncount)

        print(xyzpath)
        shutil.rmtree(xyzpath)

        return uniq,centers

    def supercell_xyz(self):
        lattice = self.data[1]
        a = lattice[0]
        b = lattice[1]
        c = lattice[2]
        alpha = lattice[3] * np.pi / 180
        beta = lattice[4] * np.pi / 180
        gamma = lattice[5] * np.pi / 180
        n = (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)
        M = np.array([a, 0, 0, b * np.cos(gamma), b * np.sin(gamma), 0, c * np.cos(beta), c * n, c * np.sqrt(np.sin(beta) ** 2 - n ** 2)])
        M = M.reshape(3, 3)
        return (M)

    def COM(self,coords):
        newcoords = []
        for line in coords:
            if line[0] == 'C':
                line[0] = 12.01
                newcoords.append(line)
            elif line[0] == 'H':
                line[0] = 1.01
                newcoords.append(line)
            elif line[0] == 'O':
                line[0] = 16.00
                newcoords.append(line)
            elif line[0] == 'N':
                line[0] = 14.01
                newcoords.append(line)
            elif line[0] == 'S':
                line[0] = 32.06
                newcoords.append(line)
            elif line[0] == 'F':
                line[0] = 19.
                newcoords.append(line)
            else:
                print('ATOM NOT CODED YET --',line[0])
        mass = 0
        vec = 0
        count = 0
        for line in newcoords:
            vec += np.multiply(line[0], line[1:4])
            mass += line[0]
            count += 1
        com = vec / mass
        return com

    def distances(self,centers):
        distances = []
        dimercom = []
        uniq = []
        mon1_id = 0
        for i in centers:
            mon1_id += 1
            mon2_id = 0
            for k in centers:
                mon2_id += 1
                if mon1_id != mon2_id:
                    dcom = i - k
                    dcom = (dcom[0]**2 + dcom[1]**2 + dcom[2]**2)**0.5
                    dcom = round(dcom,10)
                    if dcom < 5:
                        #print('DIMER??:', mon1_id,mon2_id)
                        dID = [mon1_id,mon2_id,dcom]
                        dimercom.append(dID)
                    if dcom not in distances and dcom < 18:
                        distances.append(dcom)
                        dID = [mon1_id,mon2_id,dcom]
                        uniq.append(dID)

        return dimercom,uniq
