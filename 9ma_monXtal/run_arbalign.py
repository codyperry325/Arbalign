import os
from ase.io import read,write

import genDimers
import arbalign

SSdimers = []
length = 0

Z = 4
job = 'dimer' #dimer pair -> product or monomer -> product
packing_start = [xxx,yyy,zzz] #SUpercell packing for monomer pair start
ciffile = genDimers.ReadData('9ma_monXtal.cif')
product = 'gas.xyz'
fname = 'FIgas.xyz'
reactant_path = os.path.join(os.getcwd(),fname)

cifdata = ciffile.cif()
supercell_set = genDimers.Setup(cifdata,Z,job=job,dim=[xxx,yyy,zzz],fname=fname)
SSdata = supercell_set.supercell()
lattparams = cifdata[1]
monB = read(product)
centered = monB.get_positions() - monB.get_center_of_mass()
monB.set_positions(centered)


for mfile in os.listdir(reactant_path):
    tag = mfile.replace('.xyz','')
    monA = read(os.path.join(reactant_path,mfile))
    monA_com = monA.get_center_of_mass()
    centered = monA.get_positions() - monA_com
    monA.set_positions(centered)
    monBaligned = arbalign.align_product(monA,monB)
    monBaligned_coords = monBaligned.get_positions() + monA_com
    monBaligned.set_positions(monBaligned_coords)
    # write(tag + '-aligned.xyz',monBaligned)
    SSdimers.append(monBaligned)
    length += len(monBaligned_coords)
#
replaced_ss = '{:s}_replaced_supercell.xyz'.format(fname)
with open(replaced_ss,'w') as f:
    f.write('{:d}\n\n'.format(length))
    for dim in SSdimers:
        coords = dim.get_positions()
        atoms = dim.get_chemical_symbols()
        for k,c in enumerate(coords):
            f.write('{:s}\t{:f}\t{:f}\t{:f}\n'.format(atoms[k],c[0],c[1],c[2]))

unitcell = arbalign.Structure(lattparams,replaced_ss,packing_start,ucname='{:s}.cif'.format(fname))
unitcell.make_unitcell()

