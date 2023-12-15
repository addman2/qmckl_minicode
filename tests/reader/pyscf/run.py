#!/usr/bin/env python
# coding: utf-8

# pySCF -> pyscf checkpoint file (NH3 molecule)

# load python packages
import os, sys

# load ASE modules
from ase.io import read

# load pyscf packages
from pyscf import gto, scf, mp, tools

#open boundary condition
checkpoint_file="H.chk"
pyscf_output="h"
charge=0
spin=1
basis="ccpvdz"
ecp='ccecp'
scf_method="DFT"  # HF or DFT
dft_xc="LDA_X,LDA_C_PZ" # XC for DFT


# read a structure
atom=read("molecule.xyz")
chemical_symbols=atom.get_chemical_symbols()
positions=atom.get_positions()
mol_string=""
for chemical_symbol, position in zip(chemical_symbols, positions):
    mol_string+="{:s} {:.10f} {:.10f} {:.10f} \n".format(chemical_symbol, position[0], position[1], position[2])

# build a molecule
mol = gto.Mole()
mol.atom = mol_string
mol.verbose = 10
mol.output = pyscf_output
mol.unit = 'A' # angstrom
mol.charge = charge
mol.spin = spin
mol.symmetry = True

# basis set
mol.basis = { "H": gto.basis.parse(
"""
H S
       1.000000  1.00000000
H P
       2.000000  1.00000000
""")}

# define ecp
mol.ecp = ecp

# molecular build
mol.build(cart=False)  # cart = False => use spherical basis!!

# calc type setting

if scf_method == "HF":
    # HF calculation
    if mol.spin == 0:
        print("HF kernel = RHF")
        mf = scf.RHF(mol)
        mf.chkfile = checkpoint_file
    else:
        print("HF kernel = ROHF")
        mf = scf.ROHF(mol)
        mf.chkfile = checkpoint_file

elif scf_method == "DFT":
    # DFT calculation
    if mol.spin == 0:
        print("DFT kernel = RKS")
        mf = scf.KS(mol).density_fit()
        mf.chkfile = checkpoint_file
    else:
        print("DFT kernel = ROKS")
        mf = scf.ROKS(mol)
        mf.chkfile = checkpoint_file
    mf.xc = dft_xc
else:
    raise NotImplementedError

total_energy = mf.kernel(max_cycle=1, tol=0.1)

# HF/DFT energy
print(f"Total HF/DFT energy = {total_energy}")
print("HF/DFT calculation is done.")
print("PySCF calculation is done.")
print(f"checkpoint file = {checkpoint_file}")
