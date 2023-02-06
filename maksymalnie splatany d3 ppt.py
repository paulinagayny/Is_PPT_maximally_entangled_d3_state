# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 01:24:52 2021

@author: Paulina Gayny
"""

from qutip import *
import numpy as np
#import matplotlib.pyplot as plt

print(basis(3,0))
print(basis(3,1))
print(basis(3,2))

state00 = tensor(basis(3,0), basis(3,0))
state11 = tensor(basis(3,1), basis(3,1))
state22 = tensor(basis(3,2), basis(3,2))

state01 = tensor(basis(3,0), basis(3,1))
state02 = tensor(basis(3,0), basis(3,2))
state10 = tensor(basis(3,1), basis(3,0))
#state11
state12 = tensor(basis(3,1), basis(3,2))
state20 = tensor(basis(3,2), basis(3,0))
state21 = tensor(basis(3,2), basis(3,1))
#state22

#print(state00)
# print(state01)
# print(state02)
# print(state10)
#print(state11)
# print(state12)
# print(state20)
# print(state21)
#print(state22)

#baza
A = [state00, state01, state02, state10, state11, state12, state20, state21, state22]

psi = (state00 + state11 + state22).unit()
#psi_psi = tensor(psi, psi)

#psi tensor psi
#print(psi)
psi_psi = psi * psi.trans()
#print(psi_psi)

#identity natrix
ident_matrix = state00 * state00.trans()

for k in range(1, len(A)):
    state = A[k] * A[k].trans()
    ident_matrix = ident_matrix + state
        
#print(ident_matrix)


def Is_PPT():
    eps = 0.001
    p = 1
    while p >= 0:
        #density matrix of a mixed state
        Rho = p*psi_psi  + (1-p)*(1/9)*ident_matrix
        Rho_pt = partial_transpose(Rho, [1,0])
        
        p = p - eps
        
        Rho_next = p*psi_psi  + (1-p)*(1/9)*ident_matrix
        Rho_next_pt = partial_transpose(Rho_next, [1,0])
        
        eigenenergies_array = Rho_pt.eigenenergies()
        eigenenergies_array_next = Rho_next_pt.eigenenergies()
        #print(min(eigenenergies_array) < 0)
        #print(Rho)
        #print(p)
        
        if min(eigenenergies_array) < 0 and min(eigenenergies_array_next) > 0:
            print("not PPT for p > ", p)
            #print(Rho_pt)
            #print(eigenenergies_array)
            print("Lowest eigenvalue for p = ", p)
            print(min(eigenenergies_array))
        
            print("Lowest eigenvalue for p = ", p - eps)
            print(min(eigenenergies_array_next))
            return
        

Is_PPT()

#example
#print(psi_psi)
print("Eigenvalues for p = 0,23")
rho = 0.23*psi_psi + (0.77)*(1/9)*ident_matrix
rho_pt = partial_transpose(rho, [1,0])
#print(rho_pt)
print(rho_pt.eigenenergies())

print("Eigenvalues for p = 0,27")
#print(psi_psi)
rho = 0.27*psi_psi + (0.73)*(1/9)*ident_matrix
rho_pt = partial_transpose(rho, [1,0])
#print(rho_pt)
print(rho_pt.eigenenergies())
    
    