# Cluster-DMFT Lanczos Exact Diagonalization

A Lanczos based solver for the **Cluster** Dynamical Mean-Field Theory using the N_up:N_dw implementation.  
**This code solves the normal (N_up, N_dw) case only.**

The code is based on:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

* DMFT_Tools [https://github.com/aamaricci/DMFTtools](https://github.com/aamaricci/DMFTtools) [For the driver part, see below]

The code structure is as follow:  

* The set of modules compile into a top layer named `CDMFT_ED.f90`  
* The actual implementation of the DMFT equations is case by case performed in a driver program, usually placed in the directory `drivers`. 
* In the driver code the user must includes the `CDMFT_ED` module and call the necessary procedures to solve the DMFT equations.

 The only bath type implemented to date is the 'replica' one, in which the bath consits of multiple copies of the original cluster hamiltonian. 
 The bath can be initialized in two ways: one can pass the hamiltonian of the cluster as a matrix or as a structure consisting of an array of bath parameters and an array of basis matrices coupling to each of the formers. This second way makes easier to adhere to the symmetries of the model.
 An example, solving the Hubbard model on the one-dimensional Hubbard chain, is contained in the file `drivers/cdn_hm_1dchain.f90`.

## DEVELOPMENT

### MILESTONE 1

- [x] Write code for the spin-decoupled case
- [x] Code compiles and runs without errors
- [x] Check exact diagonalization comparing results with known 2x2 plaquette data (without bath)
- [x] Check single-site cluster case against single-site DMFT code
- [x] Check 1d Hubbard chain against literature

### MILESTONE 2

- [x] Include complex routines for diagonalization and bath
- [ ] Test 2d BHZ model 

### MILESTONE 3

- [ ] Add real-space CDMFT case for finite systems
- [ ] Test Kane-Mele model with real-space CDMFT
- [ ] Test 3d BHZ model for non spin-coupling choices of cluster
- [ ] Rewrite the code for the general spin-coupled case


--

***COPYRIGHT & LICENSING***  
Copyright 2012 -  (c), Adriano Amaricci, Lorenzo Crippa, Massimo Capone.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

