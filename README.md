# Cluster-DMFT Lanczos Exact Diagonalization

A Lanczos based solver for the **Cluster** Dynamical Mean-Field Theory using the N_up:N_dw implementation.  
**This code solves the normal (N_up, N_dw) case only.**

The code is based on:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

* DMFT_Tools [https://github.com/aamaricci/DMFTtools](https://github.com/aamaricci/DMFTtools) [For the driver part, see below]

The code structure is as follow:  

* The set of modules compile into a top layer named `DMFT_ED.f90`  
* The actual implementation of the DMFT equations is case by case performed in a driver program, usually placed in the directory `drivers`. 
* In the driver code the user must includes the `DMFT_ED` module and call the necessary procedures to solve the DMFT equations.

An example, solving the Hubbard model on the one-dimensional Hubbard chain, is contained in the file `drivers/ced_hm_1dchaom.f90`.

## DEVELOPMENT

### MILESTONE 1

- [x] Write code for the spin-decoupled case
- [x] Code compiles and runs without errors
- [x] Check exact diagonalization comparing results with known 2x2 plaquette data (without bath)
- [x] Check single-site cluster case against single-site DMFT code
- [ ] Check 1d Hubbard chain against literature

### MILESTONE 2

- [ ] Include complex routines for diagonalization (and bath?)
- [ ] Test 2d BHZ model 

### MILESTONE 3

- [ ] Test 3d BHZ model for non spin-coupling choices of cluster
- [ ] Add real-space CDMFT case for finite systems
- [ ] Rewrite the code for the general spin-coupled case


--

***COPYRIGHT & LICENSING***  
Copyright 2012 -  (c), Adriano Amaricci, Lorenzo Crippa, Massimo Capone.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

