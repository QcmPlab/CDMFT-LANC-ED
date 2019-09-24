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

An example, solving the Hubbard model on the Bethe lattice, is contained in the file `drivers/ced_hm_bethe.f90`.


--

***COPYRIGHT & LICENSING***  
Copyright 2012 -  (c), Adriano Amaricci, Lorenzo Crippa, Massimo Capone.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

