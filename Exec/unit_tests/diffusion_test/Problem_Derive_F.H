
extern "C"
{

  void ca_deranalytic
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

}
