subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use bl_error_module
  use bl_constants_module
  use riemann_module
  use meth_params_module
  
  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt)         :: problo(2), probhi(2)

  real(rt)        , allocatable :: ql(:,:,:), qr(:,:,:)
  real(rt)        , allocatable :: gamcl(:,:), gamcr(:,:)
  real(rt)        , allocatable :: cav(:,:), smallc(:,:)
  real(rt)        , allocatable :: uflx(:,:,:), qint(:,:,:)

  integer :: ilo, ihi, jlo, jhi
  integer :: idir
  
  integer :: untin,i

  namelist /fortin/ rho_l, u_l, p_l, re_l, gc_l, &
                    rho_r, u_r, p_r, re_r, gc_r, &
                    cav_s, smallc_s, idir

  
  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set namelist defaults

  rho_l = 1.0
  u_l = 1.0
  p_l = 1.0
  re_l = 1.0
  gc_l = 4.0/3.0

  rho_r = 1.0
  u_r = 1.0
  p_r = 1.0
  re_r = 1.0
  gc_r = 4.0e0_rt/3.0e0_rt
  
  cav_s = 1.0
  smallc_s = 1.e-10_rt
  
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

  ! call the Riemann solver
  ilo = 1
  ihi = 1
  jlo = 1
  jhi = 1

  ! this should take place in init_godunov_indices, but can't call
  ! that from here
  ngdnv = 6                                                                                
  GDRHO = 1                                                                                
  GDU = 2                                                                                  
  GDV = 3                                                                                  
  GDW = 4                                                                                  
  GDPRES = 5                                                                               
  GDGAME = 6            
  
  allocate(ql(ilo:ihi, jlo:jhi, QVAR))
  allocate(qr(ilo:ihi, jlo:jhi, QVAR))

  allocate(gamcl(ilo:ihi, jlo:jhi))
  allocate(gamcr(ilo:ihi, jlo:jhi))

  allocate(cav(ilo:ihi, jlo:jhi))  
  allocate(smallc(ilo:ihi, jlo:jhi))

  allocate(uflx(ilo:ihi, jlo:jhi, NVAR))
  allocate(qint(ilo:ihi, jlo:jhi, ngdnv))    ! should be ngdnv, but this is larger

  ! set the Riemann arrays
  ql(:,:,:) = ZERO
  ql(ilo:ihi, jlo:jhi, QRHO) = rho_l
  ql(ilo:ihi, jlo:jhi, QU) = u_l
  ql(ilo:ihi, jlo:jhi, QPRES) = p_l
  ql(ilo:ihi, jlo:jhi, QREINT) = re_l  

  gamcl(:,:) = ZERO
  gamcl(ilo:ihi, jlo:jhi) = gc_l
  
  qr(:,:,:) = ZERO
  qr(ilo:ihi, jlo:jhi, QRHO) = rho_r
  qr(ilo:ihi, jlo:jhi, QU) = u_r
  qr(ilo:ihi, jlo:jhi, QPRES) = p_r
  qr(ilo:ihi, jlo:jhi, QREINT) = re_r  

  gamcr(:,:) = ZERO
  gamcr(ilo:ihi, jlo:jhi) = gc_r

  cav(:,:) = ZERO
  cav(ilo:ihi, jlo:jhi) = cav_s

  smallc(:,:) = ZERO
  smallc(ilo:ihi, jlo:jhi) = smallc_s
  
  ! call the Riemann solver
  idir = 1

  call riemanncg(ql, qr, ilo, jlo, ihi, jhi,  &
                 gamcl, gamcr, cav, smallc, ilo, jlo, ihi, jhi, &
                 uflx, ilo, jlo, ihi, jhi, &
                 qint, ilo, jlo, ihi, jhi, &
                 idir, ilo, ihi-1, jlo, jhi, [ilo-1, jlo-1], [ihi+1, jhi+1])
  
  ! we're done -- abort the code
  call bl_error("done with Riemann")
  
end subroutine PROBINIT


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use probdata_module
  use eos_module, only : gamma_const
  use bl_constants_module, only: M_PI, FOUR3RD
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS
  use prob_params_module, only : center
  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  real(rt)         :: xlo(2), xhi(2), time, delta(2)
  real(rt)         :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  

end subroutine ca_initdata

