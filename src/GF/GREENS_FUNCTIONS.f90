MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_GF_ELECTRON
  !
  implicit none
  private 

  public :: build_gf_lattice
  logical,save                       :: iolegend=.true.
  real(8),dimension(:,:),allocatable :: zimp,simp
contains


  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine build_gf_lattice()
    !
    call allocate_grids
    !
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    call build_sigma_normal()
    !
    if(MPIMASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G) call ed_print_impG()
       if(ed_print_G0)call ed_print_impG0()
    endif
    !
    call build_szr
    !
    call deallocate_grids
    !
  end subroutine build_gf_lattice






  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine build_szr()
    integer :: ispin,is
    real(8) :: wm1,wm2
    integer :: unit
    integer :: iorb,jorb
    if(allocated(simp))deallocate(simp)
    if(allocated(zimp))deallocate(zimp)
    allocate(simp(Nspin,Ns),zimp(Nspin,Ns))
    wm1 = pi*temp ; wm2=3d0*pi*temp
    do ispin=1,Nspin
       do is=1,Ns
          simp(ispin,is) = dimag(impSmats(ispin,is,is,1)) - &
               wm1*(dimag(impSmats(ispin,is,is,2))-dimag(impSmats(ispin,is,is,1)))/(wm2-wm1)
          zimp(ispin,is)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,is,is,1))/wm1 ))
       enddo
    enddo

    unit = free_unit()
    open(unit,file="zeta_last.ed")
    write(unit,"(90(F15.9,1X))")&
         ((zimp(ispin,is),is=1,Ns),ispin=1,Nspin)
    close(unit)         

    unit = free_unit()
    open(unit,file="sig_last.ed")
    write(unit,"(90(F15.9,1X))")&
         ((simp(ispin,is),is=1,Ns),ispin=1,Nspin)
    close(unit)         

  end subroutine build_szr


end MODULE ED_GREENS_FUNCTIONS
