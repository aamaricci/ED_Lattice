MODULE ED_CHI_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_CHI_SPIN
  !
  implicit none
  private 

  public :: build_chi_lattice
  public :: eval_chi_lattice

contains


  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine build_chi_lattice()
    !
    call build_chi_spin()
    if(MPIMASTER)&
         call write_GFmatrix(SpinChiMatrix,"ChiSpinMatrix"//str(ed_file_suffix)//".restart")
    !
  end subroutine build_chi_lattice




  subroutine eval_chi_lattice()
    !
    call allocate_grids
    !
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    !
    call eval_chi_spin()
    if(MPIMASTER)call ed_print_impChi()
    !
    call deallocate_grids
    !
  end subroutine eval_chi_lattice



end MODULE ED_CHI_FUNCTIONS
