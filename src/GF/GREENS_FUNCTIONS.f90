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
  USE ED_GF_ELECTRONS
  USE ED_GF_IMPURITIES
  !
  implicit none
  private 

  public :: build_gf_lattice
  public :: eval_gf_lattice


contains


  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine build_gf_lattice()
    !
    write(LOGfile,"(A)")"Build lattice Greens functions:"
    if(any([gf_flag(1:Norb)]))call build_gf_electrons()
    if(KondoFlag)call build_gf_impurities()
    if(MPIMASTER)&
         call write_GFmatrix(impGmatrix,"gfmatrix"//str(ed_file_suffix)//".restart")
  end subroutine build_gf_lattice



  subroutine eval_gf_lattice()
    !
    call setup_gf
    !
    impGmats=zero
    impGreal=zero
    !
    write(LOGfile,"(A)")"Eval lattice Greens functions:"
    if(any([gf_flag(1:Norb)]))call eval_gf_electrons()
    if(KondoFlag)call eval_gf_impurities()
    if(MPIMASTER.AND.ed_print_G) call ed_print_impG()
    !
    call deallocate_grids
    !
  end subroutine eval_gf_lattice






  subroutine setup_gf
    call allocate_grids
    if(allocated(impGmats))deallocate(impGmats)
    if(allocated(impGreal))deallocate(impGreal)
    allocate(impGmats(Nspin,Ns,Ns,Lmats))
    allocate(impGreal(Nspin,Ns,Ns,Lreal))
    impGmats=zero
    impGreal=zero
  end subroutine setup_gf



end MODULE ED_GREENS_FUNCTIONS




