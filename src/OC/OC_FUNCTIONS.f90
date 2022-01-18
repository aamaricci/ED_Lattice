MODULE ED_OC_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,splot
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
  USE ED_OC_ELECTRONS
  !
  implicit none
  private 

  public :: build_oc_lattice

contains


  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine build_oc_lattice()
    integer :: iorb
    real(8) :: K
    call allocate_grids
    !
    !BUILD OC
    Drude_weight = zero
    OptCond_w    = zero
    call build_oc_electrons()
    !
    K = -ed_Ekin
    do iorb=1,Norb
       Drude_weight(iorb) = -pi*(K + 2/Nsites(iorb)*Drude_weight(iorb))
    enddo
    !
    if(MPIMASTER)call ed_print_oc()
    !
    call deallocate_grids
    !
  end subroutine build_oc_lattice




  subroutine ed_print_oc
    integer           :: iorb,unit
    character(len=32) :: suffix
    !
    do iorb=1,Norb
       suffix="OptCond"//"_l"//str(iorb)
       call splot(reg(suffix)//"_realw.ed",vr,OptCond_w(iorb,:))
    enddo
    !
    unit = free_unit()
    open(unit,file="drude_weight.ed")
    write(unit,"(10(F15.9,1X))")(Drude_weight(iorb),iorb=1,Norb)
    close(unit)
    !
  end subroutine ed_print_oc


end MODULE ED_OC_FUNCTIONS
