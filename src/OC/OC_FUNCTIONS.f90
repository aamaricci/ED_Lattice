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
  public :: eval_oc_lattice

contains


  !+------------------------------------------------------------------+
  ! OC DRUDE CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine build_oc_lattice()
    integer :: iorb,unit
    real(8) :: K,D(Norb)
    character(len=32) :: suffix
    !
    write(LOGfile,"(A)")"Build lattice OC Drude:"

    call build_oc_electrons()
    !
    if(MPIMASTER)&
         call write_GFmatrix(OcMatrix,"OcMatrix"//str(ed_file_suffix)//".restart")
  end subroutine build_oc_lattice




  subroutine eval_oc_lattice()
    integer :: iorb,unit
    real(8) :: K,D(Norb),Res(Norb),temp0=0d0
    character(len=32) :: suffix
    character(len=:),allocatable :: fname
    logical :: bool
    !
    call allocate_grids
    !
    Res          = 0d0
    Drude_weight = 0d0
    OptCond_w    = zero
    !
    write(LOGfile,"(A)")"Eval lattice OC Drude:"
    call eval_oc_electrons()
    !
    D = Drude_weight

    !
    do iorb=1,Norb
       if(.not.oc_flag(iorb))cycle
       Drude_weight(iorb) = -pi/Nsites(iorb)*(ed_Ekin + 2*D(iorb))
       Res(iorb) = -1d0/Drude_weight(iorb)
    enddo
    !
    if(MPIMASTER)then
       do iorb=1,Norb
          if(.not.oc_flag(iorb))cycle
          suffix="OptCond"//"_l"//str(iorb)
          call splot(reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,OptCond_w(iorb,:))
       enddo
       !
       unit = fopen("drude.ed",append=.true.)
       write(unit,*)(Drude_weight(iorb),iorb=1,Norb)
       close(unit)
       !
       unit = fopen("resistivity.ed",append=.true.)
       write(unit,*)(Res(iorb),iorb=1,Norb)
       close(unit)
       !
       unit = fopen("jjcorr.ed",append=.true.)
       write(unit,*)(D(iorb),iorb=1,Norb)
       close(unit)
    endif
    call deallocate_grids
    !
  end subroutine eval_oc_lattice



end MODULE ED_OC_FUNCTIONS
