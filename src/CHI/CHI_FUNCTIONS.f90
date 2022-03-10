MODULE ED_CHI_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER
  USE SF_INTEGRATE, only: trapz,simps
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
    integer                      :: iorb,iimp,unit,io,isite
    real(8),dimension(Ns_imp)    :: ChiT
    real(8) :: beta
    !
    call allocate_grids
    !
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    !
    call eval_chi_spin()
    if(MPIMASTER)then
       call ed_print_impChi()
       chiT= 0d0
       beta= 1d0/temp
       if(KondoFlag)then
          do iimp=1,Nimp
             io = Ns+iimp
             if(.not.chispin_flag(io))cycle
             chiT(io) = trapz(spinChi_tau(io,io,0:),0d0,beta)
          enddo
          unit = fopen("chiT_imp.ed",append=.true.)
          write(unit,*)temp,(chiT(Ns+io),io=1,Nimp)
          close(unit)
       endif
       if(any([chispin_flag(1:Norb)]))then
          do iorb=1,Norb
             if(.not.chispin_flag(iorb))cycle
             do isite=1,Nsites(iorb)
                io  = pack_indices(isite,iorb)
                chiT(io) = trapz(spinChi_tau(io,io,:),0d0,beta)
                print*,io,ChiT(io)
             enddo
          enddo
          unit = fopen("chiT.ed",append=.true.)
          write(unit,*)temp,(chiT(io),io=1,Ns_imp)
          close(unit)
       endif
    endif
    !
    call deallocate_grids
    !
  end subroutine eval_chi_lattice



end MODULE ED_CHI_FUNCTIONS
