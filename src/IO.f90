MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  implicit none
  private


  public :: ed_get_sigma_matsubara
  public :: ed_get_sigma_realaxis
  public :: ed_get_gimp_matsubara
  public :: ed_get_gimp_realaxis
  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_eimp
  public :: ed_get_doubles
  public :: ed_get_density_matrix


  !****************************************************************************************!
  !****************************************************************************************!

  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impD
  public :: ed_print_impChi


  !****************************************************************************************!
  !****************************************************************************************!


  character(len=64)                :: suffix

  !Retrieve imp GF through routines.
  interface ed_get_gimp_matsubara
     module procedure ed_get_gimp_matsubara_main
  end interface ed_get_gimp_matsubara


  interface ed_get_gimp_realaxis
     module procedure ed_get_gimp_realaxis_main
  end interface ed_get_gimp_realaxis


  !Retrieve self-energy through routines:
  interface ed_get_sigma_matsubara
     module procedure ed_get_sigma_matsubara_main
  end interface ed_get_sigma_matsubara

  interface ed_get_sigma_realaxis
     module procedure ed_get_sigma_realaxis_main
  end interface ed_get_sigma_realaxis


  !Retrieve static common observables  
  interface ed_get_dens
     module procedure ed_get_dens_main
  end interface ed_get_dens

  interface ed_get_mag
     module procedure ed_get_mag_main
  end interface ed_get_mag

  interface ed_get_docc
     module procedure ed_get_docc_main
  end interface ed_get_docc

  interface ed_get_eimp
     module procedure :: ed_get_eimp_main
  end interface ed_get_eimp

  interface ed_get_doubles
     module procedure :: ed_get_doubles_main
  end interface ed_get_doubles

  interface ed_get_density_matrix
     module procedure :: ed_get_density_matrix_single
  end interface ed_get_density_matrix





contains



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity GF and self-energy 
  !+--------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA GREEN'S FUNCTIONS
  subroutine ed_get_gimp_matsubara_main(Gmats)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats(:,:,:,:,:) = impGmats(:,:,:,:,:)
  end subroutine ed_get_gimp_matsubara_main

  !NORMAL, REAL SELF-ENERGY
  subroutine ed_get_gimp_realaxis_main(Greal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal(:,:,:,:,:) = impGreal(:,:,:,:,:)
  end subroutine ed_get_gimp_realaxis_main


  !NORMAL, MATSUBARA SELF-ENEGRGY
  subroutine ed_get_sigma_matsubara_main(Smats)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats(:,:,:,:,:) = impSmats(:,:,:,:,:)
  end subroutine ed_get_sigma_matsubara_main


  !NORMAL, REAL SELF-ENERGY
  subroutine ed_get_sigma_realaxis_main(Sreal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal(:,:,:,:,:) = impSreal(:,:,:,:,:)
  end subroutine ed_get_sigma_realaxis_main



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+--------------------------------------------------------------------------+!
  subroutine ed_get_dens_main(dens)
    real(8),dimension(Norb) :: dens
    dens = ed_dens
  end subroutine ed_get_dens_main

  subroutine ed_get_mag_main(mag) 
    real(8),dimension(Norb) :: mag
    mag = (ed_dens_up-ed_dens_dw)
  end subroutine ed_get_mag_main

  subroutine ed_get_docc_main(docc) 
    real(8),dimension(Norb) :: docc
    docc = ed_docc
  end subroutine ed_get_docc_main

  subroutine ed_get_doubles_main(docc)
    real(8),dimension(4) :: docc
    docc = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
  end subroutine ed_get_doubles_main


  subroutine ed_get_eimp_main(eimp)
    real(8),dimension(4) :: eimp
    eimp = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
  end subroutine ed_get_eimp_main





  subroutine ed_get_density_matrix_single(dm_)
    implicit none
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:)
    !internal
    integer                                      :: unit
    integer                                      :: iorb,jorb,ispin,jspin,io,jo
    complex(8)                                   :: Tr
    complex(8),allocatable                       :: dm_custom_rot(:,:)
    real(8)                                      :: soc
    !
    if (.not.allocated(imp_density_matrix)) then
       write(LOGfile,"(A)") "imp_density_matrix is not allocated"
       stop
    endif
    !
    if(allocated(dm_))deallocate(dm_)
    allocate(dm_(Nspin*Norb,Nspin*Norb))
    dm_ = zero
    !
    ! dm in the impurity problem basis
    dm_ = nn2so_reshape(imp_density_matrix,Nspin,Norb)
    !
  end subroutine ed_get_density_matrix_single








  !+------------------------------------------------------------------+
  !                         PRINT SIGMA:
  !+------------------------------------------------------------------+  
  subroutine ed_print_impSigma
    integer                                           :: i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                                 :: suffix
    integer,dimension(:),allocatable                  :: getIorb,getJorb
    integer                                           :: totNorb,l
    !
    call allocate_grids()
    !
    totNorb=Norb*(Norb+1)/2
    allocate(getIorb(totNorb),getJorb(totNorb))
    l=0
    do iorb=1,Norb
       do jorb=iorb,Norb
          l=l+1
          getIorb(l)=iorb
          getJorb(l)=jorb
       enddo
    enddo
    if(l/=totNorb)stop "ed_print_impSigma error counting the orbitals"
    !!
    !Print the impurity functions:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          call splot("impSigma"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impSmats(ispin,ispin,iorb,jorb,:))
          call splot("impSigma"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impSreal(ispin,ispin,iorb,jorb,:))
       enddo
    enddo
    !
    call deallocate_grids()
    !
  end subroutine ed_print_impSigma




  !+------------------------------------------------------------------+
  !                         PRINT G
  !+------------------------------------------------------------------+  
  subroutine ed_print_impG
    integer                                           :: i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                                 :: suffix
    integer,dimension(:),allocatable                  :: getIorb,getJorb
    integer                                           :: totNorb,l
    !
    call allocate_grids()
    !
    totNorb=Norb*(Norb+1)/2
    allocate(getIorb(totNorb),getJorb(totNorb))
    l=0
    do iorb=1,Norb
       do jorb=iorb,Norb
          l=l+1
          getIorb(l)=iorb
          getJorb(l)=jorb
       enddo
    enddo
    if(l/=totNorb)stop "ed_print_impG error counting the orbitals"
    !!
    !Print the impurity functions:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          call splot("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impGmats(ispin,ispin,iorb,jorb,:))
          call splot("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impGreal(ispin,ispin,iorb,jorb,:))
       enddo
    enddo
    !
    call deallocate_grids()
    !
  end subroutine ed_print_impG



  !+------------------------------------------------------------------+
  !                         PRINT G0
  !+------------------------------------------------------------------+  
  subroutine ed_print_impG0
    integer                                           :: i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                                 :: suffix
    integer,dimension(:),allocatable                  :: getIorb,getJorb
    integer                                           :: totNorb,l
    !
    call allocate_grids()
    !
    totNorb=Norb*(Norb+1)/2
    allocate(getIorb(totNorb),getJorb(totNorb))
    l=0
    do iorb=1,Norb
       do jorb=iorb,Norb
          l=l+1
          getIorb(l)=iorb
          getJorb(l)=jorb
       enddo
    enddo
    if(l/=totNorb)stop "ed_print_impG0 error counting the orbitals"
    !!
    !Print the impurity functions:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          call splot("impG0"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impG0mats(ispin,ispin,iorb,jorb,:))
          call splot("impG0"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impG0real(ispin,ispin,iorb,jorb,:))
       enddo
    enddo
    !
    call deallocate_grids()
    !
  end subroutine ed_print_impG0



  !+------------------------------------------------------------------+
  !                         PRINT D (phonon Green's function)
  !+------------------------------------------------------------------+  
  subroutine ed_print_impD
    !
    call allocate_grids()
    !
    !Print the impurity functions:
    call splot("impDph_iw.ed"   ,vm,impDmats_ph(:))
    call splot("impDph_realw.ed",vr,impDreal_ph(:))
    !
    call deallocate_grids()
    !
  end subroutine ed_print_impD



  !+------------------------------------------------------------------+
  !                         PRINT CHI:
  !+------------------------------------------------------------------+  
  subroutine ed_print_impChi
    if(chispin_flag)call print_chi_spin
    if(chidens_flag)call print_chi_dens
    if(chipair_flag)call print_chi_pair
    if(chiexct_flag)call print_chi_exct
  end subroutine ed_print_impChi

  !                         SPIN-SPIN
  subroutine print_chi_spin
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("spinChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,spinChi_tau(iorb,jorb,0:))
          call splot("spinChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,spinChi_w(iorb,jorb,:))
          call splot("spinChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,spinChi_iv(iorb,jorb,:))
       enddo
    enddo
    call deallocate_grids()
  end subroutine print_chi_spin
  !                     DENSITY-DENSITY
  subroutine print_chi_dens
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("densChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,densChi_tau(iorb,jorb,0:))
          call splot("densChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,densChi_w(iorb,jorb,:))
          call splot("densChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,densChi_iv(iorb,jorb,:))
       enddo
    enddo
    call deallocate_grids()
  end subroutine print_chi_dens
  !                     PAIR-PAIR
  subroutine print_chi_pair
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("pairChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,pairChi_tau(iorb,jorb,0:))
          call splot("pairChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,pairChi_w(iorb,jorb,:))
          call splot("pairChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,pairChi_iv(iorb,jorb,:))
       enddo
    enddo
    call deallocate_grids()
  end subroutine print_chi_pair
  !                     EXCITON
  subroutine print_chi_exct
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(0,iorb,jorb,0:))
          call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(0,iorb,jorb,:))
          call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(0,iorb,jorb,:))
       enddo
    enddo
    !
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(1,iorb,jorb,0:))
          call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(1,iorb,jorb,:))
          call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(1,iorb,jorb,:))
       enddo
    enddo
    !
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(2,iorb,jorb,0:))
          call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(2,iorb,jorb,:))
          call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(2,iorb,jorb,:))
       enddo
    enddo
    !
    call deallocate_grids()
  end subroutine print_chi_exct












END MODULE ED_IO







