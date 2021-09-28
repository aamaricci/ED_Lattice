MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_LINALG
  USE SF_MISC,    only: assert_shape
  USE SF_ARRAYS,  only: linspace,arange
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
  public :: ed_get_elocal
  public :: ed_get_doubles
  public :: ed_get_density_matrix
  !
  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impChi

  character(len=128) :: suffix
  character(len=128) :: gf_suffix='.ed'
  character(len=8)   :: w_suffix






contains



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity GF and self-energy 
  !+--------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA GREEN'S FUNCTIONS
  subroutine ed_get_gimp_matsubara(Func) 
    complex(8),dimension(Nspin,Ns,Ns,Lmats),intent(inout) :: Func
    Func = impGmats
  end subroutine ed_get_gimp_matsubara

  !NORMAL, REAL SELF-ENERGY
  subroutine ed_get_gimp_realaxis(Func) 
    complex(8),dimension(Nspin,Ns,Ns,Lmats),intent(inout) :: Func
    Func = impGreal
  end subroutine ed_get_gimp_realaxis

  !NORMAL, MATSUBARA SELF-ENEGRGY
  subroutine ed_get_sigma_matsubara(Func) 
    complex(8),dimension(Nspin,Ns,Ns,Lmats),intent(inout) :: Func
    Func = impSmats
  end subroutine ed_get_sigma_matsubara

  !NORMAL, REAL SELF-ENERGY
  subroutine ed_get_sigma_realaxis(Func) 
    complex(8),dimension(Nspin,Ns,Ns,Lmats),intent(inout) :: Func
    Func = impSreal
  end subroutine ed_get_sigma_realaxis



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+--------------------------------------------------------------------------+!
  subroutine ed_get_dens(dens)
    real(8),dimension(Ns) :: dens
    dens = ed_dens
  end subroutine ed_get_dens

  subroutine ed_get_mag(mag) 
    real(8),dimension(Ns) :: mag
    mag = (ed_dens_up-ed_dens_dw)
  end subroutine ed_get_mag

  subroutine ed_get_docc(docc) 
    real(8),dimension(Ns) :: docc
    docc = ed_docc
  end subroutine ed_get_docc

  subroutine ed_get_doubles(docc)
    real(8),dimension(4) :: docc
    docc = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
  end subroutine ed_get_doubles

  subroutine ed_get_elocal(eimp)
    real(8),dimension(4) :: eimp
    eimp = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
  end subroutine ed_get_elocal

  subroutine ed_get_density_matrix(rho) 
    complex(8),dimension(Nspin,Ns,Ns),intent(inout) :: rho
    rho = imp_density_matrix
  end subroutine ed_get_density_matrix









  !+------------------------------------------------------------------+
  !                         PRINT LATTICE FUNCTIONS:
  !+------------------------------------------------------------------+
  subroutine ed_print_impSigma
    call allocate_grids()
    call ed_write_func(impSmats,"Sigma",'mats',wm,3)
    call ed_write_func(impSreal,"Sigma",'real',wr,3)
    call deallocate_grids()
  end subroutine ed_print_impSigma

  subroutine ed_print_impG
    call allocate_grids()
    call ed_write_func(impGmats,"G",'mats',wm,3)
    call ed_write_func(impGreal,"G",'real',wr,3)
    call deallocate_grids()
  end subroutine ed_print_impG

  subroutine ed_print_impG0
    call allocate_grids()
    call ed_write_func(impG0mats,"G0",'mats',wm,3)
    call ed_write_func(impG0real,"G0",'real',wr,3)
    call deallocate_grids()
  end subroutine ed_print_impG0


  !+------------------------------------------------------------------+
  !                         PRINT CHI:
  !+------------------------------------------------------------------+  
  subroutine ed_print_impChi
    if(chispin_flag)call print_chi_spin
  end subroutine ed_print_impChi

  !                         SPIN-SPIN
  subroutine print_chi_spin
    integer :: io,jo
    call allocate_grids()
    do io=1,Ns
       do jo=1,Ns
          call splot("spinChi_l"//str(io)//str(io)//"_tau.ed",tau,spinChi_tau(io,jo,0:))
          call splot("spinChi_l"//str(io)//str(jo)//"_realw.ed",vr,spinChi_w(io,jo,:))
          call splot("spinChi_l"//str(io)//str(io)//"_iv.ed",vm,spinChi_iv(io,jo,:))
       enddo
    enddo
    call deallocate_grids()
  end subroutine print_chi_spin








  !+------------------------------------------------------------------+  
  subroutine ed_write_func(Func,fname,axis,zeta,iprint)
    complex(8),dimension(:,:,:,:),intent(in) :: Func ![Nspin][Ns][Ns][Lfreq]
    character(len=*),intent(in)              :: fname
    character(len=*)                         :: axis
    real(8),dimension(:)                     :: zeta
    integer,optional                         :: iprint
    !
    integer                                  :: iprint_
    integer                                  :: Lfreq,Nfunc   
    integer                                  :: io,jo,ilat,jlat,iorb,jorb,ispin,jspin
    !
    iprint_=1;if(present(iprint))iprint_=iprint
    !
    select case(axis)
    case default;stop "ed_write_func ERROR: axis undefined. axis=[matsubara,realaxis]"
    case("matsubara","mats");w_suffix="_iw"
    case("realaxis","real") ;w_suffix="_realw"
    end select
    !
    !
    Lfreq = size(zeta)
    Nfunc = size(Func,2)
    if(Nfunc /= Ns)stop "Ed_write_func ERROR: size(func,1) /= sum(Nsites)==Ns"
    call assert_shape(Func,[Nspin,Ns,Ns,Lfreq],"Ed_write_func","Func")
    !
    !
    !1: diagonal lattice-spin-orbital
    !2: all sites, diagonal spin-orbital
    !3: all sites, all orb, same spins
    !4: all sites, all orb, all spins
    if(MpiMaster)then
       select case(iprint_)
       case default
          write(*,"(A,1x,A)")reg(fname),"ed_write_func: not written on file."
          !
          !
       case(1)
          write(*,"(A,1x,A)") reg(fname),"ed_write_func: diagonal lattice-spin-orbital."
          do ispin=1,Nspin
             do iorb=1,Norb
                do ilat=1,Nsites(iorb)
                   io = pack_indices(ilat,iorb)
                   !
                   suffix=reg(fname)//&
                        "_i"//str(ilat,site_indx_padding)//&
                        "_l"//str(iorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//str(gf_suffix)
                   call splot(reg(suffix),zeta,Func(ispin,io,io,:))
                enddo
             enddo
          enddo
          !
          !
       case(2)
          write(*,"(A,1x,A)") reg(fname),"ed_write_func:  all sites, diagonal spin-orbital."
          do ispin=1,Nspin
             do iorb=1,Norb
                do ilat=1,Nsites(iorb)
                   do jlat=1,Nsites(iorb)
                      io = pack_indices(ilat,iorb)
                      jo = pack_indices(jlat,iorb)
                      !
                      suffix=reg(fname)//&
                           "_i"//str(ilat,site_indx_padding)//"j"//str(jlat,site_indx_padding)//&
                           "_l"//str(iorb)//&
                           "_s"//str(ispin)//&
                           str(w_suffix)//str(gf_suffix)
                      call splot(reg(suffix),zeta,Func(ispin,io,jo,:))
                   enddo
                enddo
             enddo
          enddo
          !
          !
       case(3)
          write(*,"(A,1x,A)") reg(fname),"ed_write_func: all sites, all orb, same spins."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   do ilat=1,Nsites(iorb)
                      do jlat=1,Nsites(jorb)
                         io = pack_indices(ilat,iorb)
                         jo = pack_indices(jlat,jorb)
                         !
                         suffix=reg(fname)//&
                              "_i"//str(ilat,site_indx_padding)//"j"//str(jlat,site_indx_padding)//&
                              "_l"//str(iorb)//"m"//str(jorb)//&
                              "_s"//str(ispin)//&
                              str(w_suffix)//reg(gf_suffix)
                         call splot(reg(suffix),zeta,Func(ispin,io,jo,:))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
          !
       end select
    endif
  end subroutine ed_write_func


END MODULE ED_IO







