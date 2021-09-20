MODULE ED_GF_ELECTRON
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN
  implicit none
  private


  public :: build_gf_normal
  public :: build_sigma_normal


  integer                               :: istate
  integer                               :: isector,jsector
  real(8),allocatable                   :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)
  integer                               :: ialfa,jalfa
  integer                               :: ipos,jpos
  integer                               :: i,j
  integer                               :: iph,i_el
  real(8)                               :: sgn,norm2
  real(8),dimension(:),pointer          :: state_cvec
  real(8)                               :: state_e


contains



  !+------------------------------------------------------------------+
  !                        NORMAL
  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the Green's function of the impurity electrons
  subroutine build_gf_normal()
    integer :: ispin,i
    integer :: iorb,jorb
    integer :: isite,jsite
    integer :: io,jo
    !

    do ispin=1,Nspin
       do iorb=1,Norb
          do isite=1,Nsites(iorb)
             write(LOGfile,"(A)")"Get G:"//&
                  " site I"//str(isite,site_indx_padding)//&
                  " orb M"//str(iorb)//&
                  " spin"//str(ispin)
             if(MPIMASTER)call start_timer
             call lanc_build_gf_normal_diag(isite,iorb,ispin)
             if(MPIMASTER)call stop_timer(unit=LOGfile)
          enddo
       enddo
    enddo
    !
    if(offdiag_gf_flag)then     
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                do isite=1,Nsites(iorb)
                   do jsite=1,Nsites(jorb)
                      if(io==jo)cycle
                      write(LOGfile,"(A)")"Get G:"//&
                           " sites I"//str(isite,site_indx_padding)//&
                           "J"//str(jsite,site_indx_padding)//&
                           " orb M"//str(iorb)//"L"//str(jorb)//&
                           " spin"//str(ispin)
                      if(MPIMASTER)call start_timer
                      call lanc_build_gf_normal_mix(isite,jsite,iorb,jorb,ispin)
                      if(MPIMASTER)call stop_timer(unit=LOGfile)
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
       !
       !Put here off-diagonal manipulation by symmetry:
       do ispin=1,Nspin
          do io=1,Ns
             do jo=1,Ns
                if(io==jo)cycle
                impGmats(ispin,io,jo,:) = 0.5d0*(impGmats(ispin,io,jo,:) - &
                     impGmats(ispin,io,io,:) - impGmats(ispin,jo,jo,:))
                impGreal(ispin,io,jo,:) = 0.5d0*(impGreal(ispin,io,jo,:) - &
                     impGreal(ispin,io,io,:) - impGreal(ispin,jo,jo,:))
             enddo
          enddo
       enddo
    end if
    !
  end subroutine build_gf_normal






  !################################################################
  !################################################################
  !################################################################
  !################################################################





  subroutine lanc_build_gf_normal_diag(isite,iorb,ispin)
    integer,intent(in) :: isite,iorb,ispin
    integer            :: io
    type(sector)       :: sectorI,sectorJ
    !
    ialfa = 1
    ipos  = pack_indices(isite,iorb)
    !
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate) 
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c^+_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,io,io,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       endif
       !
       !REMOVE ONE PARTICLE:
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,io,io,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_diag




  !################################################################






  subroutine lanc_build_gf_normal_mix(isite,jsite,iorb,jorb,ispin)
    integer      :: isite,jsite,iorb,jorb,ispin
    integer      :: io,jo
    type(sector) :: sectorI,sectorJ
    !
    ialfa = 1
    jalfa = ialfa
    ipos  = pack_indices(isite,iorb)
    jpos  = pack_indices(jsite,jorb)
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !EVALUATE (c^+_io + c^+_jo)|gs>
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c^+_a,s + c^+_b,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             !c^+_io|gs>
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             !+c^+_jo|gs>
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,jpos,jalfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,io,jo,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE (c_io + c_jo)|gs>
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  '  apply c_a,s + c_b,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             !c_io|gs>
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             !+c_jo|gs>
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,jpos,jalfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,io,jo,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix





  !################################################################





  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,io,jo,ispin)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,io,jo,ispin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if((finiteT).and.(beta*(Ei-Egs).lt.200))then
       pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    elseif(.not.finiteT)then
       pesoBZ = vnorm2/zeta_function
    else
       pesoBZ=0.d0
    endif
    !
    !pesoBZ = vnorm2/zeta_function
    !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
    !Only the nodes in Mpi_Comm_Group did get the alanc,blanc.
    !However after delete_sectorHv MpiComm returns to be the global one
    !so we can safely Bcast the alanc,blanc (known only to the operative group)
    !to every nodes. The master is in charge of this (as a
    !participant of the operative group)
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    !
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       do i=1,Lmats
          iw=xi*wm(i)
          impGmats(ispin,io,jo,i)=impGmats(ispin,io,jo,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          impGreal(ispin,io,jo,i)=impGreal(ispin,io,jo,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_normal










  subroutine build_sigma_normal
    integer                                 :: i,ispin,iorb
    complex(8),dimension(Nspin,Ns,Ns,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nspin,Ns,Ns,Lreal) :: invG0real,invGreal
    complex(8),dimension(Ns,Ns)             :: Gtmp
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    call Hij_get_g0inv(dcmplx(0d0,wm(:)),invG0mats)
    call Hij_get_g0inv(dcmplx(wr(:),eps),invG0real)
    !

    impSmats=zero
    impSreal=zero
    !Get Gimp^-1
    do ispin=1,Nspin
       do i=1,Lmats
          Gtmp = impGmats(ispin,:,:,i)
          call inv(Gtmp)
          invGmats(ispin,:,:,i)=Gtmp
       enddo
       !
       do i=1,Lreal
          Gtmp = impGreal(ispin,:,:,i)
          call inv(Gtmp)
          invGreal(ispin,:,:,i)=Gtmp
       enddo
       !
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats(ispin,:,:,:) = invG0mats(ispin,:,:,:) - invGmats(ispin,:,:,:)
       impSreal(ispin,:,:,:) = invG0real(ispin,:,:,:) - invGreal(ispin,:,:,:)
    enddo
    !
    !Get G0and:
    call Hij_get_g0func(dcmplx(0d0,wm(:)),impG0mats)
    call Hij_get_g0func(dcmplx(wr(:),eps),impG0real)
    !
  end subroutine build_sigma_normal


END MODULE ED_GF_ELECTRON











