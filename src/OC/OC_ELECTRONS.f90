MODULE ED_OC_ELECTRONS
  USE SF_CONSTANTS, only:one,xi,zero,pi,isnan,isinfty
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  implicit none
  private


  public :: build_oc_electrons
  public :: eval_oc_electrons

  integer                                 :: istate,iorb,jorb,ispin
  integer                                 :: isector
  complex(8),allocatable                  :: vvinit(:)
  real(8),allocatable                     :: alfa_(:),beta_(:)
  integer                                 :: ialfa
  integer                                 :: jalfa
  integer                                 :: ipos,jpos
  integer                                 :: i,j
  real(8)                                 :: norm2
  complex(8),dimension(:),pointer         :: state_cvec
  real(8)                                 :: state_e
  complex(8),dimension(:,:,:),allocatable :: Hij



contains


  subroutine build_oc_electrons()
    integer :: iorb
    !
    call deallocate_GFmatrix(OcMatrix)
    !
    select case(ed_method)
    case ('lapack','full')
       !do nothing
       return
    case default
       allocate(Hij(Nspin,Ns,Ns))
       call Hij_get(Hij)
       do iorb=1,Norb
          if(.not.oc_flag(iorb))cycle
          if(MPIMASTER)write(LOGfile,"(A)")"Build OC:"//" orb M"//str(iorb)
          if(MPIMASTER)call start_timer
          call allocate_GFmatrix(OcMatrix(iorb),Nstate=state_list%size)
          call lanc_ed_build_oc(iorb)
          if(MPIMASTER)call stop_timer(unit=LOGfile)
       enddo
       deallocate(Hij)
    end select
    !
  end subroutine build_oc_electrons


  subroutine eval_oc_electrons()
    integer :: iorb
    !
    do iorb=1,Norb
       if(.not.oc_flag(iorb))cycle
       if(MPIMASTER)write(LOGfile,"(A)")"Eval OC:"//" orb M"//str(iorb)
       if(MPIMASTER)call start_timer
       select case(ed_method)
       case default
          call lanc_ed_eval_oc(iorb)
       case ('lapack','full')
          call full_ed_eval_oc(iorb)
       end select
       if(MPIMASTER)call stop_timer(unit=LOGfile)
    enddo
    !
  end subroutine eval_oc_electrons





  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine lanc_ed_build_oc(iorb)
    integer,intent(in) :: iorb
    integer            :: iup,idw,jup,jdw
    integer            :: mup,mdw
    integer            :: io,jo
    integer            :: k1,k2
    real(8)            :: sg1,sg2
    integer            :: nup(Ns),ndw(Ns)
    integer            :: isite,jsite
    type(sector)       :: sectorI
    !
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(OcMatrix(iorb),istate,Nchan=1)
       !
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
               'Apply J  :',isector,sectorI%Nups,sectorI%Ndws
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          !
          do i=1,sectorI%Dim
             iup = iup_index(i,sectorI%DimUp)
             idw = idw_index(i,sectorI%DimUp)
             !
             mup = sectorI%H(1)%map(iup)
             mdw = sectorI%H(2)%map(idw)
             !
             nup = bdecomp(mup,Ns)
             ndw = bdecomp(mdw,Ns)
             !
             !Apply J_{iorb} =  -xi*\sum_sigma \sum_i
             !                  H(sigma,iorb,i,i+1) Cdg_{sigma,iorb,i}C_{sigma,iorb,i+1} -
             !                  H*(sigma,iorb,i+1,i) Cdg_{sigma,iorb,i+1}C_{sigma,iorb,i}
             !     == Jdg_{iorb}
             do isite=1,Nsites(iorb)
                !H(sigma,iorb,i,i+1) Cdg_{sigma,iorb,i}C_{sigma,iorb,i+1} ==
                jsite = isite+1
                if(isite==Nsites(iorb))jsite=1
                io = pack_indices(isite,iorb)
                jo = pack_indices(jsite,iorb)
                ! UP
                if((nup(jo)==1) .AND. (nup(io)==0) .AND. (Hij(1,io,jo)/=zero) )then
                   call c(jo,mup,k1,sg1)
                   call cdg(io,k1,k2,sg2)
                   jup= binary_search(sectorI%H(1)%map,k2)
                   j  = jup + (idw-1)*sectorI%DimUp                   
                   vvinit(j) = vvinit(j) + Hij(1,io,jo)*sg1*sg2*state_cvec(i)
                endif
                ! DW
                if((ndw(jo)==1) .AND. (ndw(io)==0).AND. (Hij(Nspin,io,jo)/=zero))then
                   call c(jo,mdw,k1,sg1)
                   call cdg(io,k1,k2,sg2)
                   jdw= binary_search(sectorI%H(2)%map,k2)
                   j  = iup + (jdw-1)*sectorI%DimUp                   
                   vvinit(j) = vvinit(j) + Hij(Nspin,io,jo)*sg1*sg2*state_cvec(i)
                endif
                !
                !
                !-H*(sigma,iorb,i+1,i) Cdg_{sigma,iorb,i+1}C_{sigma,iorb,i} == 
                jsite = isite+1
                if(isite==Nsites(iorb))jsite=1
                io = pack_indices(jsite,iorb)
                jo = pack_indices(isite,iorb)
                ! UP
                if((nup(jo)==1) .AND. (nup(io)==0) .AND. (Hij(1,io,jo)/=zero) )then
                   call c(jo,mup,k1,sg1)
                   call cdg(io,k1,k2,sg2)
                   jup= binary_search(sectorI%H(1)%map,k2)
                   j  = jup + (idw-1)*sectorI%DimUp                   
                   vvinit(j) = vvinit(j) - conjg(Hij(1,io,jo))*sg1*sg2*state_cvec(i)
                endif
                ! DW
                if((ndw(jo)==1) .AND. (ndw(io)==0).AND. (Hij(Nspin,io,jo)/=zero))then
                   call c(jo,mdw,k1,sg1)
                   call cdg(io,k1,k2,sg2)
                   jdw= binary_search(sectorI%H(2)%map,k2)
                   j  = iup + (jdw-1)*sectorI%DimUp                   
                   vvinit(j) = vvinit(j) - conjg(Hij(Nspin,io,jo))*sg1*sg2*state_cvec(i)
                endif
             enddo
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=zero
       endif
       vvinit = -xi*vvinit
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_oc(norm2,state_e,alfa_,beta_,iorb,ichan=1,istate=istate)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
    enddo
    return
  end subroutine lanc_ed_build_oc




  !################################################################



  subroutine add_to_lanczos_oc(vnorm2,Ei,alanc,blanc,io,ichan,istate)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: io,ichan,istate
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2 
    ! pesoBZ = 1d0/zeta_function
    ! if(finiteT)pesoBZ = exp(-(Ei-Egs)/temp)/zeta_function
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(OcMatrix(io),istate,ichan,Nlanc)
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB!*pesoBZ
       !
       OcMatrix(io)%state(istate)%channel(ichan)%weight(j) = peso
       OcMatrix(io)%state(istate)%channel(ichan)%poles(j)  = de

    enddo
  end subroutine add_to_lanczos_oc




  !################################################################
  !################################################################
  !################################################################
  !################################################################



  subroutine lanc_ed_eval_oc(iorb)
    integer,intent(in) :: iorb
    integer            :: Nstates,istate
    integer            :: Nchannels,ichan
    integer            :: Nexcs,iexc
    real(8)            :: peso,de,pesoBZ,beta,Ei,Egs
    !
    !    
    if(.not.allocated(OcMatrix(iorb)%state)) then
       print*, "GF_NORMAL WARNING: OcMatrix%state not allocated. Nothing to do"
       return
    endif
    !
    beta= 1d0/temp
    Egs = state_list%emin
    pesoBZ = 1d0/zeta_function
    !
    !this is the total number of available states  == state_list%size
    Nstates = size(OcMatrix(iorb)%state) 
    !Get trimmed state for the actual value of temp == state_list%trimd_size
    call es_trim_size(state_list,temp,cutoff) 
    do istate=1,state_list%trimd_size     !Nstates
       if(.not.allocated(OcMatrix(iorb)%state(istate)%channel))cycle
       Ei =  es_return_energy(state_list,istate)
       if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))/zeta_function
       Nchannels = size(OcMatrix(iorb)%state(istate)%channel)
       do ichan=1,Nchannels
          Nexcs  = size(OcMatrix(iorb)%state(istate)%channel(ichan)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = OcMatrix(iorb)%state(istate)%channel(ichan)%weight(iexc)
             peso  = peso*pesoBZ
             dE    = OcMatrix(iorb)%state(istate)%channel(ichan)%poles(iexc)
             !
             if( (isnan(peso/de)).OR.(isinfty(peso/de)) )cycle
             if(abs(dE)<1d-12)cycle
             Drude_weight(iorb) = Drude_weight(iorb) + peso/dE
             do i=1,Lreal
                OptCond_w(iorb,i)=OptCond_w(iorb,i) + peso/dE*eps/( (vr(i)-dE)**2 + eps**2 )
             enddo
          enddo
       enddo
    enddo
    return
  end subroutine lanc_ed_eval_oc



  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine full_ed_eval_oc(iorb)
    integer      :: isite,jsite,iorb,jorb
    integer      :: io,jo,mup,mdw,iup,idw,jup,jdw
    type(sector) :: sectorI,sectorJ
    real(8)      :: sg1,sg2
    complex(8)   :: Chio
    integer      :: Nups(Ns_Ud)
    integer      :: Ndws(Ns_Ud)
    integer      :: i,j,ll,m,isector,k1,k2,li,rj
    integer      :: idim,ia,nup(Ns),ndw(Ns)
    real(8)      :: Ei,Ej,cc,peso,pesotot
    real(8)      :: expterm,de,w0,it
    complex(8)   :: iw 
    !
    !
    !
    do isector=1,Nsectors !loop over <i| total particle number
       call get_Nup(isector,nups)
       call get_Ndw(isector,ndws)
       if(ed_filling/=0 .AND. (sum(Nups)+sum(Ndws)/=ed_filling) )cycle
       !
       call eta(isector,Nsectors,LOGfile)
       call build_sector(isector,sectorI)
       !
       do i=1,sectorI%Dim 
          do j=1,sectorI%Dim
             Chio=zero
             expterm=exp(-espace(isector)%e(i)/temp)+exp(-espace(isector)%e(j)/temp)
             if(expterm<cutoff)cycle
             do li=1,sectorI%Dim
                iup = iup_index(li,sectorI%DimUp)
                idw = idw_index(li,sectorI%DimUp)
                !
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                !
                nup = bdecomp(mup,Ns)
                ndw = bdecomp(mdw,Ns)
                !
                !Apply J_{iorb} =  xi*\sum_sigma \sum_i
                !                  H(sigma,iorb,i,i+1) Cdg_{sigma,iorb,i}C_{sigma,iorb,i+1} -
                !                  H*(sigma,iorb,i+1,i) Cdg_{sigma,iorb,i+1}C_{sigma,iorb,i}
                !     == J^+_{iorb}
                do isite=1,Nsites(iorb)-1
                   ! i->i+1
                   io = pack_indices(isite,iorb)
                   jo = pack_indices(isite+1,iorb)
                   !H(sigma,iorb,i,i+1) Cdg_{sigma,iorb,i}C_{sigma,iorb,i+1} ==
                   !H(sigma,iorb,io,jo) Cdg_{sigma,iorb,io}C_{sigma,iorb,jo}
                   ! UP
                   if( (nup(jo)==1) .AND. (nup(io)==0) .AND. (Hij(1,io,jo)/=zero) )then
                      call c(jo,mup,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jup =binary_search(sectorI%H(1)%map,k2)
                      rj  = jup + (idw-1)*sectorI%DimUp                   
                      Chio= Chio + conjg(espace(isector)%M(rj,j))*Hij(1,io,jo)*sg1*sg2*espace(isector)%M(li,i)
                   endif
                   ! DW
                   if((ndw(jo)==1) .AND. (ndw(io)==0).AND. (Hij(Nspin,io,jo)/=zero))then
                      call c(jo,mdw,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jdw = binary_search(sectorI%H(2)%map,k2)
                      rj  = iup + (jdw-1)*sectorI%DimUp                   
                      Chio= Chio + conjg(espace(isector)%M(rj,j))*Hij(Nspin,io,jo)*sg1*sg2*espace(isector)%M(li,i)
                   endif
                   !
                   !
                   ! i+1->i
                   io = pack_indices(isite+1,iorb)
                   jo = pack_indices(isite,iorb)
                   !-H*(sigma,iorb,i+1,i) Cdg_{sigma,iorb,i+1}C_{sigma,iorb,i} == 
                   !-H*(sigma,iorb,io,jo) Cdg_{sigma,iorb,io}C_{sigma,iorb,jo}
                   ! UP
                   if((nup(jo)==1) .AND. (nup(io)==0) .AND. (Hij(1,io,jo)/=zero) )then
                      call c(jo,mup,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jup =binary_search(sectorI%H(1)%map,k2)
                      rj  = jup + (idw-1)*sectorI%DimUp                   
                      Chio= Chio - conjg(espace(isector)%M(rj,j))*Hij(1,io,jo)*sg1*sg2*espace(isector)%M(li,i)
                   endif
                   ! DW
                   if((ndw(jo)==1) .AND. (ndw(io)==0).AND. (hij(Nspin,io,jo)/=zero))then
                      call c(jo,mdw,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jdw = binary_search(sectorI%H(2)%map,k2)
                      rj  = iup + (jdw-1)*sectorI%DimUp                   
                      Chio= Chio - conjg(espace(isector)%M(rj,j))*hij(Nspin,io,jo)*sg1*sg2*espace(isector)%M(li,i)
                   endif
                enddo
             enddo
             Chio = xi*Chio
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             dE=Ei-Ej
             peso = abs(Chio)**2/zeta_function
             !
             if( (isnan(peso/de)).OR.(isinfty(peso/de)) )cycle
             if(abs(dE)<1d-12)cycle
             Drude_weight(iorb) = Drude_weight(iorb) + peso/de
             !Real-frequency: Retarded = Commutator = response function
             do m=1,Lreal
                iw=vr(m)-dE
                OptCond_w(iorb,m)=OptCond_w(iorb,m)+peso/dE*eps/(iw**2+eps**2)
             enddo
             !
          enddo
       enddo
       call delete_sector(sectorI)
    enddo

  end subroutine full_ed_eval_oc





END MODULE ED_OC_ELECTRONS
























