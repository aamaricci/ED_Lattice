MODULE ED_CHI_SPIN_ELECTRONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
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


  public :: build_chi_spin_electrons
  public :: eval_chi_spin_electrons



  integer                         :: istate,iorb,jorb,ispin
  integer                         :: isector
  complex(8),allocatable          :: vvinit(:)
  real(8),allocatable             :: alfa_(:),beta_(:)
  integer                         :: ipos,jpos
  integer                         :: i,j
  real(8)                         :: sgn,norm2
  real(8)                         :: state_e




contains


  !+------------------------------------------------------------------+
  !                            SPIN
  !PURPOSE  : Evaluate the Spin susceptibility \Chi_spin for a 
  ! single orbital: \chi = <S_a(\tau)S_a(0)>
  ! note: as S_a is hermitian particle and holes contributions
  ! are identical so work out only one lanczos tridiag. work out the 
  ! reduction for both values of isign in the same call.
  !+------------------------------------------------------------------+
  subroutine build_chi_spin_electrons()
    integer :: ispin,i,iimp
    integer :: iorb,jorb
    integer :: isite,jsite
    integer :: io,jo
    !
    select case(ed_method)
    case ('lapack','full')
       return
    case default
       !Diagonal
       do iorb=1,Norb
          if(.not.chispin_flag(iorb))cycle
          do isite=1,Nsites(iorb)
             io  = pack_indices(isite,iorb)
             if(MPIMASTER)write(LOGfile,"(A)")"Build spinChi:"//&
                  " site I"//str(isite,site_indx_padding)//&
                  " orb M"//str(iorb)
             if(MPIMASTER)call start_timer
             call allocate_GFmatrix(SpinChiMatrix(io,io),Nstate=state_list%size)
             call lanc_build_spinChi_diag(isite,iorb)
             if(MPIMASTER)call stop_timer(unit=LOGfile)
          enddo
       enddo
       !Off-diagonal
       if(offdiag_chispin_flag.AND.Norb>1)then
          do iorb=1,Norb
             do jorb=1,Norb
                if(.not.chispin_flag(iorb).AND..not.chispin_flag(jorb))cycle
                do isite=1,Nsites(iorb)
                   do jsite=1,Nsites(jorb)
                      io  = pack_indices(isite,iorb)
                      jo  = pack_indices(jsite,jorb)
                      if(io==jo)cycle
                      if(MPIMASTER)write(LOGfile,"(A)")"Build spinChi:"//&
                           " sites I"//str(isite,site_indx_padding)//&
                           "J"//str(jsite,site_indx_padding)//&
                           " orb M"//str(iorb)//"L"//str(jorb)
                      if(MPIMASTER)call start_timer
                      call allocate_GFmatrix(SpinChiMatrix(io,jo),Nstate=state_list%size)
                      call lanc_build_spinChi_mix(isite,jsite,iorb,jorb)
                      if(MPIMASTER)call stop_timer(unit=LOGfile)
                   enddo
                enddo
             enddo
          enddo
       endif
       !
    end select
  end subroutine build_chi_spin_electrons

  subroutine eval_chi_spin_electrons()
    integer :: ispin,i,iimp
    integer :: iorb,jorb
    integer :: isite,jsite
    integer :: io,jo
    !    
    do iorb=1,Norb
       if(.not.chispin_flag(iorb))cycle
       do isite=1,Nsites(iorb)
          if(MPIMASTER)write(LOGfile,"(A)")"Eval spinChi:"//&
               " site I"//str(isite,site_indx_padding)//&
               " orb M"//str(iorb)
          io  = pack_indices(isite,iorb)
          if(MPIMASTER)call start_timer
          select case(ed_method)
          case default
             call lanc_eval_spinChi_electrons(isite,isite,iorb,iorb)
          case ('lapack','full')
             call full_eval_spinChi_electrons(isite,isite,iorb,iorb)
          end select
          if(MPIMASTER)call stop_timer(unit=LOGfile)
       enddo
    enddo
    !
    if(offdiag_chispin_flag.AND.Norb>1)then
       do iorb=1,Norb
          do jorb=1,Norb
             if(.not.chispin_flag(iorb).AND..not.chispin_flag(jorb))cycle
             do isite=1,Nsites(iorb)
                do jsite=1,Nsites(jorb)
                   io  = pack_indices(isite,iorb)
                   jo  = pack_indices(jsite,jorb)
                   if(io==jo)cycle
                   if(MPIMASTER)write(LOGfile,"(A)")"Eval spinChi:"//&
                        " sites I"//str(isite,site_indx_padding)//&
                        "J"//str(jsite,site_indx_padding)//&
                        " orb M"//str(iorb)//"L"//str(jorb)
                   if(MPIMASTER)call start_timer
                   select case(ed_method)
                   case default
                      call lanc_eval_spinChi_electrons(isite,jsite,iorb,jorb)
                   case ('lapack','full')
                      call full_eval_spinChi_electrons(isite,isite,iorb,iorb)
                   end select
                   if(MPIMASTER)call stop_timer(unit=LOGfile)
                enddo
             enddo
          enddo
       enddo
       !
       !
       select case(ed_method)
       case default
          do iorb=1,Norb
             do jorb=1,Norb
                if(.not.chispin_flag(iorb).AND..not.chispin_flag(jorb))cycle
                do isite=1,Nsites(iorb)
                   do jsite=1,Nsites(jorb)
                      io  = pack_indices(isite,iorb)
                      jo  = pack_indices(jsite,jorb)
                      if(io==jo)cycle
                      spinChi_w(io,jo,:)   = 0.5d0*(spinChi_w(io,jo,:) - spinChi_w(io,io,:) - spinChi_w(jo,jo,:))
                      spinChi_tau(io,jo,:) = 0.5d0*(spinChi_tau(io,jo,:) - spinChi_tau(io,io,:) - spinChi_tau(jo,jo,:))
                      spinChi_iv(io,jo,:)  = 0.5d0*(spinChi_iv(io,jo,:) - spinChi_iv(io,io,:) - spinChi_iv(jo,jo,:))
                      !
                      spinChi_w(jo,io,:)   = spinChi_w(io,jo,:)
                      spinChi_tau(jo,io,:) = spinChi_tau(io,jo,:)
                      spinChi_iv(jo,io,:)  = spinChi_iv(io,jo,:)
                   enddo
                enddo
             enddo
          enddo
       case ('lapack','full')
          continue
       end select
    endif
    !
  end subroutine eval_chi_spin_electrons




  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_build_spinChi_diag(isite,iorb)
    integer,intent(in)                  :: isite,iorb
    integer                             :: io
    type(sector)                        :: sectorI,sectorJ
    complex(8),dimension(:),allocatable :: state_cvec
    !
    io    = pack_indices(isite,iorb)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(SpinChiMatrix(io,io),istate,Nchan=1)
       !
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec) 
       else
          call es_return_cvector(state_list,istate,state_cvec) 
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")'Apply Sz_iorb  :',isector,sectorI%Nups,sectorI%Ndws
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          do i=1,sectorI%Dim
             call apply_op_Sz(i,sgn,io,sectorI)            
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=zero
       endif
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,io,io,ichan=1,istate=istate)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(state_cvec))deallocate(state_cvec)       
    enddo
    return
  end subroutine lanc_build_spinChi_diag




  !################################################################





  subroutine lanc_build_spinChi_mix(isite,jsite,iorb,jorb)
    integer                             :: isite,jsite,iorb,jorb
    integer                             :: io,jo
    type(sector)                        :: sectorI,sectorJ
    real(8)                             :: Siorb,Sjorb
    complex(8),dimension(:),allocatable :: state_cvec
    !
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(SpinChiMatrix(io,jo),istate,Nchan=1)
       !
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec) 
       else
          call es_return_cvector(state_list,istate,state_cvec) 
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
#endif
       !
       !EVALUATE (Sz_jorb + Sz_iorb)|gs> = Sz_jorb|gs> + Sz_iorb|gs>
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,2I4)")'Apply (Sz_jorb + Sz_iorb):',isector,sectorI%Nups,sectorI%Ndws
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          do i=1,sectorI%Dim
             call apply_op_Sz(i,Siorb,io,sectorI)
             call apply_op_Sz(i,Sjorb,jo,sectorI)
             sgn       = Siorb + Sjorb
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=zero
       endif
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,io,jo,ichan=1,istate=istate)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(state_cvec))deallocate(state_cvec)
    enddo
    return
  end subroutine lanc_build_spinChi_mix





  !################################################################





  subroutine add_to_lanczos_spinChi(vnorm2,Ei,alanc,blanc,io,jo,ichan,istate)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso,beta
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: io,jo,ichan,istate
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
    call allocate_GFmatrix(SpinChiMatrix(io,jo),istate,ichan,Nlanc)
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB
       !
       SpinChiMatrix(io,jo)%state(istate)%channel(ichan)%weight(j) = peso
       SpinChiMatrix(io,jo)%state(istate)%channel(ichan)%poles(j)  = de
       !
       ! ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! ! This ensures that the correct null contribution is obtained.
       ! ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! ! we do not include the contribution (because we are in the situation described above).
       ! ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       ! if(beta*dE > 1d-3)spinChi_iv(io,jo,0)=spinChi_iv(io,jo,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       ! do i=1,Lmats
       !    spinChi_iv(io,jo,i)=spinChi_iv(io,jo,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       ! enddo
       ! !Symmetrize for low-T /large-beta, mostly occurring for zero T calculations
       ! do i=0,Ltau/2-1
       !    spinChi_tau(io,jo,i)=spinChi_tau(io,jo,i) + peso*exp(-tau(i)*dE)
       ! enddo
       ! spinChi_tau(io,jo,Ltau/2)=spinChi_tau(io,jo,Ltau/2) + peso*0.5d0*(exp(-tau(Ltau/2)*dE)+exp(-(beta-tau(Ltau/2))*dE))
       ! do i=Ltau/2+1,Ltau
       !    spinChi_tau(io,jo,i)=spinChi_tau(io,jo,i) + peso*exp(-(beta-tau(i))*dE)
       ! enddo
       ! !
       ! do i=1,Lreal
       !    spinChi_w(io,jo,i)=spinChi_w(io,jo,i) - peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       ! enddo
    enddo
  end subroutine add_to_lanczos_spinChi








  !################################################################
  !################################################################
  !################################################################
  !################################################################



  subroutine lanc_eval_spinChi_electrons(isite,jsite,iorb,jorb)
    integer,intent(in)                  :: isite,jsite,iorb,jorb
    integer                             :: Nstates,istate
    integer                             :: Nchannels,ichan
    integer                             :: Nexcs,iexc,io,jo
    real(8)                             :: peso,de,pesoBZ,beta,Ei,Egs
    real(8),dimension(Ns,Ns,0:Ltau)     :: spinChi_tau_tmp
    complex(8),dimension(Ns,Ns,Lreal)   :: spinChi_w_tmp
    complex(8),dimension(Ns,Ns,0:Lmats) :: spinChi_iv_tmp
    !
    spinChi_tau_tmp=0d0
    spinChi_w_tmp=zero
    spinChi_iv_tmp=zero
    !
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
    !    
    if(.not.allocated(SpinChiMatrix(io,jo)%state)) then
       print*, "CHI_SPIN WARNING: SpinChiMatrix%state not allocated. Nothing to do"
       return
    endif
    !
    beta= 1d0/temp
    Egs = state_list%emin
    pesoBZ = 1d0/zeta_function
    !
    !this is the total number of available states  == state_list%size
    Nstates = size(SpinChiMatrix(io,jo)%state) 
    !Get trimmed state for the actual value of temp == state_list%trimd_size
    call es_trim_size(state_list,temp,cutoff) 
    do istate=1+MpiRank,state_list%trimd_size,MpiSize
       if(.not.allocated(SpinChiMatrix(io,jo)%state(istate)%channel))cycle
       Ei =  es_return_energy(state_list,istate)
       if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))/zeta_function
       Nchannels = size(SpinChiMatrix(io,jo)%state(istate)%channel)
       do ichan=1,Nchannels
          Nexcs  = size(SpinChiMatrix(io,jo)%state(istate)%channel(ichan)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = SpinChiMatrix(io,jo)%state(istate)%channel(ichan)%weight(iexc)
             peso  = peso*pesoBZ
             dE    = SpinChiMatrix(io,jo)%state(istate)%channel(ichan)%poles(iexc)
             ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
             ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
             ! This ensures that the correct null contribution is obtained.
             ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
             ! we do not include the contribution (because we are in the situation described above).
             ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
             ! if(beta*dE > 1d-3)spinChi_iv_tmp(io,jo,0)=spinChi_iv_tmp(io,jo,0) + 2*peso*(1d0-exp(-beta*dE))/dE 
             ! do i=1,Lmats
             !    spinChi_iv_tmp(io,jo,i)=spinChi_iv_tmp(io,jo,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
             ! enddo
             !Symmetrize for low-T /large-beta, mostly occurring for zero T calculations
             do i=0,Ltau/2-1
                spinChi_tau_tmp(io,jo,i)=spinChi_tau_tmp(io,jo,i) + peso*exp(-tau(i)*dE)
             enddo
             spinChi_tau_tmp(io,jo,Ltau/2)=spinChi_tau_tmp(io,jo,Ltau/2) + peso*0.5d0*(exp(-tau(Ltau/2)*dE)+exp(-(beta-tau(Ltau/2))*dE))
             do i=Ltau/2+1,Ltau
                spinChi_tau_tmp(io,jo,i)=spinChi_tau_tmp(io,jo,i) + peso*exp(-(beta-tau(i))*dE)
             enddo
             !
             do i=1,Lreal
                spinChi_w_tmp(io,jo,i)=spinChi_w_tmp(io,jo,i) - peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
             enddo
          enddo
       enddo
    enddo
#ifdef _MPI
    spinChi_tau=0d0
    spinChi_w=zero
    spinChi_iv=zero
    if(MpiStatus)then
       call AllReduce_Mpi(MpiComm,spinChi_tau_tmp,spinChi_tau)
       call AllReduce_Mpi(MpiComm,spinChi_w_tmp,spinChi_w)
       call AllReduce_Mpi(MpiComm,spinChi_iv_tmp,spinChi_iv) 
    else
       spinChi_tau=spinChi_tau_tmp
       spinChi_w=spinChi_w_tmp
       spinChi_iv=spinChi_iv_tmp
    endif
#else
    spinChi_tau=spinChi_tau_tmp
    spinChi_w=spinChi_w_tmp
    spinChi_iv=spinChi_iv_tmp
#endif
    return
  end subroutine lanc_eval_spinChi_electrons








  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################




  subroutine full_eval_spinChi_electrons(isite,jsite,iorb,jorb)
    integer      :: isite,jsite,iorb,jorb
    integer      :: io,jo
    type(sector) :: sectorI,sectorJ
    integer      :: Nups(1)
    integer      :: Ndws(1)
    real(8)      :: Chio,Chjo,Sio,Sjo
    integer      :: i,j,ll,m,isector
    integer      :: idim,ia
    real(8)      :: Ei,Ej,cc,peso,pesotot,beta
    real(8)      :: expterm,de,w0,it
    complex(8)   :: iw 
    !
    !
    !Spin susceptibility \X(tau). |<i|S_z|j>|^2
    !
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
    !
    beta= 1d0/temp
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
             Chio=0d0
             Chjo=0d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,sectorI%Dim
                call apply_op_Sz(i,Sio,io,sectorI)
                Chio   = Chio + espace(isector)%M(ll,i)*Sio*conjg(espace(isector)%M(ll,j))
                call apply_op_Sz(i,Sjo,jo,sectorI)
                Chjo   = Chjo + espace(isector)%M(ll,i)*Sjo*conjg(espace(isector)%M(ll,j))
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ei-Ej
             peso = Chio*Chjo/zeta_function
             !
             !Matsubara (bosonic) frequency
             if(beta*dE > 1d-3)spinChi_iv(io,jo,0)=spinChi_iv(io,jo,0) + peso*2*exp(-beta*Ej)*(1d0-exp(-beta*dE))/dE
             do m=1,Lmats
                spinChi_iv(io,jo,m)=spinChi_iv(io,jo,m)+ peso*exp(-beta*Ej)*2*dE/(vm(m)**2 + de**2)
             enddo
             !
             !Imaginary time: V
             do m=0,Ltau 
                it=tau(m)
                spinChi_tau(io,jo,m)=spinChi_tau(io,jo,m) + exp(-it*Ei)*exp(-(beta-it)*Ej)*peso
             enddo
             !
             !Real-frequency: Retarded = Commutator = response function
             do m=1,Lreal
                iw=dcmplx(vr(m),eps)
                spinChi_w(io,jo,m)=spinChi_w(io,jo,m)-peso*(exp(-beta*Ei) - exp(-beta*Ej))/(iw+de)
             enddo
             !
          enddo
       enddo
       call delete_sector(sectorI)
    enddo
  end subroutine full_eval_spinChi_electrons





END MODULE ED_CHI_SPIN_ELECTRONS
























