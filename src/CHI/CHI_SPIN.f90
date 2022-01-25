MODULE ED_CHI_SPIN
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


  public :: build_chi_spin

  integer                         :: istate,iorb,jorb,ispin
  integer                         :: isector
  complex(8),allocatable          :: vvinit(:)
  real(8),allocatable             :: alfa_(:),beta_(:)
  integer                         :: ialfa
  integer                         :: jalfa
  integer                         :: ipos,jpos
  integer                         :: i,j
  real(8)                         :: sgn,norm2
  complex(8),dimension(:),pointer :: state_cvec
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
  subroutine build_chi_spin()
    integer :: ispin,i
    integer :: iorb,jorb
    integer :: isite,jsite
    integer :: io,jo
    !
    do iorb=1,Norb
       do isite=1,Nsites(iorb)
          write(LOGfile,"(A)")"Get spinChi:"//&
               " site I"//str(isite,site_indx_padding)//&
               " orb M"//str(iorb)
          if(MPIMASTER)call start_timer
          select case(ed_method)
          case default
             call lanc_ed_build_spinChi_diag(isite,iorb)
          case ('lapack','full')
             call full_ed_build_spinChi_main(isite,isite,iorb,iorb)
          end select
          if(MPIMASTER)call stop_timer(unit=LOGfile)
       enddo
    enddo
    !
    if(offdiag_chispin_flag.AND.Norb>1)then
       do iorb=1,Norb
          do jorb=1,Norb
             do isite=1,Nsites(iorb)
                do jsite=1,Nsites(jorb)
                   io  = pack_indices(isite,iorb)
                   jo  = pack_indices(jsite,jorb)
                   if(io==jo)cycle
                   write(LOGfile,"(A)")"Get spinChi:"//&
                        " sites I"//str(isite,site_indx_padding)//&
                        "J"//str(jsite,site_indx_padding)//&
                        " orb M"//str(iorb)//"L"//str(jorb)
                   if(MPIMASTER)call start_timer
                   select case(ed_method)
                   case default
                      call lanc_ed_build_spinChi_mix(isite,jsite,iorb,jorb)
                   case ('lapack','full')
                      call full_ed_build_spinChi_main(isite,isite,iorb,iorb)
                   end select
                   if(MPIMASTER)call stop_timer(unit=LOGfile)
                enddo
             enddo
          enddo
       enddo
       !
       !
       do io=1,Ns
          do jo=1,Ns
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
    endif
    !
  end subroutine build_chi_spin






  !################################################################
  !################################################################
  !################################################################
  !################################################################






     subroutine lanc_ed_build_spinChi_diag(isite,iorb)
       integer,intent(in) :: isite,iorb
       integer            :: io
       type(sector)       :: sectorI,sectorJ
       !
       ialfa = 1
       io    = pack_indices(isite,iorb)
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
                  'Apply Sz  :',isector,sectorI%Nups,sectorI%Ndws
             allocate(vvinit(sectorI%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_Sz(i,sgn,io,ialfa,sectorI)            
                vvinit(i) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorI)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,io,io)
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
     end subroutine lanc_ed_build_spinChi_diag




  !################################################################




  subroutine lanc_ed_build_spinChi_mix(isite,jsite,iorb,jorb)
    integer      :: isite,jsite,iorb,jorb
    integer      :: io,jo
    type(sector) :: sectorI,sectorJ
    real(8)      :: Siorb,Sjorb
    !
    ialfa = 1
    jalfa = ialfa
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
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
       !EVALUATE (Sz_jorb + Sz_iorb)|gs> = Sz_jorb|gs> + Sz_iorb|gs>
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
          if(ed_verbose==3)write(LOGfile,"(A,I15)")'Apply (Sz_jorb + Sz_iorb):',isector
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          do i=1,sectorI%Dim
             call apply_op_Sz(i,Siorb,io,ialfa,sectorI)
             call apply_op_Sz(i,Sjorb,jo,jalfa,sectorI)
             sgn       = Siorb + Sjorb
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=zero
       endif
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,io,jo)
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
  end subroutine lanc_ed_build_spinChi_mix





  !################################################################


  subroutine add_to_lanczos_spinChi(vnorm2,Ei,alanc,blanc,io,jo)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: io,jo
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
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
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)spinChi_iv(io,jo,0)=spinChi_iv(io,jo,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          spinChi_iv(io,jo,i)=spinChi_iv(io,jo,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       !Symmetrize for low-T /large-beta, mostly occurring for zero T calculations
       do i=0,Ltau/2-1
          spinChi_tau(io,jo,i)=spinChi_tau(io,jo,i) + peso*exp(-tau(i)*dE)
       enddo
       spinChi_tau(io,jo,Ltau/2)=spinChi_tau(io,jo,Ltau/2) + peso*0.5d0*(exp(-tau(Ltau/2)*dE)+exp(-(beta-tau(Ltau/2))*dE))
       do i=Ltau/2+1,Ltau
          spinChi_tau(io,jo,i)=spinChi_tau(io,jo,i) + peso*exp(-(beta-tau(i))*dE)
       enddo
       !
       do i=1,Lreal
          spinChi_w(io,jo,i)=spinChi_w(io,jo,i) - peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_spinChi





  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine full_ed_build_spinChi_main(isite,jsite,iorb,jorb)
    integer      :: isite,jsite,iorb,jorb
    integer      :: io,jo
    type(sector) :: sectorI,sectorJ
    integer            :: Nups(Ns_Ud)
    integer            :: Ndws(Ns_Ud)
    real(8)      :: Chio,Chjo,Sio,Sjo
    integer      :: i,j,ll,m,isector
    integer      :: idim,ia
    real(8)      :: Ei,Ej,cc,peso,pesotot
    real(8)      :: expterm,de,w0,it
    complex(8)   :: iw 
    !
    !
    !Spin susceptibility \X(tau). |<i|S_z|j>|^2
    !
    ialfa = 1
    jalfa = ialfa
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
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
                call apply_op_Sz(i,Sio,io,ialfa,sectorI)
                Chio   = Chio + espace(isector)%M(ll,i)*Sio*conjg(espace(isector)%M(ll,j))
                call apply_op_Sz(i,Sjo,jo,jalfa,sectorI)
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
  end subroutine full_ed_build_spinChi_main





END MODULE ED_CHI_SPIN
























