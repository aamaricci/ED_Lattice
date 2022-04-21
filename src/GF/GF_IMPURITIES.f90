MODULE ED_GF_IMPURITIES
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


  public :: build_gf_impurities
  public :: eval_gf_impurities


  integer                             :: istate
  integer                             :: isector,jsector
  complex(8),allocatable              :: vvinit(:)
  real(8),allocatable                 :: alfa_(:),beta_(:)
  integer                             :: i,j
  real(8)                             :: sgn,norm2
  real(8)                             :: state_e


contains



  !PURPOSE  : Build and Store the Green's functions weights and poles structure
  subroutine build_gf_impurities()
    integer :: ispin,i,iimp
    integer :: iorb,jorb
    integer :: isite,jsite
    integer :: io,jo
    !
    call deallocate_GFmatrix(impGmatrix)
    !
    select case(ed_method)
    case ('lapack','full')
       !do nothing
       return
    case default
       !Impurity GF
       if(KondoFlag.AND.gf_flag(Norb+1))then
          do ispin=1,Nspin
             do iimp=1,iNs
                if(MPIMASTER)call start_timer
                if(MPIMASTER)write(LOGfile,"(A)")"Build G:"//" imp "//str(iimp)//&
                     " spin"//str(ispin)
                call allocate_GFmatrix(impGmatrix(ispin,eNs+iimp,eNs+iimp),Nstate=state_list%size)
                call lanc_build_gf_impurity_diag(iimp,ispin)
                if(MPIMASTER)call stop_timer(unit=LOGfile)
             enddo
          enddo
       endif
    end select
  end subroutine build_gf_impurities






  !Evaluate the Green's functions from knowledge of the weights+poles for a given case
  subroutine eval_gf_impurities()
    integer :: ispin,i,iimp
    integer :: iorb,jorb
    integer :: isite,jsite
    integer :: io,jo
    !
    if(KondoFlag)then
       do ispin=1,Nspin
          do iimp=1,iNs
             if(MPIMASTER)write(LOGfile,"(A)")"Eval G:"//" imp"//str(iimp)//&
                  " spin"//str(ispin)
             if(MPIMASTER)call start_timer
             select case(ed_method)
             case default
                call lanc_eval_gf_impurity(iimp,iimp,ispin)
             case ('lapack','full')
                call full_eval_gf_impurity(iimp,ispin)
             end select
             if(MPIMASTER)call stop_timer(unit=LOGfile)
          enddo
       enddo
    endif
    !
  end subroutine eval_gf_impurities






  !################################################################
  !################################################################
  !################################################################
  !################################################################





  subroutine lanc_build_gf_impurity_diag(iimp,ispin)
    integer,intent(in)                  :: iimp,ispin
    integer                             :: io,ipos
    type(sector)                        :: sectorI,sectorJ
    complex(8),dimension(:),allocatable :: state_cvec
    !
    io    = eNs + iimp
    ipos  = ket_imp_index(iimp,ispin)
    !
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(ispin,io,io),istate,Nchan=2) !2=particle/hole exc
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
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c^+_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,ipos,1,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,io,io,ispin,ichan=1,istate=istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,io,io),istate,1,Nexc=0)
       endif
       !
       !REMOVE ONE PARTICLE:
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,ipos,1,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,io,io,ispin,ichan=2,istate=istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,io,io),istate,2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    return
  end subroutine lanc_build_gf_impurity_diag



  !################################################################

  

  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,io,jo,ispin,ichan,istate)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,io,jo,ispin,ichan,istate
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
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
    call allocate_GFmatrix(impGmatrix(ispin,io,jo),istate,ichan,Nlanc)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = Z(1,j)*Z(1,j)*vnorm2
       impGmatrix(ispin,io,jo)%state(istate)%channel(ichan)%weight(j) = peso
       impGmatrix(ispin,io,jo)%state(istate)%channel(ichan)%poles(j)  = isign*de
    enddo
  end subroutine add_to_lanczos_gf_normal




  !################################################################
  !################################################################
  !################################################################
  !################################################################


  subroutine lanc_eval_gf_impurity(iimp,jimp,ispin)
    integer,intent(in) :: iimp,jimp,ispin
    integer            :: Nstates,istate
    integer            :: Nchannels,ichan
    integer            :: Nexcs,iexc,io,jo
    real(8)            :: peso,de,pesoBZ,Ei,Egs,beta
    !
    io = eNs + iimp
    jo = eNs + jimp
    !    
    if(.not.allocated(impGmatrix(ispin,io,jo)%state)) then
       print*, "GF_NORMAL WARNING: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    beta= 1d0/temp
    Egs = state_list%emin
    pesoBZ = 1d0/zeta_function
    !
    !this is the total number of available states  == state_list%size
    Nstates = size(impGmatrix(ispin,io,jo)%state) 
    !Get trimmed state for the actual value of temp == state_list%trimd_size
    call es_trim_size(state_list,temp,cutoff) 
    do istate=1,state_list%trimd_size     !Nstates
       if(.not.allocated(impGmatrix(ispin,io,jo)%state(istate)%channel))cycle
       Ei =  es_return_energy(state_list,istate)
       if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))/zeta_function
       Nchannels = size(impGmatrix(ispin,io,jo)%state(istate)%channel)
       do ichan=1,Nchannels
          Nexcs  = size(impGmatrix(ispin,io,jo)%state(istate)%channel(ichan)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = impGmatrix(ispin,io,jo)%state(istate)%channel(ichan)%weight(iexc)
             de    = impGmatrix(ispin,io,jo)%state(istate)%channel(ichan)%poles(iexc)
             impGmats(ispin,io,jo,:)=impGmats(ispin,io,jo,:) + pesoBZ*peso/(dcmplx(0d0,wm(:))-de)
             impGreal(ispin,io,jo,:)=impGreal(ispin,io,jo,:) + pesoBZ*peso/(dcmplx(wr(:),eps)-de)
          enddo
       enddo
    enddo
    return
  end subroutine lanc_eval_gf_impurity


  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################



  subroutine full_eval_gf_impurity(iimp,ispin)
    integer,intent(in) :: iimp,ispin
    integer            :: io,ipos
    type(sector)       :: sectorI,sectorJ
    integer            :: Nups(1)
    integer            :: Ndws(1)
    complex(8)         :: op_mat(2)
    real(8)            :: spectral_weight
    real(8)            :: sgn_cdg,sgn_c
    integer            :: m,i,j,li,rj
    real(8)            :: Ei,Ej
    real(8)            :: expterm,peso,de,w0
    complex(8)         :: iw
    !    
    !
    io    = eNs+iimp
    ipos  = ket_imp_index(iimp,ispin)!2*eNs + iimp + (ispin-1)*iNs
    !
    do isector=1,Nsectors
       call get_Nup(isector,nups)
       call get_Ndw(isector,ndws)
       if(ed_filling/=0 .AND. (abs(sum(Nups)+sum(Ndws)-ed_filling)>1) )cycle
       !
       jsector=getCDGsector(1,ispin,isector)
       if(jsector==0)cycle
       !
       call build_sector(isector,sectorI)
       call build_sector(jsector,sectorJ)
       !
       do i=1,sectorI%Dim          !loop over the states in the i-th sect.
          do j=1,sectorJ%Dim       !loop over the states in the j-th sect.
             !
             expterm=exp(-espace(isector)%e(i)/temp)+exp(-espace(jsector)%e(j)/temp)
             if(expterm < cutoff)cycle
             !
             op_mat=0d0
             !
             do li=1,sectorI%Dim              !loop over the component of |I> (IN state!)
                call apply_op_CDG(li,rj,sgn_cdg,ipos,1,sectorI,sectorJ)
                if(sgn_cdg==0d0.OR.rj==0)cycle
                !
                op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
             enddo
             !
             do rj=1,sectorJ%Dim
                call apply_op_C(rj,li,sgn_c,ipos,1,sectorJ,sectorI)
                if(sgn_c==0d0.OR.li==0)cycle
                !
                op_mat(2)=op_mat(2) + conjg(espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
             enddo
             !
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei;
             peso=expterm/zeta_function
             spectral_weight=peso*product(op_mat)
             !
             do m=1,Lmats
                iw=xi*wm(m)
                impGmats(ispin,io,io,m)=impGmats(ispin,io,io,m)+spectral_weight/(iw-de)
             enddo
             !
             do m=1,Lreal
                w0=wr(m);iw=cmplx(w0,eps)
                impGreal(ispin,io,io,m)=impGreal(ispin,io,io,m)+spectral_weight/(iw-de)
             enddo
             !
          enddo
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    enddo
  end subroutine full_eval_gf_impurity



END MODULE ED_GF_IMPURITIES



















