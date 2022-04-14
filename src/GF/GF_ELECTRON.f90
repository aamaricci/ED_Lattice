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


  public :: build_gf
  public :: eval_gf


  integer                             :: istate
  integer                             :: isector,jsector
  complex(8),allocatable              :: vvinit(:)
  real(8),allocatable                 :: alfa_(:),beta_(:)
  integer                             :: ialfa,jalfa
  integer                             :: i,j
  real(8)                             :: sgn,norm2
  real(8)                             :: state_e


contains



  !PURPOSE  : Build and Store the Green's functions weights and poles structure
  subroutine build_gf()
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
       if(KondoFlag.AND.gf_flag(Ns+1))then
          do ispin=1,Nspin
             do iimp=1,iNs
                if(MPIMASTER)call start_timer
                if(MPIMASTER)write(LOGfile,"(A)")"Build G:"//" imp "//str(iimp)//&
                     " spin"//str(ispin)
                call allocate_GFmatrix(impGmatrix(ispin,Ns+iimp,Ns+iimp),Nstate=state_list%size)
                call lanc_build_gf_impurity_diag(iimp,ispin)
                if(MPIMASTER)call stop_timer(unit=LOGfile)
             enddo
          enddo
       endif
       !Diagonal GF
       do ispin=1,Nspin
          do iorb=1,Norb
             if(.not.gf_flag(iorb))cycle
             do isite=1,Nsites(iorb)
                io    = pack_indices(isite,iorb)
                if(MPIMASTER)call start_timer
                if(MPIMASTER)write(LOGfile,"(A)")"Build G:"//&
                     " site I"//str(isite,site_indx_padding)//&
                     " orb M"//str(iorb)//&
                     " spin"//str(ispin)
                call allocate_GFmatrix(impGmatrix(ispin,io,io),Nstate=state_list%size)
                call lanc_build_gf_electrons_diag(isite,iorb,ispin)
                if(MPIMASTER)call stop_timer(unit=LOGfile)
             enddo
          enddo
       enddo
       !Off-diagonal GF
       if(offdiag_gf_flag)then
          do ispin=1,Nspin
             do iorb=1,Norb
                if(.not.gf_flag(iorb))cycle
                do jorb=1,Norb
                   if(.not.gf_flag(jorb))cycle
                   do isite=1,Nsites(iorb)
                      do jsite=1,Nsites(jorb)
                         io  = pack_indices(isite,iorb)
                         jo  = pack_indices(jsite,jorb)
                         if(io==jo)cycle
                         write(LOGfile,"(A)")"Build G:"//&
                              " sites I"//str(isite,site_indx_padding)//&
                              "J"//str(jsite,site_indx_padding)//&
                              " orb M"//str(iorb)//"L"//str(jorb)//&
                              " spin"//str(ispin)
                         if(MPIMASTER)call start_timer
                         call allocate_GFmatrix(impGmatrix(ispin,io,jo),Nstate=state_list%size)
                         call lanc_build_gf_electrons_mix(isite,jsite,iorb,jorb,ispin)
                         if(MPIMASTER)call stop_timer(unit=LOGfile)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       end if
       !
    end select
  end subroutine build_gf





  !Evaluate the Green's functions from knowledge of the weights+poles for a given case
  subroutine eval_gf()
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

    do ispin=1,Nspin
       do iorb=1,Norb
          do isite=1,Nsites(iorb)
             if(MPIMASTER)write(LOGfile,"(A)")"Eval G:"//&
                  " site I"//str(isite,site_indx_padding)//&
                  " orb M"//str(iorb)//&
                  " spin"//str(ispin)
             io    = pack_indices(isite,iorb)
             if(MPIMASTER)call start_timer
             select case(ed_method)
             case default
                call lanc_eval_gf_electrons(isite,isite,iorb,iorb,ispin)
             case ('lapack','full')
                call full_eval_gf_electrons_diag(isite,iorb,ispin)
             end select
             if(MPIMASTER)call stop_timer(unit=LOGfile)
          enddo
       enddo
    enddo
    !
    if(offdiag_gf_flag)then     
       do ispin=1,Nspin
          do iorb=1,Norb
             if(.not.gf_flag(iorb))cycle
             do jorb=1,Norb
                if(.not.gf_flag(jorb))cycle
                do isite=1,Nsites(iorb)
                   do jsite=1,Nsites(jorb)
                      io  = pack_indices(isite,iorb)
                      jo  = pack_indices(jsite,jorb)
                      if(io==jo)cycle
                      write(LOGfile,"(A)")"Eval G:"//&
                           " sites I"//str(isite,site_indx_padding)//&
                           "J"//str(jsite,site_indx_padding)//&
                           " orb M"//str(iorb)//"L"//str(jorb)//&
                           " spin"//str(ispin)
                      if(MPIMASTER)call start_timer
                      select case(ed_method)
                      case default
                         call lanc_eval_gf_electrons(isite,jsite,iorb,jorb,ispin)
                      case ('lapack','full')
                         call full_eval_gf_electrons_mix(isite,jsite,iorb,jorb,ispin)
                      end select
                      if(MPIMASTER)call stop_timer(unit=LOGfile)
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
       !
       !Put here off-diagonal manipulation by symmetry:
       select case(ed_method)
       case default
          do ispin=1,Nspin
             do iorb=1,Norb
                if(.not.gf_flag(iorb))cycle
                do jorb=1,Norb
                   if(.not.gf_flag(jorb))cycle
                   do isite=1,Nsites(iorb)
                      do jsite=1,Nsites(jorb)
                         io  = pack_indices(isite,iorb)
                         jo  = pack_indices(jsite,jorb)
                         if(io==jo)cycle
                         impGmats(ispin,io,jo,:) = 0.5d0*(impGmats(ispin,io,jo,:) - impGmats(ispin,io,io,:) - impGmats(ispin,jo,jo,:))
                         impGreal(ispin,io,jo,:) = 0.5d0*(impGreal(ispin,io,jo,:) - impGreal(ispin,io,io,:) - impGreal(ispin,jo,jo,:))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       case ('lapack','full')
          continue
       end select
    end if
    !
  end subroutine eval_gf






  !################################################################
  !################################################################
  !################################################################
  !################################################################





  subroutine lanc_build_gf_electrons_diag(isite,iorb,ispin)
    integer,intent(in)                  :: isite,iorb,ispin
    integer                             :: io
    type(sector)                        :: sectorI,sectorJ
    complex(8),dimension(:),allocatable :: state_cvec
    !
    ialfa = 1
    io    = pack_indices(isite,iorb)
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
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c^+_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,io,ialfa,ispin,sectorI,sectorJ)
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
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,io,ialfa,ispin,sectorI,sectorJ)
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
  end subroutine lanc_build_gf_electrons_diag




  !################################################################




  subroutine lanc_build_gf_impurity_diag(iimp,ispin)
    integer,intent(in)                  :: iimp,ispin
    integer                             :: io,ipos
    type(sector)                        :: sectorI,sectorJ
    complex(8),dimension(:),allocatable :: state_cvec
    !
    ialfa = 1
    io    = Ns + iimp
    ipos  = 2*Ns + iimp + (ispin-1)*iNs
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
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c^+_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,ipos,ialfa,1,sectorI,sectorJ)
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
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,ipos,ialfa,1,sectorI,sectorJ)
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



  subroutine lanc_build_gf_electrons_mix(isite,jsite,iorb,jorb,ispin)
    integer                             :: isite,jsite,iorb,jorb,ispin
    integer                             :: io,jo
    type(sector)                        :: sectorI,sectorJ
    complex(8),dimension(:),allocatable :: state_cvec
    !
    ialfa = 1
    jalfa = ialfa
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(ispin,io,jo),istate,Nchan=2)
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
                call apply_op_CDG(i,j,sgn,io,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             !+c^+_jo|gs>
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,jo,jalfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,io,jo,ispin,ichan=1,istate=istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,io,jo),istate,1,Nexc=0)
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
                call apply_op_C(i,j,sgn,io,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             !+c_jo|gs>
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,jo,jalfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,io,jo,ispin,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,io,jo),istate,2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    return
  end subroutine lanc_build_gf_electrons_mix





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






  subroutine lanc_eval_gf_electrons(isite,jsite,iorb,jorb,ispin)
    integer,intent(in) :: isite,jsite,iorb,jorb,ispin
    integer            :: Nstates,istate
    integer            :: Nchannels,ichan
    integer            :: Nexcs,iexc,io,jo
    real(8)            :: peso,de,pesoBZ,Ei,Egs,beta
    !
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
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
  end subroutine lanc_eval_gf_electrons





  subroutine lanc_eval_gf_impurity(iimp,jimp,ispin)
    integer,intent(in) :: iimp,jimp,ispin
    integer            :: Nstates,istate
    integer            :: Nchannels,ichan
    integer            :: Nexcs,iexc,io,jo
    real(8)            :: peso,de,pesoBZ,Ei,Egs,beta
    !
    io = Ns + iimp
    jo = Ns + jimp
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






  subroutine full_eval_gf_electrons_diag(isite,iorb,ispin)
    integer,intent(in) :: isite,iorb,ispin
    integer            :: io
    type(sector)       :: sectorI,sectorJ
    integer            :: Nups(Ns_Ud)
    integer            :: Ndws(Ns_Ud)
    complex(8)         :: op_mat(2)
    real(8)            :: spectral_weight
    real(8)            :: sgn_cdg,sgn_c
    integer            :: m,i,j,li,rj
    real(8)            :: Ei,Ej
    real(8)            :: expterm,peso,de,w0
    complex(8)         :: iw
    !    
    !
    ialfa = 1
    io    = pack_indices(isite,iorb)
    !
    do isector=1,Nsectors
       call get_Nup(isector,nups)
       call get_Ndw(isector,ndws)
       if(ed_filling/=0 .AND. (abs(sum(Nups)+sum(Ndws)-ed_filling)>1) )cycle
       !
       jsector=getCDGsector(ialfa,ispin,isector)
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
                call apply_op_CDG(li,rj,sgn_cdg,io,ialfa,ispin,sectorI,sectorJ)
                if(sgn_cdg==0d0.OR.rj==0)cycle
                !
                op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
             enddo
             !
             do rj=1,sectorJ%Dim
                call apply_op_C(rj,li,sgn_c,io,ialfa,ispin,sectorJ,sectorI)
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
  end subroutine full_eval_gf_electrons_diag


  subroutine full_eval_gf_impurity(iimp,ispin)
    integer,intent(in) :: iimp,ispin
    integer            :: io,ipos
    type(sector)       :: sectorI,sectorJ
    integer            :: Nups(Ns_Ud)
    integer            :: Ndws(Ns_Ud)
    complex(8)         :: op_mat(2)
    real(8)            :: spectral_weight
    real(8)            :: sgn_cdg,sgn_c
    integer            :: m,i,j,li,rj
    real(8)            :: Ei,Ej
    real(8)            :: expterm,peso,de,w0
    complex(8)         :: iw
    !    
    !
    ialfa = 1
    io    = Ns+iimp
    ipos  = 2*Ns + iimp + (ispin-1)*iNs
    !
    do isector=1,Nsectors
       call get_Nup(isector,nups)
       call get_Ndw(isector,ndws)
       if(ed_filling/=0 .AND. (abs(sum(Nups)+sum(Ndws)-ed_filling)>1) )cycle
       !
       jsector=getCDGsector(ialfa,ispin,isector)
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
                call apply_op_CDG(li,rj,sgn_cdg,ipos,ialfa,1,sectorI,sectorJ)
                if(sgn_cdg==0d0.OR.rj==0)cycle
                !
                op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
             enddo
             !
             do rj=1,sectorJ%Dim
                call apply_op_C(rj,li,sgn_c,ipos,ialfa,1,sectorJ,sectorI)
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



  subroutine full_eval_gf_electrons_mix(isite,jsite,iorb,jorb,ispin)
    integer      :: isite,jsite,iorb,jorb,ispin
    integer      :: io,jo
    type(sector) :: sectorI,sectorJ
    integer      :: Nups(Ns_Ud)
    integer      :: Ndws(Ns_Ud)
    complex(8)   :: op_mat(2)
    complex(8)   :: spectral_weight
    real(8)      :: sgn_cdg,sgn_c
    integer      :: m,i,j,li,rj
    real(8)      :: Ei,Ej
    real(8)      :: expterm,peso,de,w0
    complex(8)   :: iw
    !
    ialfa = 1
    jalfa = ialfa
    io  = pack_indices(isite,iorb)
    jo  = pack_indices(jsite,jorb)
    !
    do isector=1,Nsectors
       call get_Nup(isector,nups)
       call get_Ndw(isector,ndws)
       !
       if(ed_filling/=0 .AND. (sum(Nups)+sum(Ndws)/=ed_filling) )cycle
       !
       jsector=getCDGsector(ialfa,ispin,isector)
       if(jsector==0)cycle
       !
       call build_sector(isector,sectorI)
       call build_sector(jsector,sectorJ)
       !
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
                call apply_op_CDG(li,rj,sgn_cdg,io,ialfa,ispin,sectorI,sectorJ)
                if(sgn_cdg==0d0.OR.rj==0)cycle
                !
                op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
             enddo
             !
             do rj=1,sectorJ%Dim
                call apply_op_C(rj,li,sgn_c,jo,jalfa,ispin,sectorJ,sectorI)
                if(sgn_c==0d0.OR.li==0)cycle
                !
                op_mat(2)=op_mat(2) + conjg(espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
             enddo
             !
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei
             peso=expterm/zeta_function
             spectral_weight=peso*product(op_mat)
             !
             do m=1,Lmats
                iw=xi*wm(m)
                impGmats(ispin,io,jo,m)=impGmats(ispin,io,jo,m)+spectral_weight/(iw-de)
             enddo
             !
             do m=1,Lreal
                w0=wr(m);iw=cmplx(w0,eps)
                impGreal(ispin,io,jo,m)=impGreal(ispin,io,jo,m)+spectral_weight/(iw-de)
             enddo
             !
          enddo
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    enddo
  end subroutine full_eval_gf_electrons_mix






  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################










END MODULE ED_GF_ELECTRON













! subroutine eval_sigma_normal
!   integer                                 :: i,ispin,iorb
!   complex(8),dimension(Nspin,Ns,Ns,Lmats) :: invG0mats,invGmats
!   complex(8),dimension(Nspin,Ns,Ns,Lreal) :: invG0real,invGreal
!   complex(8),dimension(Ns,Ns)             :: Gtmp
!   !
!   invG0mats = zero
!   invGmats  = zero
!   invG0real = zero
!   invGreal  = zero
!   !
!   !Get G0^-1
!   call Hij_get_g0inv(dcmplx(0d0,wm(:)),invG0mats)
!   call Hij_get_g0inv(dcmplx(wr(:),eps),invG0real)
!   !
!   impSmats=zero
!   impSreal=zero
!   !Get Gimp^-1
!   do ispin=1,Nspin
!      do i=1,Lmats
!         Gtmp = impGmats(ispin,:,:,i)
!         call inv(Gtmp)
!         invGmats(ispin,:,:,i)=Gtmp
!      enddo
!      !
!      do i=1,Lreal
!         Gtmp = impGreal(ispin,:,:,i)
!         call inv(Gtmp)
!         invGreal(ispin,:,:,i)=Gtmp
!      enddo
!      !
!      !Get Sigma functions: Sigma= G0^-1 - G^-1
!      impSmats(ispin,:,:,:) = invG0mats(ispin,:,:,:) - invGmats(ispin,:,:,:)
!      impSreal(ispin,:,:,:) = invG0real(ispin,:,:,:) - invGreal(ispin,:,:,:)
!   enddo
!   !
!   !Get G0and:
!   call Hij_get_g0func(dcmplx(0d0,wm(:)),impG0mats)
!   call Hij_get_g0func(dcmplx(wr(:),eps),impG0real)
!   !
! end subroutine eval_sigma_normal











