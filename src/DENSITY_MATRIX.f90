MODULE ED_DENSITY_MATRIX
  USE SF_CONSTANTS, only:zero,pi,xi
  !USE SF_IOTOOLS, only:free_unit,reg,txtfy
  !USE SF_ARRAYS, only: arange
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN
  implicit none
  private
  !
  public :: eval_dm_lattice


  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:,:),allocatable :: sz2,n2
  real(8),dimension(:,:),allocatable :: dens_ImpUp,dens_ImpDw
  real(8)                            :: dens_ph
  real(8)                            :: s2tot
  real(8)                            :: Egs
  real(8)                            :: Ei
  !
  integer                            :: iorb,jorb,iorb1,jorb1
  integer                            :: ispin,jspin
  integer                            :: isite,jsite
  integer                            :: io,jo,is,js
  integer                            :: ibath
  integer                            :: r,m,k,k1,k2,k3,k4
  integer                            :: iup,idw
  integer                            :: jup,jdw
  integer                            :: mup,mdw
  integer                            :: isectorDim
  real(8)                            :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  !
  integer                            :: i,j,ii
  integer                            :: isector,jsector
  !
  logical                            :: Jcondition
  !
  type(sector)                       :: sectorI,sectorJ


contains 


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the lattice density matrix
  !+-------------------------------------------------------------------+
  subroutine eval_dm_lattice()
    if(MPIMASTER)then
       write(LOGfile,"(A)")"Get lattice density matrix:"
       call start_timer()
    endif
    select case(ed_method)
    case default
       call lanc_density_matrix()
    case ("lapack","full")
       call full_density_matrix()
    end select
    if(MPIMASTER)call stop_timer(unit=LOGfile)
  end subroutine eval_dm_lattice






  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine lanc_density_matrix()
    integer                             :: iprob,istate,Nud(2,Ns),iud(2),jud(2),val,iimp
    integer,dimension(2)                :: Indices,Jndices
    integer,dimension(1,Ns)             :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: IbUp,IbDw,JbUp,JbDw  ![Ns]
    real(8),dimension(Ns)               :: nup,ndw,Sz,nt
    complex(8),dimension(:),allocatable :: state_cvec
    real(8)                             :: boltzman_weight
    real(8)                             :: state_weight
    !
    allocate(dens(Ns),dens_up(Ns),dens_dw(Ns))
    allocate(docc(Ns))
    allocate(magz(Ns),sz2(Ns,Ns),n2(Ns,Ns))
    allocate(dens_ImpUp(2,eNs))
    allocate(dens_ImpDw(2,eNs))
    !
    ed_dm_lattice = zero
    !
    call es_trim_size(state_list,temp,cutoff)
    do istate=1,state_list%trimd_size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
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
       !
       boltzman_weight = 1.d0 ; if(finiteT)boltzman_weight=exp(-(Ei-Egs)/temp)
       boltzman_weight = boltzman_weight/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i = 1,sectorI%Dim
             call build_op_Ns(i,IbUp,Ibdw,sectorI)
             do j = 1,sectorI%Dim
                call build_op_Ns(j,JbUp,JbDw,sectorI)
                ed_dm_lattice(i,j) = ed_dm_lattice(i,j) + &
                  state_cvec(i)*conjg(state_cvec(j))*boltzman_weight
             enddo
          enddo
          !
          call delete_sector(sectorI)
       endif
       !
       if(allocated(state_cvec))deallocate(state_cvec)      
       !
    enddo
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dm_lattice)
    endif
#endif
    !
  end subroutine lanc_density_matrix





  subroutine full_density_matrix()
    integer                         :: iprob,i,j
    integer                         :: izero,istate
    integer                         :: isector,jsector
    integer                         :: idim,jdim
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite
    integer                         :: r,m,k,val
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: boltzman_weight
    real(8)                         :: state_weight
    real(8)                         :: weight
    real(8)                         :: Ei
    real(8)                         :: temp_
    integer                         :: Nud(2,Ns),iud(2),jud(2)
    integer,dimension(2)      :: Indices,Jndices
    integer,dimension(1)        :: Nups,Ndws
    integer,dimension(Ns)       :: IbUp,IbDw,JbUp,JbDw ![Ns]
    real(8),dimension(Ns)       :: nup,ndw,Sz,nt
    complex(8),dimension(:),pointer :: evec
    !
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Ns),dens_up(Ns),dens_dw(Ns))
    allocate(docc(Ns))
    allocate(magz(Ns),sz2(Ns,Ns),n2(Ns,Ns))

    !
    egs     = gs_energy
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    temp_ = temp
    if(.not.finiteT)temp_=0.001d0
    !
    do isector=1,Nsectors
       call get_Nup(isector,nups)
       call get_Ndw(isector,ndws)
       if(ed_filling/=0 .AND. (sum(Nups)+sum(Ndws)/=ed_filling) )cycle

       call build_sector(isector,sectorI)
       !
       do istate=1,sectorI%Dim
          Ei=espace(isector)%e(istate)
          boltzman_weight=exp(-Ei/temp_)/zeta_function
          if(boltzman_weight < cutoff)cycle
          !
          evec => espace(isector)%M(:,istate)
          !
          do i = 1,sectorI%Dim
             call build_op_Ns(i,IbUp,Ibdw,sectorI)
             do j = 1,sectorI%Dim
                call build_op_Ns(j,JbUp,JbDw,sectorI)
                ed_dm_lattice(i,j) = ed_dm_lattice(i,j) + &
                  evec(i)*conjg(evec(j))*boltzman_weight
             enddo
          enddo
          !
       enddo
       call delete_sector(sectorI)
       if(associated(evec))nullify(evec)
    enddo
    !
  end subroutine full_density_matrix










  !####################################################################
  !                           I/O ROUTINES
  !####################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : write density matrices to file
  !+-------------------------------------------------------------------+



end MODULE ED_DENSITY_MATRIX
















