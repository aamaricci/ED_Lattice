! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
  public  :: directMatVec_kondo
#ifdef _MPI
  public  :: directMatVec_MPI_main
  public  :: directMatVec_MPI_kondo
#endif



contains


  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: Nup,Ndw
    real(8),dimension(Ns)               :: Sz 
    complex(8),dimension(Nspin,Ns,Ns)   :: Hij,Hloc
    complex(8),dimension(Nspin,Ns)      :: Hdiag

    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = dreal(diagonal(Hloc(ispin,:,:)))
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct/HxV_local.f90"
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    include "direct/HxV_non_local.f90"
    !
    !UP HAMILTONIAN TERMS
    include "direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "direct/HxV_dw.f90"
    !

    !-----------------------------------------------!
    !
    return
  end subroutine directMatVec_main



  subroutine directMatVec_kondo(Nloc,vin,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    complex(8),dimension(Nspin,Ns,Ns)   :: Hij,Hloc
    real(8),dimension(Nspin,Ns)         :: Hdiag
    integer,dimension(2*Ns_imp)         :: ib
    integer,dimension(Ns)               :: Nup,Ndw
    real(8),dimension(Ns)               :: Sz
    integer,dimension(Nimp)             :: NpUp,NpDw
    real(8),dimension(Nimp)             :: Szp
    integer                             :: i,j,io_up,io_dw,imp_up,imp_dw
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = dreal(diagonal(Hloc(ispin,:,:)))
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !
    do j=MpiIstart,MpiIend
       m   = Hsector%H(1)%map(j)
       ib  = bdecomp(m,2*Ns_imp)
       Nup = ib(1:Ns)
       Ndw = ib(Ns+1:2*Ns)
       NpUp= ib(2*Ns+1:2*Ns+Nimp)
       NpDw= ib(2*Ns+Nimp+1:2*Ns+2*Nimp)
       Sz  = 0.5d0*(Nup-Ndw)
       Szp = 0.5d0*(NpUp-NpDw)

       !LOCAL HAMILTONIAN TERMS
       include "direct/HxV_diag.f90"
       !
       !NON-LOCAL INTERACTION HAMILTONIAN TERMS
       include "direct/HxV_se_ph.f90"
       !
       !KONDO COUPLING HAMILTONIAN TERMS
       include "direct/HxV_kondo.f90"
       !
       !HOPPING TERMS
       include "direct/HxV_hop.f90"
       !
    enddo
    !-----------------------------------------------!
    !
    !
  end subroutine directMatVec_kondo




#ifdef _MPI
  subroutine directMatVec_MPI_main(Nloc,vin,Hv)
    integer                           :: Nloc,N
    complex(8),dimension(Nloc)           :: Vin
    complex(8),dimension(Nloc)           :: Hv
    complex(8),dimension(:),allocatable  :: vt,Hvt
    !
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)             :: Nup,Ndw
    real(8),dimension(Ns)             :: Sz 
    complex(8),dimension(Nspin,Ns,Ns) :: Hij,Hloc
    real(8),dimension(Nspin,Ns)       :: Hdiag
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index    
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = dreal(diagonal(Hloc(ispin,:,:)))
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct_mpi/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "direct_mpi/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    allocate(vt(mpiQup*DimDw))
    allocate(Hvt(mpiQup*DimDw))
    vt=zero
    Hvt=zero
    !
    call vector_transpose_MPI(DimUp,MpiQdw,Vin(1:DimUp*MpiQdw),DimDw,MpiQup,vt) !Vin^T --> Vt
    include "direct_mpi/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=zero         !reallocate Vt
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)  !Hvt^T --> Vt
    Hv(1:DimUp*MpiQdw) = Hv(1:DimUp*MpiQdw) + Vt
    !
    deallocate(vt,Hvt)
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = zero
       call allgather_vector_MPI(MpiComm,vin,vt)
       !
       include "direct_mpi/HxV_non_local.f90"
       !
       deallocate(Vt)
    endif
    !-----------------------------------------------!
    !
    return
  end subroutine directMatVec_MPI_main



  subroutine directMatVec_MPI_kondo(Nloc,v,Hv)
    integer                             :: Nloc,N
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vin
    !
    integer,dimension(2*Ns_imp)         :: ib
    integer,dimension(Ns)               :: Nup,Ndw
    real(8),dimension(Ns)               :: Sz
    integer,dimension(Nimp)             :: NpUp,NpDw
    real(8),dimension(Nimp)             :: Szp
    integer                             :: io_up,io_dw,imp_up,imp_dw
    complex(8),dimension(Nspin,Ns,Ns)   :: Hij,Hloc
    real(8),dimension(Nspin,Ns)         :: Hdiag
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index    
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = dreal(diagonal(Hloc(ispin,:,:)))
    enddo
    !
    !MPI part:
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    allocate(vin(N)) ; vin = zero
    call allgather_vector_MPI(MpiComm,v,vin)
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !
    states: do j=MpiIstart,MpiIend
       m   = Hsector%H(1)%map(j)
       ib  = bdecomp(m,2*Ns_imp)
       Nup = ib(1:Ns)
       Ndw = ib(Ns+1:2*Ns)
       NpUp= ib(2*Ns+1:2*Ns+Nimp)
       NpDw= ib(2*Ns+Nimp+1:2*Ns+2*Nimp)
       Sz  = 0.5d0*(Nup-Ndw)
       Szp = 0.5d0*(NpUp-NpDw)

       !LOCAL HAMILTONIAN TERMS
       include "direct/HxV_diag.f90"
       !
       !NON-LOCAL INTERACTION HAMILTONIAN TERMS
       include "direct/HxV_se_ph.f90"
       !
       !KONDO COUPLING HAMILTONIAN TERMS
       include "direct/HxV_kondo.f90"
       !
       !HOPPING TERMS
       include "direct/HxV_hop.f90"
       !
    enddo states
    !
    !-----------------------------------------------!
    deallocate(vin)
    return
  end subroutine directMatVec_MPI_kondo
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
