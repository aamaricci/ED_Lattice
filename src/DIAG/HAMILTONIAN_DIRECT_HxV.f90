! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec
#ifdef _MPI
  public  :: directMatVec_MPI
#endif



contains


  subroutine directMatVec(Nloc,vin,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    integer,dimension(2)                :: Indices,Jndices ![2-2*Norb]
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
  end subroutine directMatVec




#ifdef _MPI
  subroutine directMatVec_MPI(Nloc,vin,Hv)
    integer                             :: Nloc,N
    complex(8),dimension(Nloc)          :: Vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    integer,dimension(2)                :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns)               :: Nup,Ndw
    real(8),dimension(Ns)               :: Sz 
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
  end subroutine directMatVec_MPI
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
