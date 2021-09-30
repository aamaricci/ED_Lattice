! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
#ifdef _MPI
  public  :: directMatVec_MPI_main
#endif



contains


  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                           :: Nloc
    real(8),dimension(Nloc)           :: vin
    real(8),dimension(Nloc)           :: Hv
    real(8),dimension(:),allocatable  :: vt,Hvt
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)             :: Nup,Ndw
    complex(8),dimension(Nspin,Ns,Ns) :: Hij,Hloc
    complex(8),dimension(Nspin,Ns)    :: Hdiag

    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = diagonal(Hloc(ispin,:,:))
    enddo
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct/HxV_local.f90"
    if(ed_filling==0)then
       include "direct/Hxv_muhf.f90"
    endif
    !
    !UP HAMILTONIAN TERMS
    include "direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "direct/HxV_dw.f90"
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       include "direct/HxV_non_local.f90"
    endif
    !-----------------------------------------------!
    !
    return
  end subroutine directMatVec_main


#ifdef _MPI
  subroutine directMatVec_MPI_main(Nloc,vin,Hv)
    integer                           :: Nloc,N
    real(8),dimension(Nloc)           :: Vin
    real(8),dimension(Nloc)           :: Hv
    real(8),dimension(:),allocatable  :: vt,Hvt
    !
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)             :: Nup,Ndw
    complex(8),dimension(Nspin,Ns,Ns) :: Hij,Hloc
    complex(8),dimension(Nspin,Ns)    :: Hdiag
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
       Hdiag(ispin,:) = diagonal(Hloc(ispin,:,:))
    enddo
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct_mpi/HxV_local.f90"
    if(ed_filling==0)then
       include "direct_mpi/Hxv_muhf.f90"
    endif
    !
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "direct_mpi/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    allocate(vt(mpiQup*DimDw))
    allocate(Hvt(mpiQup*DimDw))
    vt=0d0
    Hvt=0d0
    !
    call vector_transpose_MPI(DimUp,MpiQdw,Vin(1:DimUp*MpiQdw),DimDw,MpiQup,vt) !Vin^T --> Vt
    include "direct_mpi/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=0d0         !reallocate Vt
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
    Hv(1:DimUp*MpiQdw) = Hv(1:DimUp*MpiQdw) + Vt
    !
    deallocate(vt,Hvt)
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = 0d0
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
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
