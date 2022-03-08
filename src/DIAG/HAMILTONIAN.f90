MODULE ED_HAMILTONIAN
  USE ED_HAMILTONIAN_COMMON
  USE ED_HAMILTONIAN_STORED_HxV
  USE ED_HAMILTONIAN_DIRECT_HxV
  !
  implicit none
  private


  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector
  public  :: delete_Hv_sector
  public  :: vecDim_Hv_sector

  !> Tridiag sparse Hamiltonian of the sector
  public  :: tridiag_Hv_sector

  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_main
  public  :: spMatVec_kondo
#ifdef _MPI
  public  :: spMatVec_MPI_main
  public  :: spMatVec_MPI_kondo
#endif


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
  public  :: directMatVec_kondo
#ifdef _MPI
  public  :: directMatVec_MPI_main
  public  :: directMatVec_MPI_kondo
#endif




contains



  !####################################################################
  !                 MAIN ROUTINES: BUILD SECTOR
  !####################################################################
  subroutine build_Hv_sector(isector,Hmat)
    integer                            :: isector
    complex(8),dimension(:,:),optional :: Hmat   
    select case(KondoFlag)       
    case(.true.)
       if(present(Hmat))then
          call build_Hv_sector_kondo(isector,Hmat)
       else
          call build_Hv_sector_kondo(isector)
       endif
    case(.false.)
       if(present(Hmat))then
          call build_Hv_sector_main(isector,Hmat)
       else
          call build_Hv_sector_main(isector)
       endif
    end select
  end subroutine build_Hv_sector



  !####################################################################
  !                 MAIN ROUTINES: DELETE SECTOR
  !####################################################################
  subroutine delete_Hv_sector()
    select case(KondoFlag)
    case (.true.)
       call delete_Hv_sector_kondo()
    case(.false.)
       call delete_Hv_sector_main()
    end select
  end subroutine delete_Hv_sector



  !####################################################################
  !                 MAIN ROUTINES: DELETE SECTOR
  !####################################################################
  function vecDim_Hv_sector(isector) result(vecDim)
    integer :: isector
    integer :: vecDim
    select case(KondoFlag)
    case (.true.)
       vecDim = vecDim_Hv_sector_kondo(isector)
    case(.false.)
       vecDim = vecDim_Hv_sector_main(isector)
    end select
  end function vecDim_Hv_sector


  !####################################################################
  !                 MAIN ROUTINES: DELETE SECTOR
  !####################################################################
  subroutine tridiag_Hv_sector(isector,vvinit,alanc,blanc,norm2)
    integer                          :: isector
    complex(8),dimension(:)          :: vvinit
    real(8),dimension(:),allocatable :: alanc,blanc
    real(8)                          :: norm2
    integer                          :: vecDim
    select case(KondoFlag)
    case (.true.)
       call tridiag_Hv_sector_kondo(isector,vvinit,alanc,blanc,norm2)
    case(.false.)
       call tridiag_Hv_sector_main(isector,vvinit,alanc,blanc,norm2)
    end select
  end subroutine tridiag_Hv_sector










  !####################################################################
  !####################################################################
  !                          MAIN 
  !####################################################################
  !####################################################################

  subroutine build_Hv_sector_main(isector,Hmat)
    integer                            :: isector
    complex(8),dimension(:,:),optional :: Hmat   
    integer                            :: irank,ierr
    integer                            :: i,iup,idw
    integer                            :: j,jup,jdw
    !
    call build_sector(isector,Hsector)
    !
    !This is not really needed but it eases the writing:
    allocate(DimUps(Ns_Ud))
    allocate(DimDws(Ns_Ud))
    Dim    = Hsector%Dim
    DimUp  = Hsector%DimUp
    DimDw  = Hsector%DimDw
    DimUps = Hsector%DimUps
    DimDws = Hsector%DimDws
    !
    !
    !#################################
    !          MPI SETUP
    !#################################
    mpiAllThreads=.true.
    !>PREAMBLE: check that split of the DW is performed with the minimum #cpu: no idle cpus allowed (with zero elements)
#ifdef _MPI
    if(MpiStatus)then
       if(DimDw < MpiSize)then
          if(MpiMaster.AND.ed_verbose>4)write(*,*)"Reducing N_cpu to DimDw:",DimDw,MpiSize-DimDw
          allocate(MpiMembers(0:DimDw-1))
          forall(irank=0:DimDw-1)MpiMembers(irank)=irank       
          call Mpi_Group_Incl(MpiGroup_Global,DimDw,MpiMembers,MpiGroup,ierr)
          call Mpi_Comm_create(MpiComm_Global,MpiGroup,MpiComm,ierr)
          deallocate(MpiMembers)
          mpiAllThreads=.false.
          call Barrier_MPI(MpiComm_Global)
#ifdef _DEBUG
          if(ed_verbose>4)then
             if(MpiMaster)write(LOGfile,*)&
                  "       mpiRank,   MpiComm, Comm_Global, Comm_World, Comm_Null, Undefined"
             do i=0,MpiSize-1
                call Barrier_MPI(MpiComm_Global)
                if(MpiRank==i)write(*,*)i,MpiComm,MpiComm_Global,Mpi_Comm_World,Mpi_comm_null,Mpi_Undefined
             enddo
             call Barrier_MPI(MpiComm_Global)
          endif
#endif
       endif
       if( MpiComm /= MPI_COMM_NULL )then
          MpiRank = Get_Rank_MPI(MpiComm)
          MpiSize = Get_Size_MPI(MpiComm)
       endif
    endif
#endif
    !
    !Dw split:    
    mpiQdw = DimDw/MpiSize
    mpiRdw = mod(DimDw,MpiSize)
    if(MpiRank < mod(DimDw,MpiSize))then
       !Total split: split DW \times UP 
       mpiRdw = 0
       MpiQdw = MpiQdw+1
    endif
    !
    !Total split: split DW \times UP
    mpiQ = DimUp*mpiQdw
    mpiR = DimUp*mpiRdw
    mpiIstart = 1 + MpiRank*mpiQ+mpiR
    mpiIend   = (MpiRank+1)*mpiQ+mpiR
    mpiIshift = MpiRank*mpiQ+mpiR
    !
    !
#ifdef _MPI
#ifdef _DEBUG
    if(MpiStatus.AND.ed_verbose>4.AND.(MpiComm/=Mpi_Comm_Null).AND.MpiSize>=1)then
       if(MpiMaster)write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,      mpi_Qdw,      mpiR_dw,  mpi_Istart,  mpi_Iend,  Iend-Istart,  Comm, Comm_Global"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)write(*,*)MpiRank,MpiQ,MpiR,mpiQdw,MpiRdw,MpiIstart,MpiIend,MpiIend-MpiIstart+1,MpiComm,MpiComm_Global
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
#endif
    !
    !
    !#################################
    !          HxV SETUP
    !#################################
    if(present(Hmat))then
       spHtimesV_p => null()
       call ed_buildh_main(Hmat)          
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       spHtimesV_p => spMatVec_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => spMatVec_MPI_main
#endif
       call ed_buildh_main()
    case (.false.)
       spHtimesV_p => directMatVec_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => directMatVec_MPI_main
#endif
    end select
    !
  end subroutine build_Hv_sector_main



  subroutine delete_Hv_sector_main()
    integer :: iud,ierr,i
    call delete_sector(Hsector)
    if(allocated(DimUps))deallocate(DimUps)
    if(allocated(DimDws))deallocate(DimDws)
    Dim    = 0
    DimUp  = 0
    DimDw  = 0
    !
    !There is no difference here between Mpi and serial version, as Local part was removed.
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0d)
       if(Jhflag)call sp_delete_matrix(MpiComm,spH0nd)
    else
       call sp_delete_matrix(spH0d)
       if(Jhflag)call sp_delete_matrix(spH0nd)
    endif
#else
    call sp_delete_matrix(spH0d)
    if(Jhflag)call sp_delete_matrix(spH0nd)
#endif
    do iud=1,Ns_Ud
       call sp_delete_matrix(spH0ups(iud))
       call sp_delete_matrix(spH0dws(iud))
    enddo
    !
    spHtimesV_p => null()
    !
#ifdef _MPI
    if(MpiStatus)then
       if(MpiGroup/=Mpi_Group_Null)call Mpi_Group_free(MpiGroup,ierr)
       if(MpiComm/=Mpi_Comm_Null.AND.MpiComm/=Mpi_Comm_World)call Mpi_Comm_Free(MpiComm,ierr)
       MpiComm = MpiComm_Global
       MpiSize = get_Size_MPI(MpiComm_Global)
       MpiRank = get_Rank_MPI(MpiComm_Global)
       mpiQup=0
       mpiRup=0
       mpiQdw=0
       mpiRdw=0
       mpiQ=0
       mpiR=0
       mpiIstart=0
       mpiIend=0
       mpiIshift=0
    endif
#endif
    !
  end subroutine delete_Hv_sector_main


  function vecDim_Hv_sector_main(isector) result(vecDim)
    integer :: isector
    integer :: vecDim
    integer :: mpiQdw
    integer :: DimUps(Ns_Ud),DimUp
    integer :: DimDws(Ns_Ud),DimDw
    !
    call get_DimUp(isector,DimUps) ; DimUp = product(DimUps)
    call get_DimDw(isector,DimDws) ; DimDw = product(DimDws)
    !
#ifdef _MPI
    if(MpiStatus)then
       !Dw split:
       mpiQdw = DimDw/MpiSize
       if(MpiRank < mod(DimDw,MpiSize) ) MpiQdw = MpiQdw+1
    else
       mpiQdw = DimDw
    endif
#else
    mpiQdw = DimDw
#endif
    !
    vecDim=DimUp*mpiQdw
    !
  end function vecDim_Hv_sector_main


  subroutine tridiag_Hv_sector_main(isector,vvinit,alanc,blanc,norm2)
    integer                             :: isector
    complex(8),dimension(:)             :: vvinit
    real(8),dimension(:),allocatable    :: alanc,blanc
    real(8)                             :: norm2
    complex(8),dimension(:),allocatable :: vvloc
    integer                             :: vecDim
    !
    if(MpiMaster)then
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
    endif
#ifdef _MPI
    if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
    call build_Hv_sector_main(isector)
    allocate(alanc(Hsector%Nlanc),blanc(Hsector%Nlanc))
    alanc=0d0 ; blanc=0d0
    if(norm2/=0d0)then
#ifdef _MPI
       if(MpiStatus)then
          vecDim = vecDim_Hv_sector_main(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alanc,blanc)
       else
          call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
#endif
    endif
    call delete_Hv_sector_main()
  end subroutine tridiag_Hv_sector_main








  !####################################################################
  !####################################################################
  !                          KONDO
  !####################################################################
  !####################################################################

  subroutine build_Hv_sector_kondo(isector,Hmat)
    integer                            :: isector
    complex(8),dimension(:,:),optional :: Hmat
    integer                            :: irank
    integer                            :: i,j
    !
    !
    call build_sector(isector,Hsector)
    !
    Dim    = Hsector%Dim
    !
    !#################################
    !          MPI SETUP
    !#################################
    mpiAllThreads=.true.
    MpiQ = Dim/MpiSize
    MpiR = 0
    if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    !
    MpiIshift = MpiRank*mpiQ
    MpiIstart = MpiRank*mpiQ + 1
    MpiIend   = (MpiRank+1)*mpiQ + mpiR
    !
#ifdef _MPI
#ifdef _DEBUG
    if(MpiStatus.AND.ed_verbose>4)then
       write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)then
             write(LOGfile,*)MpiRank,MpiQ,MpiR,MpiIstart,MpiIend,MpiIend-MpiIstart+1
          endif
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
#endif
    !
    !
    !#################################
    !          HxV SETUP
    !#################################
    if(present(Hmat))then
       spHtimesV_p => null()
       call ed_buildh_kondo(Hmat)          
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       spHtimesV_p => spMatVec_kondo
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => spMatVec_MPI_kondo
#endif
       call ed_buildh_kondo()
    case (.false.)
       spHtimesV_p => directMatVec_kondo
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => directMatVec_MPI_kondo
#endif
    end select
    !
  end subroutine build_Hv_sector_kondo



  subroutine delete_Hv_sector_kondo()
    !
    call delete_sector(Hsector)
    Dim    = 0
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0d)
    else
       call sp_delete_matrix(spH0d)
    endif
#else
    call sp_delete_matrix(spH0d)
#endif
    !
    spHtimesV_p => null()
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiComm = MpiComm_Global
       MpiSize = get_Size_MPI(MpiComm_Global)
       MpiRank = get_Rank_MPI(MpiComm_Global)
       mpiQ=0
       mpiR=0
       mpiIstart=0
       mpiIend=0
       mpiIshift=0
    endif
#endif
    !
  end subroutine delete_Hv_sector_kondo


  function vecDim_Hv_sector_kondo(isector) result(vecDim)
    integer :: isector
    integer :: vecDim
    integer :: Dim
    !
    Dim  = getdim(isector)
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiQ = Dim/MpiSize
       MpiR = 0
       if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    else
       MpiQ = Dim
       MpiR = 0
    endif
#else
    MpiQ = Dim
    MpiR = 0
#endif
    !
    vecDim=MpiQ + MpiR
    !
  end function vecDim_Hv_sector_kondo

  subroutine tridiag_Hv_sector_kondo(isector,vvinit,alanc,blanc,norm2)
    integer                             :: isector
    complex(8),dimension(:)             :: vvinit
    real(8),dimension(:),allocatable    :: alanc,blanc
    real(8)                             :: norm2
    complex(8),dimension(:),allocatable :: vvloc
    integer                             :: vecDim
    !
    !
    if(MpiMaster)then
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
    endif
#ifdef _MPI
    if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
    call build_Hv_sector_kondo(isector)
    allocate(alanc(Hsector%Nlanc),blanc(Hsector%Nlanc))
    alanc=0d0 ; blanc=0d0
    if(norm2/=0d0)then
#ifdef _MPI
       if(MpiStatus)then
          vecDim = vecDim_Hv_sector_kondo(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alanc,blanc)
       else
          call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
#endif
    endif
    call delete_Hv_sector_kondo()
  end subroutine tridiag_Hv_sector_kondo



end MODULE ED_HAMILTONIAN
