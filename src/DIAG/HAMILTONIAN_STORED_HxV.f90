! > BUILD SPARSE HAMILTONIAN of the SECTOR
MODULE ED_HAMILTONIAN_STORED_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Matric constructors
  public :: ed_buildh_main

  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_main
#ifdef _MPI
  public  :: spMatVec_MPI_main
#endif


contains


  subroutine ed_buildh_main(Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    integer                               :: isector,ispin,i,j
    complex(8),dimension(:,:),allocatable :: Htmp_up,Htmp_dw,Hrdx,Hmat_tmp
    integer,dimension(2*Ns_Ud)            :: Indices    ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)       :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                 :: Nup,Ndw
    real(8),dimension(Ns)                 :: Sz 
    complex(8),dimension(Nspin,Ns,Ns)     :: Hij,Hloc
    real(8),dimension(Nspin,Ns)           :: Hdiag
    !
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = dreal(diagonal(Hloc(ispin,:,:)))
    enddo
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,DimUp*DimDw)
       if(Jhflag)then
          call sp_set_mpi_matrix(MpiComm,spH0nd,mpiIstart,mpiIend,mpiIshift)
          call sp_init_matrix(MpiComm,spH0nd,DimUp*DimDw)
       endif
    else
       call sp_init_matrix(spH0d,DimUp*DimDw)
       if(Jhflag)call sp_init_matrix(spH0nd,DimUp*DimDw)
    endif
#else
    call sp_init_matrix(spH0d,DimUp*DimDw)
    if(Jhflag)call sp_init_matrix(spH0nd,DimUp*DimDw)
#endif
    call sp_init_matrix(spH0dws(1),DimDw)
    call sp_init_matrix(spH0ups(1),DimUp)
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
    include "stored/H_local.f90"
    if(ed_filling==0)then
       include "stored/H_muhf.f90"
    endif
    !
    !NON-LOCAL HAMILTONIAN TERMS
    if(jhflag)then
       include "stored/H_non_local.f90"
    endif
    !
    !UP TERMS
    include "stored/H_up.f90"
    !
    !DW TERMS
    include "stored/H_dw.f90"
    !
    !-----------------------------------------------!
    if(present(Hmat))then
       Hmat = zero
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=zero
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=zero
       allocate(Hmat_tmp(DimUp*DimDw,DimUp*DimDw));Hmat_tmp=zero
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat_tmp)
       else
          call sp_dump_matrix(spH0d,Hmat_tmp)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat_tmp)
#endif
       !
       if(Jhflag)then
          allocate(Hrdx(DimUp*DimDw,DimUp*DimDw));Hrdx=zero
#ifdef _MPI
          if(MpiStatus)then
             call sp_dump_matrix(MpiComm,spH0nd,Hrdx)
          else
             call sp_dump_matrix(spH0nd,Hrdx)
          endif
#else
          call sp_dump_matrix(spH0nd,Hrdx)
#endif
          Hmat_tmp = Hmat_tmp + Hrdx
          deallocate(Hrdx)
       endif
       !
       call sp_dump_matrix(spH0ups(1),Htmp_up)
       call sp_dump_matrix(spH0dws(1),Htmp_dw)
       Hmat_tmp = Hmat_tmp + kronecker_product(Htmp_dw,zeye(DimUp))
       Hmat_tmp = Hmat_tmp + kronecker_product(zeye(DimDw),Htmp_up)
       !
       Hmat = Hmat_tmp
       !
       deallocate(Htmp_up,Htmp_dw,Hmat_tmp)
    endif
    !
    return
    !
  end subroutine ed_buildh_main

















  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial
  ! - MPI
  !+------------------------------------------------------------------+
  subroutine spMatVec_main(Nloc,v,Hv)
    integer                    :: Nloc
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    complex(8)                 :: val
    integer                    :: i,iup,idw,j,jup,jdw,jj
    !
    !
    Hv=zero
    !
    !Local:
    do i = 1,Nloc
       do jj=1,spH0d%row(i)%Size
          val = spH0d%row(i)%vals(jj)
          j = spH0d%row(i)%cols(jj)
          Hv(i) = Hv(i) + val*v(j)
       enddo
    enddo
    !
    !DW:
    do iup=1,DimUp
       do idw=1,DimDw
          i = iup + (idw-1)*DimUp
          do jj=1,spH0dws(1)%row(idw)%Size
             jup = iup
             jdw = spH0dws(1)%row(idw)%cols(jj)
             val = spH0dws(1)%row(idw)%vals(jj)
             j     = jup +  (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !UP:
    do idw=1,DimDw
       !
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j =  jup + (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !

    !Non-Local:
    if(jhflag)then
       do i = 1,Nloc         
          do jj=1,spH0nd%row(i)%Size
             val = spH0nd%row(i)%vals(jj)
             j = spH0nd%row(i)%cols(jj)
             Hv(i) = Hv(i) + val*v(j)
          enddo
       enddo
    endif
    !
  end subroutine spMatVec_main



#ifdef _MPI
  subroutine spMatVec_mpi_main(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    !
    integer                             :: N
    complex(8),dimension(:),allocatable :: vt,Hvt
    complex(8),dimension(:),allocatable :: vin
    complex(8)                          :: val
    integer                             :: i,iup,idw,j,jup,jdw,jj
    integer                             :: irank
    !
    ! if(MpiComm==Mpi_Comm_Null)return
    ! if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=zero
    do i=1,Nloc                 !==spH0%Nrow
       do jj=1,spH0d%row(i)%Size
          val = spH0d%row(i)%vals(jj)
          Hv(i) = Hv(i) + val*v(i)
       enddo
    end do
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do idw=1,MpiQdw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          hxv_up: do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j   = jup + (idw-1)*DimUp
             Hv(i) = Hv(i) + val*v(j)
          end do hxv_up
       enddo
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    allocate(vt(mpiQup*DimDw))
    allocate(Hvt(mpiQup*DimDw))
    vt=0d0
    Hvt=0d0
    call vector_transpose_MPI(DimUp,MpiQdw,v(1:DimUp*MpiQdw),DimDw,MpiQup,vt)
    do idw=1,MpiQup             !<= Transposed order:  column-wise DW <--> UP  
       do iup=1,DimDw           !<= Transposed order:  column-wise DW <--> UP
          i = iup + (idw-1)*DimDw
          hxv_dw: do jj=1,spH0dws(1)%row(iup)%Size
             jup = spH0dws(1)%row(iup)%cols(jj)
             jdw = idw             
             j   = jup + (jdw-1)*DimDw
             val = spH0dws(1)%row(iup)%vals(jj)
             Hvt(i) = Hvt(i) + val*vt(j)
          end do hxv_dw
       enddo
    end do
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=0d0
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
    Hv(1:DimUp*MpiQdw) = Hv(1:DimUp*MpiQdw) + Vt
    deallocate(vt,Hvt)
    !
    !
    !Non-Local:
    if(jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       ! 
       allocate(vt(N)) ; vt = 0d0
       call allgather_vector_MPI(MpiComm,v,vt)
       !
       do i=1,Nloc
          matmul: do jj=1,spH0nd%row(i)%Size
             val = spH0nd%row(i)%vals(jj)
             j = spH0nd%row(i)%cols(jj)
             Hv(i) = Hv(i) + val*Vt(j)
          enddo matmul
       enddo
       deallocate(Vt)
    endif
    !
  end subroutine spMatVec_mpi_main
#endif



end MODULE ED_HAMILTONIAN_STORED_HXV







