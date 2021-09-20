module ED_GRAPH_MATRIX
  USE SF_IOTOOLS
  USE SF_TIMER,     only: start_timer,stop_timer
  USE SF_CONSTANTS, only: zero,one
  USE SF_LINALG,    only: eye,inv,diag
  USE SF_MISC,      only: assert_shape
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


  type link
     integer            :: siteI,siteJ
     integer            :: orbI,orbJ
     integer            :: spin
     complex(8)         :: Hval
     type(link),pointer :: next !link to next box (chain)
  end type link

  type Hij_structure
     type(link),pointer                      :: root !head/root of the list\== list itself
     character(len=:),allocatable            :: file     !Name of the W90 file 
     integer                                 :: Size=0   !Number of hopping elements
     integer,allocatable,dimension(:)        :: Nsites    !Number of sites per orbital and spin
     integer                                 :: Ns=0     !total Number of lattice sites per spin
     integer                                 :: Norb=0   !Number of orbitals
     integer                                 :: Nspin=0  !Number of spin channels (1 or 2)
     complex(8),allocatable,dimension(:,:,:) :: Hij    !The H(Ri,Rj)_ab^sr Hamiltonian [Ns,Ns]
     complex(8),allocatable,dimension(:,:,:) :: Hloc   !The local part of H(Ri,Ri)_ab^sr
     logical                                 :: built  =.false. !Hij is built
     logical                                 :: status =.false. !Allocated
  end type hij_structure

  interface Hij_get
     module procedure :: Hij_get_main
     module procedure :: Hij_get_ispin
     module procedure :: Hij_get_whole
  end interface Hij_get

  interface Hij_local
     module procedure :: Hij_local_main
     module procedure :: Hij_local_ispin
     module procedure :: Hij_local_whole
  end interface Hij_local


  public :: Hij_structure
  !
  public :: Hij_init
  public :: Hij_delete
  public :: Hij_info
  public :: Hij_add_link
  public :: Hij_read
  public :: Hij_build
  public :: Hij_get
  public :: Hij_local
  public :: Hij_write
  public :: Hij_get_g0inv
  public :: Hij_get_g0func


  integer :: mpi_ierr
  integer :: mpi_rank
  integer :: mpi_size
  logical :: mpi_master


  type(Hij_structure),public :: Hmatrix


contains



  !< Setup/Delete the default structure
  subroutine Hij_init(Nsites,Nspin)
    integer,dimension(:)      :: Nsites
    integer,optional          :: Nspin
    integer                   :: Nspin_
    integer                   :: Ns
    integer                   :: Norb
    !
    Nspin_ = 1; if(present(Nspin))Nspin_=Nspin
    !
    if(Hmatrix%status)call Hij_delete()
    !
    Norb = size(Nsites)
    Ns   = sum(Nsites)
    !
    allocate(Hmatrix%root)
    Hmatrix%root%next => null()
    !
    allocate(Hmatrix%Nsites(Norb))
    Hmatrix%Nsites  = Nsites
    Hmatrix%file   = ""
    Hmatrix%size   = 0
    Hmatrix%Ns     = Ns
    Hmatrix%Norb   = Norb
    Hmatrix%Nspin  = Nspin_
    !
    allocate( Hmatrix%Hij(Nspin_,Ns,Ns) )
    Hmatrix%Hij    = zero
    allocate( Hmatrix%Hloc(Nspin_,Ns,Ns) )
    Hmatrix%Hloc   = zero
    Hmatrix%built  =.false.
    Hmatrix%status =.true.
    return
  end subroutine Hij_init
  !
  subroutine Hij_delete()
    type(link),pointer  :: p,c
    !
    if(.not.Hmatrix%status)return
    !
    do
       p => Hmatrix%root
       c => p%next    !current is the first node (root's next)
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next => null()
       call del_link(c)
       deallocate(c)
    end do
    deallocate(Hmatrix%root)
    deallocate(Hmatrix%Nsites)
    deallocate(Hmatrix%file)
    deallocate(Hmatrix%Hij)
    deallocate(Hmatrix%Hloc)
    Hmatrix%size   = 0
    Hmatrix%Ns     = 0
    Hmatrix%Norb   = 0
    Hmatrix%Nspin  = 0
    Hmatrix%built  =.false.
    Hmatrix%status =.false.
  end subroutine Hij_delete

  subroutine Hij_info()
    !
    if(.not.Hmatrix%status)then
       write(*,"(A)")"Hij_info WARNING: Hmatrix not allocated"
       return
    endif
    !
    write(*,"(A,10I8)")"Nsites  =",Hmatrix%Nsites
    write(*,"(A,I8)")"Ns      =",Hmatrix%Ns
    write(*,"(A,I8)")"Nspin   =",Hmatrix%Nspin
    write(*,"(A,I8)")"Norb    =",Hmatrix%Norb
    write(*,"(A,I8)")"# link  =",Hmatrix%size
    if(Hmatrix%file /= "")write(*,"(A,1x,A)")"file    =",Hmatrix%file
    write(*,"(A,L8)")"built H =",Hmatrix%built
    return
  end subroutine Hij_info





  !Add/Remove link to a given structure
  subroutine Hij_add_link(siteI,siteJ,orbI,orbJ,spin,Hval)
    integer                           :: siteI,siteJ
    integer                           :: orbI,orbJ
    integer                           :: spin
    complex(8) ,intent(in)            :: Hval
    type(link),pointer                :: p,c
    integer                           :: k
    !
    !
    if(.not.Hmatrix%status)stop "Hij_add_link ERROR: Hmatrix not allocated"
    !
    !
    if(spin>Hmatrix%Nspin .OR. spin<=0)&
         stop "Hij_add_link ERROR: spinI or spinJ either > Nspin OR < 0"
    if(orbI>Hmatrix%Norb  .OR. orbJ>Hmatrix%Norb  .OR. orbI<=0 .OR. orbJ<=0)&
         stop "Hij_add_link ERROR: orbI or orbJ either > Norb OR < 0"
    if(siteI>Hmatrix%Nsites(orbI) .OR. siteJ>Hmatrix%Nsites(orbJ) .OR. siteI<=0 .OR. siteI<=0 )&
         stop "Hij_add_link ERROR: siteI or siteJ either > Nlat OR < 0"
    !
    p => Hmatrix%root
    c => p%next
    do k=1,Hmatrix%size                    !traverse the Hmatrix
       if(.not.associated(c))exit !beginning of the Hmatrix
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the Hmatrix
    !
    p%next%siteI = siteI
    p%next%siteJ = siteJ
    p%next%orbI  = orbI
    p%next%orbJ  = orbJ
    p%next%spin  = spin
    p%next%Hval  = Hval
    Hmatrix%size=Hmatrix%size+1
    !
    if(.not.associated(c))then !end of the Hmatrix special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine Hij_add_link




  !Set/Get/Delete link: INTERNAL USE
  subroutine set_link(self,siteI,siteJ,orbI,orbJ,spin,Hval)
    type(link),intent(inout) :: self
    integer,intent(in)       :: siteI,siteJ
    integer,intent(in)       :: orbI,orbJ
    integer,intent(in)       :: spin
    complex(8),intent(in)    :: Hval
    self%siteI = siteI
    self%siteJ = siteJ
    self%orbI  = orbI
    self%orbJ  = orbJ
    self%spin  = spin
    self%Hval  = Hval
  end subroutine set_link
  !

  subroutine get_link(self,siteI,siteJ,orbI,orbJ,spin,Hval)
    type(link),intent(inout) :: self
    integer,intent(inout)    :: siteI,siteJ
    integer,intent(inout)    :: orbI,orbJ
    integer,intent(inout)    :: spin
    complex(8),intent(inout) :: Hval
    siteI = self%siteI
    siteJ = self%siteJ
    orbI  = self%orbI
    orbJ  = self%orbJ
    spin  = self%spin
    Hval  = self%Hval
  end subroutine get_link
  !
  subroutine del_link(self)
    type(link),intent(inout) :: self
    self%siteI = 0
    self%siteJ = 0
    self%orbI  = 0
    self%orbJ  = 0
    self%spin  = 0
    self%Hval  = zero
    self%next  => null()
  end subroutine del_link
  !
  subroutine print_link(self,unit)
    type(link),intent(in) :: self
    integer               :: unit
    write(unit,"(5I6,2F12.4)")&
         self%siteI,self%siteJ,&
         self%orbI,self%orbJ,&
         self%spin,dreal(self%Hval),dimag(self%Hval)
  end subroutine print_link
  !
  function indices_link(siteI,orbI,spin) result(Indx)
    integer,intent(in)          :: siteI,orbI
    integer,intent(in),optional :: spin
    integer                     :: Indx
    select case(orbI)
    case(1)
       Indx = siteI
    case default
       Indx = siteI + (orbI-1)*Hmatrix%Nsites(orbI-1)
    end select
    if(present(spin))Indx = Indx + (spin-1)*Hmatrix%Ns
  end function indices_link









  !< Read a given structure from file
  subroutine Hij_read(file)
    character(len=*)          :: file
    integer                   :: Nsize
    integer                   :: unitIO
    integer                   :: ih,i,j,a,b,spin
    real(8)                   :: re,im
    !
    if(.not.Hmatrix%status)stop "Hij_read ERROR: Hmatrix not allocated"
    !
    !Read from file initial info OR reconstruct them is Header is missing
    open(free_unit(unitIO),file=reg(file),status="old",action="read")
    Nsize = file_length(reg(file))
    Hmatrix%file = reg(file)
    Hmatrix%Size = Nsize
    !
    do ih=1,Hmatrix%Size
       read(unitIO,*)i,j,a,b,spin,re,im
       call Hij_add_link(i,j,a,b,spin,dcmplx(re,im))
    enddo
    close(unitIO)
    return
  end subroutine Hij_read



  !< Build Hij from a given structure
  subroutine Hij_build()
    type(link),pointer :: c           
    integer            :: i,j,a,b,spin
    integer            :: io,jo
    complex(8)         :: val
    !
    if(.not.Hmatrix%status)stop "Hij_build ERROR: Hmatrix not allocated"
    !
    c => Hmatrix%root%next
    do
       if(.not.associated(c))exit
       call get_link(c,i,j,a,b,spin,val)
       io = indices_link(i,a)!,spin)
       jo = indices_link(j,b)!,spin)
       Hmatrix%Hij(spin,io,jo) = val
       if(i==j)Hmatrix%Hloc(spin,io,jo) = val
       c => c%next
    enddo
    do spin=1,Hmatrix%Nspin
       if(.not.herm_check(Hmatrix%Hij(spin,:,:),1d-6))stop "Hij_build ERROR: Hij is not Hermitian"
    enddo
    Hmatrix%built=.true.
    return
  end subroutine Hij_build




  subroutine Hij_get_main(Hij)
    complex(8),dimension(:,:,:) :: Hij
    integer                     :: Nlso,Nspin
    if(.not.Hmatrix%status)stop "Hij_get ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Ns
    Nspin= Hmatrix%Nspin
    call assert_shape(Hij,[Nspin,Nlso,Nlso],"Hij_get","Hij")
    if(.not.Hmatrix%built)call Hij_build()
    Hij = Hmatrix%Hij
  end subroutine Hij_Get_Main
  !
  subroutine Hij_get_ispin(Hij,spin)
    complex(8),dimension(:,:) :: Hij
    integer                   :: spin
    integer                   :: Nlso
    if(.not.Hmatrix%status)stop "Hij_get ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Ns
    call assert_shape(Hij,[Nlso,Nlso],"Hij_get","Hij")
    if(.not.Hmatrix%built)call Hij_build()
    Hij = Hmatrix%Hij(spin,:,:)
  end subroutine Hij_Get_Ispin
  !
  subroutine Hij_get_whole(Hij)
    complex(8),dimension(:,:) :: Hij
    integer                   :: spin
    integer                   :: Nlso
    if(.not.Hmatrix%status)stop "Hij_get ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Nspin*Hmatrix%Ns
    call assert_shape(Hij,[Nlso,Nlso],"Hij_get","Hij")
    if(.not.Hmatrix%built)call Hij_build()
    Hij = pack_matrix(Hmatrix%Hij,Hmatrix%Ns,Hmatrix%Nspin)
  end subroutine Hij_Get_Whole




  subroutine Hij_local_main(Hloc)
    complex(8),dimension(:,:,:) :: Hloc
    integer                     :: Nlso,Nspin
    if(.not.Hmatrix%status)stop "Hij_local ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Ns
    Nspin= Hmatrix%Nspin
    call assert_shape(Hloc,[Nspin,Nlso,Nlso],"Hij_local","Hloc")
    if(.not.Hmatrix%built)call Hij_build()
    Hloc = Hmatrix%Hloc
  end subroutine Hij_local_main
  !
  subroutine Hij_local_ispin(Hloc,spin)
    complex(8),dimension(:,:) :: Hloc
    integer                   :: spin
    integer                   :: Nlso
    if(.not.Hmatrix%status)stop "Hij_local ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Ns
    call assert_shape(Hloc,[Nlso,Nlso],"Hij_local","Hloc")
    if(.not.Hmatrix%built)call Hij_build()
    Hloc = Hmatrix%Hloc(spin,:,:)
  end subroutine Hij_Local_Ispin
  !
  subroutine Hij_local_whole(Hloc)
    complex(8),dimension(:,:) :: Hloc
    integer                   :: spin
    integer                   :: Nlso
    if(.not.Hmatrix%status)stop "Hij_local ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Nspin*Hmatrix%Ns
    call assert_shape(Hloc,[Nlso,Nlso],"Hij_get","Hij")
    if(.not.Hmatrix%built)call Hij_build()
    Hloc = pack_matrix(Hmatrix%Hloc,Hmatrix%Ns,Hmatrix%Nspin)
  end subroutine Hij_Local_Whole






  subroutine Hij_write(unit,file)
    character(len=*),optional :: file
    integer,optional          :: unit
    integer                   :: unit_
    type(link),pointer        :: c           
    logical                   :: mpi_master
    !
    if(.not.Hmatrix%status)stop "Hij_write ERROR: Hmatrix not allocated"
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(.not.Hmatrix%status)stop "write_Hmatrix ERROR: Hmatrix structure not allocated. Call setup_Hmatrix first."
    !
    if(mpi_master)then
       unit_ = 6 ; if(present(unit))unit_=unit
       if(present(file))open(free_unit(unit_),file=reg(file))
       c => Hmatrix%root%next
       do
          if(.not.associated(c))exit
          call print_link(c,unit_)
          c => c%next
       enddo
       if(present(file))close(unit_)
    endif
    return
  end subroutine Hij_write








  subroutine Hij_get_g0inv(zeta,Gloc)
    complex(8),dimension(:)                     :: zeta
    complex(8),dimension(:,:,:,:),intent(inout) :: Gloc      ![Nspin][Ns][Ns][Lfreq]
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable     :: Hij       ![Nspin][Ns][Ns]
    integer                                     :: Lfreq
    integer                                     :: Ns,Nspin
    integer                                     :: i,is,ispin
    !
    !
    if(.not.Hmatrix%status)stop "Hij_get_gfunc ERROR: Hmatrix not allocated"
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    Lfreq=size(zeta)
    !
    Ns    = Hmatrix%Ns
    Nspin = Hmatrix%Nspin
    Lfreq = size(zeta)
    !Testing part:
    call assert_shape(Gloc,[Nspin,Ns,Ns,Lfreq],'dmft_get_gloc_normal_gij',"Gloc")
    !
    allocate(Hij(Nspin,Ns,Ns))
    call Hij_get(Hij)
    !
    !pass each Z_site to the routines that invert (Z-Hij) for each k-point 
    do ispin=1,Nspin
       do i=1,Lfreq
          Gloc(ispin,:,:,i) = zeta(i)*eye(Ns) - Hij(ispin,:,:)
       enddo
    enddo
  end subroutine Hij_get_g0inv



  subroutine Hij_get_g0func(zeta,Gloc)
    complex(8),dimension(:)                     :: zeta
    complex(8),dimension(:,:,:,:),intent(inout) :: Gloc      ![Nspin][Ns][Ns][Lfreq]
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable     :: Hij       ![Nspin][Ns][Ns]
    complex(8),dimension(:,:),allocatable       :: csi ![Ns][Lfreq]
    integer                                     :: Lfreq
    integer                                     :: Ns,Nspin
    integer                                     :: i,is,ispin
    !
    !
    if(.not.Hmatrix%status)stop "Hij_get_gfunc ERROR: Hmatrix not allocated"
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    Lfreq=size(zeta)
    !
    Ns    = Hmatrix%Ns
    Nspin = Hmatrix%Nspin
    Lfreq = size(zeta)
    !Testing part:
    call assert_shape(Gloc,[Nspin,Ns,Ns,Lfreq],'dmft_get_gloc_normal_gij',"Gloc")
    !
    allocate(Hij(Nspin,Ns,Ns))
    call Hij_get(Hij)
    !
    allocate(csi(Ns,Lfreq))
    forall(is=1:Ns)csi(is,:) = zeta(:)
    !
    !pass each Z_site to the routines that invert (Z-Hij) for each k-point 
    if(mpi_master)call start_timer
    do ispin=1,Nspin
       call invert_gij_mpi(csi,Hij(ispin,:,:),Gloc(ispin,:,:,:))
    enddo
    if(mpi_master)call stop_timer
  end subroutine Hij_get_g0func







  !##################################################################
  !                   COMPUTATIOAL ROUTINES
  !##################################################################
  !PARALLEL ON FREQ:
  subroutine invert_gij_mpi(zeta,Hij,Gout)
    complex(8),dimension(:,:),intent(in)      :: zeta    ![Ns][Lfreq]
    complex(8),dimension(:,:),intent(in)      :: Hij     ![Ns][Ns]
    complex(8),dimension(:,:,:),intent(inout) :: Gout    ![Ns][Ns][Lfreq]
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable   :: Gtmp
    complex(8),dimension(:,:),allocatable     :: Gmatrix
    integer                                   :: Lfreq
    integer                                   :: Ns
    integer                                   :: i
    !
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    Ns    = size(zeta,1)
    Lfreq = size(Gout,3)
    call assert_shape(zeta,[Ns,Lfreq],"invert_gij_mpi","zeta")
    call assert_shape(Hij,[Ns,Ns],"invert_gij_mpi","Hij")
    call assert_shape(Gout,[Ns,Ns,Lfreq],"invert_gij_mpi","Gout")
    !
    allocate(Gtmp(Ns,Ns,Lfreq))
    allocate(Gmatrix(Ns,Ns))
    Gtmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = diag(zeta(:,i)) - Hij
       call inv(Gmatrix)
       Gtmp(:,:,i) = Gmatrix
    enddo
    !
#ifdef _MPI    
    if(check_MPI())then
       Gout=zero
       call MPI_AllReduce(Gtmp, Gout, size(Gout), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Gout=Gtmp
    endif
#else
    Gout=Gtmp
#endif
    return
  end subroutine invert_gij_mpi




  function herm_check(A,err) result(bool)
    complex(8),intent(in) :: A(:,:)
    real(8),optional      :: err
    real(8)               :: err_
    logical               :: bool
    err_ = 0d0;if(present(err))err_=err
    bool = .true.
    if( any( abs(A-conjg(transpose(A))) > err_ ) )bool=.false.
  end function herm_check




  function pack_matrix(Hin,Ns,Nspin) result(Hout)
    complex(8),dimension(Nspin,Ns,Ns)       :: Hin
    complex(8),dimension(Nspin*Ns,Nspin*Ns) :: Hout
    integer                                 :: Ns,Nspin
    integer                                 :: ispin,is,js
    integer                                 :: io,jo
    do ispin=1,Nspin
       do is=1,Ns
          do js=1,Ns
             io = is + (ispin-1)*Ns
             jo = js + (ispin-1)*Ns
             Hout(io,jo) = Hin(ispin,is,js)
          enddo
       enddo
    enddo
  end function pack_matrix
  !
  function unpack_matrix(Hin,Ns,Nspin) result(Hout)
    complex(8),dimension(Nspin*Ns,Nspin*Ns) :: Hin
    complex(8),dimension(Nspin,Ns,Ns)       :: Hout
    integer                                 :: Ns,Nspin
    integer                                 :: ispin,is,js
    integer                                 :: io,jo
    do ispin=1,Nspin
       do is=1,Ns
          do js=1,Ns
             io = is + (ispin-1)*Ns
             jo = js + (ispin-1)*Ns
             Hout(ispin,is,js) = Hin(io,jo)
          enddo
       enddo
    enddo
  end function unpack_matrix






END MODULE ED_GRAPH_MATRIX
















