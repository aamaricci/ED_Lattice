MODULE ED_SECTOR
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private




  public :: build_sector
  public :: delete_sector
  public :: show_sector
  !
  public :: apply_op_C
  public :: apply_op_CDG
  public :: apply_op_Sz
  public :: apply_op_N
  public :: build_op_Ns

  public :: get_Sector
  public :: get_QuantumNumbers
  public :: get_Nup
  public :: get_Ndw
  public :: get_DimUp
  public :: get_DimDw
  !
  public :: indices2state
  public :: state2indices
  public :: iup_index
  public :: idw_index
  !
  public :: twin_sector_order
  public :: get_twin_sector
  public :: flip_state


  interface apply_op_Sz
     module procedure :: apply_op_Sz_single
     module procedure :: apply_op_Sz_range
  end interface apply_op_Sz

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate



contains


  !##################################################################
  !##################################################################
  !BUILD SECTORS
  !##################################################################
  !##################################################################
  subroutine build_sector(isector,self)
    integer,intent(in)                  :: isector
    type(sector)                        :: self
    integer                             :: i,iup,idw,ipup,ipdw
    integer                             :: nup_,ndw_
    integer                             :: dim,iel,iimp
    !
    if(self%status)call delete_sector(self)
    !
    self%index = isector
    !
    if(KondoFlag)then
       !
       allocate(self%H(1))
       allocate(self%Nups(1))
       allocate(self%Ndws(1))
       call get_Nup(isector,self%Nups);self%Nup=self%Nups(1)
       call get_Ndw(isector,self%Ndws);self%Ndw=self%Ndws(1)
       self%DimEl=getDim(isector)
       self%Dim=self%DimEl
       !
       call map_allocate(self%H(1),self%Dim)
       dim = 0
       do ipdw=0,2**iNs-1
          do ipup=0,2**iNs-1
             if(popcnt(ipup)+popcnt(ipdw) /= Nimp)cycle
             do idw=0,2**eNs-1
                if(popcnt(idw) + popcnt(ipdw) /= self%Ndw) cycle
                do iup=0,2**eNs-1
                   if(popcnt(iup) + popcnt(ipup) /= self%Nup) cycle
                   dim        = dim+1
                   iel        = iup + idw*2**eNs
                   iimp       = ipup+ipdw*2**iNs
                   self%H(1)%map(dim) =  iel + iimp*2**(2*eNs)
                enddo
             enddo
          enddo
       enddo
       !
    else
       !
       allocate(self%H(2))
       allocate(self%DimUps(1))
       allocate(self%DimDws(1))
       allocate(self%Nups(1))
       allocate(self%Ndws(1))
       !
       call get_Nup(isector,self%Nups);self%Nup=sum(self%Nups)
       call get_Ndw(isector,self%Ndws);self%Ndw=sum(self%Ndws)
       call get_DimUp(isector,self%DimUps);self%DimUp=product(self%DimUps)
       call get_DimDw(isector,self%DimDws);self%DimDw=product(self%DimDws)
       self%DimEl=self%DimUp*self%DimDw
       self%Dim=self%DimEl
       !
       call map_allocate(self%H,[self%DimUps,self%DimDws])
       !UP    
       dim=0
       do iup=0,2**Ns-1
          nup_ = popcnt(iup)
          if(nup_ /= self%Nups(1))cycle
          dim  = dim+1
          self%H(1)%map(dim) = iup
       enddo
       !DW
       dim=0
       do idw=0,2**Ns-1
          ndw_= popcnt(idw)
          if(ndw_ /= self%Ndws(1))cycle
          dim = dim+1
          self%H(2)%map(dim) = idw
       enddo
       !
    endif
    !
    self%Nlanc = min(self%Dim,lanc_nGFiter)
    self%status=.true.
    !
  end subroutine build_sector


  subroutine delete_sector(self)
    type(sector) :: self
    call map_deallocate(self%H)
    if(allocated(self%H))deallocate(self%H)
    if(allocated(self%DimUps))deallocate(self%DimUps)
    if(allocated(self%DimDws))deallocate(self%DimDws)
    if(allocated(self%Nups))deallocate(self%Nups)
    if(allocated(self%Ndws))deallocate(self%Ndws)
    self%index=0
    self%DimUp=0
    self%DimDw=0
    self%Dim=0
    self%DimEl=0
    self%Nup=0
    self%Ndw=0
    self%Nlanc=0
    self%status=.false.
  end subroutine delete_sector



  subroutine show_sector(isector)
    integer,intent(in)                  :: isector
    integer                             :: i,iup,idw,ipup,ipdw
    integer                             :: nup_,ndw_,Nups(1),Ndws(1),Nup,Ndw
    integer                             :: dim,iel,iimp,DimEl,DimUps(1),DimDws(1),DimUp,DimDw
    !
    if(KondoFlag)then
       !
       call get_Nup(isector,Nups);Nup=Nups(1)
       call get_Ndw(isector,Ndws);Ndw=Ndws(1)
       DimEl=getDim(isector)
       !
       dim = 0
       do ipdw=0,2**iNs-1
          do ipup=0,2**iNs-1
             if(popcnt(ipup)+popcnt(ipdw) /= Nimp)cycle
             do idw=0,2**eNs-1
                if(popcnt(idw) + popcnt(ipdw) /= Ndw) cycle
                do iup=0,2**eNs-1
                   if(popcnt(iup) + popcnt(ipup) /= Nup) cycle
                   dim        = dim+1
                   iel        = iup + idw*2**eNs
                   iimp       = ipup+ipdw*2**iNs
                   write(*,"(4I6,A1,2I12,A1,I12)",advance='no')&
                        iup,idw,ipup,ipdw,"-",&
                        iup + idw*2**eNs,ipup+ipdw*2**iNs,"-",iel + iimp*2**(2*eNs)
                   call print_conf(iup,eNs,.false.)
                   call print_conf(idw,eNs,.false.)
                   call print_conf(ipup,iNs,.false.)
                   call print_conf(ipdw,iNs,.true.)
                enddo
             enddo
          enddo
       enddo
       !
    else
       !
       call get_Nup(isector,Nups);Nup=sum(Nups)
       call get_Ndw(isector,Ndws);Ndw=sum(Ndws)
       call get_DimUp(isector,DimUps);DimUp=product(DimUps)
       call get_DimDw(isector,DimDws);DimDw=product(DimDws)
       DimEl=DimUp*DimDw
       !
       !UP
       do idw=0,2**Ns-1
          if(popcnt(idw) /= ndw) cycle
          do iup=0,2**Ns-1
             if(popcnt(iup) /= nup) cycle
             write(*,"(2I6,A1,I12)",advance='no')&
                  iup,idw,"-",iup + idw*2**Ns
             call print_conf(iup,Ns,.false.)
             call print_conf(idw,Ns,.true.)
          enddo
       enddo
    endif
    !
  end subroutine show_sector





  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer                             :: i
    do i=1,size(H)
       call map_allocate_scalar(H(i),N(i))
    enddo
  end subroutine map_allocate_vector






  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector









  subroutine apply_op_C(i,j,sgn,ipos,ispin,sectorI,sectorJ) 
    integer, intent(in)     :: i,ipos,ispin
    type(sector),intent(in) :: sectorI,sectorJ
    integer,intent(out)     :: j
    real(8),intent(out)     :: sgn
    integer                 :: r
    integer                 :: i_el,ii,iorb
    integer,dimension(2)    :: Indices
    integer,dimension(2)    :: Jndices
    integer,dimension(eNs)  :: Nud
    integer                 :: Iud
    integer,dimension(2*Ns) :: ib
    !
    j=0
    sgn=0d0
    !
    if(KondoFlag)then
       i_el = sectorI%H(1)%map(i)
       ib   = bdecomp(i_el,2*Ns)
       if(ib(ipos)/=1)return
       call c(ipos,i_el,r,sgn)
       j    = binary_search(sectorJ%H(1)%map,r)
    else
       call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
       iud = sectorI%H(ispin)%map(Indices(ispin))
       nud = Bdecomp(iud,eNs)
       if(nud(ipos)/=1)return
       call c(ipos,iud,r,sgn)
       Jndices        = Indices
       Jndices(ispin) = binary_search(sectorJ%H(ispin)%map,r)
       call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
    endif
  end subroutine apply_op_C


  subroutine apply_op_CDG(i,j,sgn,ipos,ispin,sectorI,sectorJ) 
    integer, intent(in)     :: i,ipos,ispin
    type(sector),intent(in) :: sectorI,sectorJ
    integer,intent(out)     :: j
    real(8),intent(out)     :: sgn
    integer                 :: r
    integer                 :: i_el,ii,iorb
    integer,dimension(2)    :: Indices
    integer,dimension(2)    :: Jndices
    integer,dimension(eNs)  :: Nud
    integer                 :: Iud
    integer,dimension(2*Ns) :: ib
    !
    j=0
    sgn=0d0
    !
    if(KondoFlag)then
       i_el = sectorI%H(1)%map(i)
       ib   = bdecomp(i_el,2*Ns)
       if(ib(ipos)/=0)return
       call cdg(ipos,i_el,r,sgn)
       j    = binary_search(sectorJ%H(1)%map,r)
    else
       call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
       iud = sectorI%H(ispin)%map(Indices(ispin))
       nud = Bdecomp(iud,eNs)
       if(nud(ipos)/=0)return
       call cdg(ipos,iud,r,sgn)
       Jndices        = Indices
       Jndices(ispin) = binary_search(sectorJ%H(ispin)%map,r)
       call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
    endif
  end subroutine apply_op_CDG


  subroutine apply_op_Sz_single(i,sgn,ipos,sectorI) 
    integer, intent(in)     :: i,ipos
    type(sector),intent(in) :: sectorI
    real(8),intent(out)     :: sgn
    integer                 :: i_el,ii,iorb,iup,idw
    integer,dimension(2)    :: Indices
    integer,dimension(2)    :: Jndices
    integer,dimension(2*Ns) :: ib
    integer,dimension(eNs)  :: Nup,Ndw
    integer,dimension(iNs)  :: Npup,Npdw
    integer,dimension(Ns)   :: Nups,Ndws
    !
    sgn=0d0
    !
    if(KondoFlag)then
       i_el  = sectorI%H(1)%map(i)
       ib    = bdecomp(i_el,2*Ns)
       nup   = ib(1:eNs)
       ndw   = ib(eNs+1:2*eNs)
       npup  = ib(2*eNs+1:2*eNs+iNs)
       npdw  = ib(2*eNs+iNs+1:2*eNs+2*iNs)
       Nups  = [Nup,Npup]
       Ndws  = [Ndw,Npdw]
    else
       call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
       iup   = sectorI%H(1)%map(Indices(1))
       idw   = sectorI%H(2)%map(Indices(2))
       Nups  = Bdecomp(iup,Ns)
       Ndws  = Bdecomp(idw,Ns)       
    endif
    sgn = Nups(ipos)-Ndws(ipos)
    sgn = sgn/2d0
  end subroutine apply_op_Sz_single


  subroutine apply_op_Sz_range(i,sgn,ipos,sectorI) 
    integer, intent(in)     :: i,ipos(2)
    type(sector),intent(in) :: sectorI
    real(8),intent(out)     :: sgn
    integer                 :: i_el,ii,iorb,iup,idw
    integer,dimension(2)    :: Indices
    integer,dimension(2)    :: Jndices
    integer,dimension(2*Ns) :: ib
    integer,dimension(eNs)  :: Nup,Ndw
    integer,dimension(iNs)  :: Npup,Npdw
    integer,dimension(Ns)   :: Nups,Ndws
    !
    sgn=0d0
    !
    if(KondoFlag)then
       i_el  = sectorI%H(1)%map(i)
       ib    = bdecomp(i_el,2*Ns)
       nup   = ib(1:eNs)
       ndw   = ib(eNs+1:2*eNs)
       npup  = ib(2*eNs+1:2*eNs+iNs)
       npdw  = ib(2*eNs+iNs+1:2*eNs+2*iNs)
       Nups  = [Nup,Npup]
       Ndws  = [Ndw,Npdw]
    else
       call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
       iup   = sectorI%H(1)%map(Indices(1))
       idw   = sectorI%H(2)%map(Indices(2))
       Nups  = Bdecomp(iup,Ns)
       Ndws  = Bdecomp(idw,Ns)       
    endif
    !
    ! sgn = sum(Nups(ipos(1):ipos(2))-Ndws(ipos(1):ipos(2)))
    do ii=ipos(1),ipos(2)
       sgn = sgn + Nups(ii)-Ndws(ii)
    enddo
    sgn = sgn/2d0
  end subroutine apply_op_Sz_range



  subroutine apply_op_N(i,sgn,ipos,sectorI) 
    integer, intent(in)     :: i,ipos
    type(sector),intent(in) :: sectorI
    real(8),intent(out)     :: sgn
    integer                 :: i_el,ii,iorb,iup,idw
    integer,dimension(2)    :: Indices
    integer,dimension(2)    :: Jndices
    integer,dimension(2*Ns) :: ib
    integer,dimension(eNs)  :: Nup,Ndw
    integer,dimension(iNs)  :: Npup,Npdw
    integer,dimension(Ns)   :: Nups,Ndws
    !
    sgn=0d0
    !
    if(KondoFlag)then
       i_el = sectorI%H(1)%map(i)
       ib   = bdecomp(i_el,2*Ns)
       nup   = ib(1:eNs)
       ndw   = ib(eNs+1:2*eNs)
       npup  = ib(2*eNs+1:2*eNs+iNs)
       npdw  = ib(2*eNs+iNs+1:2*eNs+2*iNs)
       Nups  = [Nup,Npup]
       Ndws  = [Ndw,Npdw]
    else
       call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
       iup   = sectorI%H(1)%map(Indices(1))
       idw   = sectorI%H(2)%map(Indices(2))
       Nups  = Bdecomp(iup,Ns)
       Ndws  = Bdecomp(idw,Ns)       
    endif
    sgn = Nups(ipos)+Ndws(ipos)
  end subroutine apply_op_N



  subroutine build_op_Ns(i,Nups,Ndws,sectorI) 
    integer, intent(in)     :: i
    type(sector),intent(in) :: sectorI
    integer                 :: i_el,ii,iorb,iup,idw
    integer,dimension(2)    :: Indices
    integer,dimension(2*Ns) :: ib
    integer,dimension(eNs)  :: Nup,Ndw
    integer,dimension(iNs)  :: Npup,Npdw
    integer,dimension(Ns)   :: Nups,Ndws
    !
    if(KondoFlag)then
       i_el  = sectorI%H(1)%map(i)
       ib    = bdecomp(i_el,2*Ns)
       nup   = ib(1:eNs)
       ndw   = ib(eNs+1:2*eNs)
       npup  = ib(2*eNs+1:2*eNs+iNs)
       npdw  = ib(2*eNs+iNs+1:2*eNs+2*iNs)
       Nups  = [Nup,Npup]
       Ndws  = [Ndw,Npdw]
    else
       call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
       iup   = sectorI%H(1)%map(Indices(1))
       idw   = sectorI%H(2)%map(Indices(2))
       Nups  = Bdecomp(iup,Ns)
       Ndws  = Bdecomp(idw,Ns)
    endif
  end subroutine build_op_Ns










  subroutine get_Sector(QN,N,isector)
    integer,dimension(:) :: QN
    integer              :: N
    integer              :: isector
    integer              :: i,Nind,factor,Nup,Ndw
    if(KondoFlag)then
       Nup = QN(1)
       Ndw = QN(2)
       isector=getSector(Nup,Ndw)
       if(isector==0)then
          print*,Nup,Ndw
          stop "get_Sector ERROR: KondoFlag, looking for inexistent sector"
       endif
    else
       Nind = size(QN)
       Factor = N+1
       isector = 1
       do i=Nind,1,-1
          isector = isector + QN(i)*(Factor)**(Nind-i)
       enddo
    endif
  end subroutine get_Sector


  subroutine get_QuantumNumbers(isector,N,QN)
    integer                     :: isector,N
    integer,dimension(:)        :: QN
    integer                     :: i,count,Dim
    integer,dimension(size(QN)) :: QN_
    !
    if(KondoFlag)then
       QN(1)=getNup(isector)
       QN(2)=getNdw(isector)
    else
       Dim = size(QN)
       if(mod(Dim,2)/=0)stop "get_QuantumNumbers error: Dim%2 != 0"
       count=isector-1
       do i=1,Dim
          QN_(i) = mod(count,N+1)
          count      = count/(N+1)
       enddo
       QN = QN_(Dim:1:-1)
    endif
  end subroutine get_QuantumNumbers


  subroutine get_Nup(isector,Nup)
    integer              :: isector,Nup(1)
    integer              :: i,count
    integer,dimension(2) :: indices_
    if(KondoFlag)then
       Nup(1) = getNup(isector)
    else
       count=isector-1
       do i=1,2
          indices_(i) = mod(count,Ns+1)
          count      = count/(Ns+1)
       enddo
       Nup = indices_(2:2:-1)
    endif
  end subroutine get_Nup


  subroutine get_Ndw(isector,Ndw)
    integer              :: isector,Ndw(1)
    integer              :: i,count
    integer,dimension(2) :: indices_
    if(KondoFlag)then
       Ndw(1) = getNdw(isector)
    else
       count=isector-1
       do i=1,2
          indices_(i) = mod(count,Ns+1)
          count      = count/(Ns+1)
       enddo
       Ndw = indices_(1:1:-1)
    endif
  end subroutine get_Ndw


  subroutine  get_DimUp(isector,DimUps)
    integer                :: isector,DimUps(1)
    integer                :: Nups(1)
    call get_Nup(isector,Nups)
    DimUps(1) = binomial(Ns,Nups(1))
  end subroutine get_DimUp


  subroutine get_DimDw(isector,DimDws)
    integer                :: isector,DimDws(1)
    integer                :: Ndws(1)
    call get_Ndw(isector,Ndws)
    DimDws(1) = binomial(Ns,Ndws(1))
  end subroutine get_DimDw


  subroutine indices2state(ivec,Nvec,istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end subroutine indices2state

  subroutine state2indices(istate,Nvec,ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end subroutine state2indices


  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index


  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index












  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw.
  !
  !- build the map from the A-sector to \HHH
  !- get the list of states in \HHH corresponding to sector B twin of A
  !- return the ordering of B-states in \HHH with respect to those of A
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                       :: isector
    integer,dimension(:)          :: order
    type(sector)                  :: sectorH
    type(sector_map),dimension(2) :: H
    integer,dimension(2)          :: Indices,Istates
    integer,dimension(1)          :: DimUps,DimDws
    integer                       :: Dim,DimUp,DimDw
    integer                       :: i,iud
    !
    Dim = GetDim(isector)
    if(size(Order)/=Dim)stop "twin_sector_order error: wrong dimensions of *order* array"
    call get_DimUp(isector,DimUps)
    call get_DimDw(isector,DimDws)
    DimUp = product(DimUps)
    DimDw = product(DimDws)
    !
    call build_sector(isector,sectorH)
    do i=1,sectorH%Dim
       call state2indices(i,[sectorH%DimUps,sectorH%DimDws],Indices)
       forall(iud=1:2)Istates(iud) = sectorH%H(iud)%map(Indices(iud))
       Order(i) = flip_state( Istates )!flipped electronic state (GLOBAL state number {1:2^2Ns})
    enddo
    call delete_sector(sectorH)
    !
    call sort_array(Order) !sorted and changed the values from the global state numbers to the ones of the sector {1:DimUp*DimDw}
    !
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  !+------------------------------------------------------------------+
  function flip_state(istate) result(j)
    integer,dimension(2) :: istate
    integer              :: j
    integer,dimension(1) :: jups,jdws
    integer,dimension(2) :: dims
    !
    jups = istate(2:2)
    jdws = istate(1:1)
    dims = 2**Ns
    call indices2state([jups,jdws],Dims,j)
    !
  end function flip_state


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in)   :: isector
    integer              :: jsector
    integer,dimension(1) :: Iups,Idws
    call get_Nup(isector,iups)
    call get_Ndw(isector,idws)
    call get_Sector([idws,iups],Ns,jsector)
  end function get_twin_sector

















  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array





  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  subroutine print_conf(i,Ntot,advance)
    integer :: dim,i,j,Ntot
    logical :: advance
    integer :: ivec(Ntot)
    ivec = bdecomp(i,Ntot)
    write(LOGfile,"(A1)",advance="no")"|"
    write(LOGfile,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    if(advance)then
       write(LOGfile,"(A1)",advance="yes")">"
    else
       write(LOGfile,"(A1)",advance="no")">"
    endif
  end subroutine print_conf




end MODULE ED_SECTOR
