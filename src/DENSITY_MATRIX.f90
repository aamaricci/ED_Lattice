MODULE ED_DENSITY_MATRIX
   USE SF_CONSTANTS, only:zero,pi,xi
   USE SF_IOTOOLS, only:free_unit,reg,str
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
   public :: ed_get_density_matrix

   interface ed_get_density_matrix
      !! Generic interface to retrieve all density matrices: if you omit the
      !! number of sites you want, you'll get the dm for the full lattice.
      !! Otherwise you'll get a reduced density matrix for Nred sites. There
      !! is currently no way to specify which exact sites you want to keep.
      module procedure :: get_lattice_density_matrix !(DM[,print_flag])
      module procedure :: get_reduced_density_matrix !(DM,Nred[,print])
   end interface ed_get_density_matrix

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
      complex(8),allocatable :: dm(:,:)
      integer                :: Nr
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
      if(MPIMASTER)then
         call ed_get_density_matrix(dm,doprint=ed_print_dm)
         call stop_timer(unit=LOGfile)
      endif
      !
      if(MPIMASTER.and.rdm_flag)then
         write(LOGfile,"(A)")"Get reduced density matrices:"
         call start_timer()
         do Nr=Nsites(1)-1,1,-1
            call ed_get_density_matrix(dm,Nr,doprint=ed_print_dm)
         enddo
         call stop_timer(unit=LOGfile)
      endif
   end subroutine eval_dm_lattice






   !+-------------------------------------------------------------------+
   !PURPOSE  : Lanc method
   !+-------------------------------------------------------------------+
   subroutine lanc_density_matrix()
      integer                             :: istate,idm,jdm
      integer,dimension(2)                :: Indices,Jndices
      integer,dimension(1,Ns)             :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
      integer,dimension(Ns)               :: IbUp,IbDw,JbUp,JbDw  ![Ns]
      real(8),dimension(Ns)               :: nup,ndw,Sz,nt
      complex(8),dimension(:),allocatable :: state_cvec
      real(8)                             :: boltzman_weight
      real(8)                             :: state_weight
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
               idm = bjoin([IbUp,IbDw],2*Ns) + 1
               do j = 1,sectorI%Dim
                  call build_op_Ns(j,JbUp,JbDw,sectorI)
                  jdm = bjoin([JbUp,JbDw],2*Ns) + 1
                  ed_dm_lattice(idm,jdm) = ed_dm_lattice(idm,jdm) + &
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




   !+-------------------------------------------------------------------+
   !PURPOSE  : Full diagonalization (lapack)
   !+-------------------------------------------------------------------+
   subroutine full_density_matrix()
      integer                         :: istate,idm,jdm
      real(8)                         :: boltzman_weight
      real(8)                         :: state_weight
      real(8)                         :: temp_
      integer,dimension(2)            :: Indices,Jndices
      integer,dimension(1)            :: Nups,Ndws
      integer,dimension(Ns)           :: IbUp,IbDw,JbUp,JbDw ![Ns]
      complex(8),dimension(:),pointer :: evec
      !
      !
      temp_ = temp
      if(.not.finiteT)temp_=0.001d0
      !
      do isector=1,Nsectors
         !
         call get_Nup(isector,nups)
         call get_Ndw(isector,ndws)
         if(ed_filling/=0 .AND. (sum(Nups)+sum(Ndws)/=ed_filling) )cycle
         !
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
               idm = bjoin([IbUp,IbDw],2*Ns) + 1
               do j = 1,sectorI%Dim
                  call build_op_Ns(j,JbUp,JbDw,sectorI)
                  jdm = bjoin([JbUp,JbDw],2*Ns) + 1
                  ed_dm_lattice(idm,jdm) = ed_dm_lattice(idm,jdm) + &
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





   !+---------------------------------------------------------------------+
   !PURPOSE : retrieve the lattice density matrix from the internal scope
   !+---------------------------------------------------------------------+
   subroutine get_lattice_density_matrix(dm,doprint)
      complex(8),allocatable,intent(out)           :: dm(:,:)
      logical               ,intent(in) ,optional  :: doprint
      logical                                      :: doprint_
      !
      doprint_=.false.; if(present(doprint)) doprint_=doprint
      !
      !Check if density matrix is present in the global scope
      if(.not.allocated(ed_dm_lattice))then
         stop "ERROR: lattice density matrix is not allocated"
      endif
      !
      dm = ed_dm_lattice
      !
      !Print to file (if requested)
      if(doprint_)then
         call print_density_matrix(dm,4**Ns)
      endif
      !
   end subroutine get_lattice_density_matrix



   !+---------------------------------------------------------------------+
   !PURPOSE : reduce the density matrix by tracing out (Nsites-Nred) sites
   !+---------------------------------------------------------------------+
   ! NB: the spin-factorization of the loops is paramount,
   !     in order to trace along the proper matrix elements
   subroutine get_reduced_density_matrix(rdm,Nred,doprint)
      complex(8),dimension(:,:),allocatable,intent(out)          :: rdm
      integer                              ,intent(in)           :: Nred
      logical                              ,intent(in),optional  :: doprint
      logical                                                    :: doprint_
      integer    :: i,j,io,jo,iUP,iDW,jUP,jDW
      integer    :: iBITup,iBITdw,jBITup,jBITdw
      integer    :: iREDup,iREDdw,jREDup,jREDdw
      integer    :: iTrUP,iTrDW,jTrUP,jTrDW,Ntrace,Nrdm
      !
      if(Nred>Nsites(1))then
         stop "ERROR: cannot request a density matrix reduced to more sites then available in the lattice"
      endif
      !
      doprint_=.false.; if(present(doprint)) doprint_=doprint
      !
      !Check if density matrix is present in the global scope
      if(.not.allocated(ed_dm_lattice))then
         stop "ERROR: lattice density matrix is not allocated"
      endif
      !
      associate (dm=>ed_dm_lattice) ! alias for readability
         !
         Nrdm=Norb*Nred !Number of bits associated to the reduced system
         Ntrace=Ns-Nrdm !Number of bits to trace away from the state
         allocate(rdm(4**Nrdm,4**Nrdm))
         rdm=zero
         !
         !Trace the dm to the requested number of sites
         do iUP = 1,2**Ns
            do iDW = 1,2**Ns
               i = iUP + (iDW-1)*2**Ns
               iBITup = iup-1
               iBITdw = idw-1
               iREDup = Ibits(iBITup,0,Nrdm)
               iREDdw = Ibits(iBITdw,0,Nrdm)
               iTrUP  = Ibits(iBITup,Nrdm,Ntrace)
               iTrDW  = Ibits(iBITdw,Nrdm,Ntrace)
               do jUP = 1,2**Ns
                  do jDW = 1,2**Ns
                     j = jUP + (jDW-1)*2**Ns
                     jBITup = jup-1
                     jBITdw = jdw-1
                     jREDup = Ibits(jBITup,0,Nrdm)
                     jREDdw = Ibits(jBITdw,0,Nrdm)
                     jTrUP  = Ibits(jBITup,Nrdm,Ntrace)
                     jTrDW  = Ibits(jBITdw,Nrdm,Ntrace)
                     if(jTrUP/=iTrUP.or.jTrDW/=iTrDW)cycle
                     io = (iREDup+1) + iREDdw*2**Nrdm
                     jo = (jREDup+1) + jREDdw*2**Nrdm
                     rdm(io,jo) = rdm(io,jo) + dm(i,j)
                  enddo
               enddo
            enddo
         enddo
         !
      end associate
      !
      !Print to file (if requested)
      if(doprint_)then
         call print_density_matrix(rdm,4**Nrdm)
      endif
      !
   end subroutine get_reduced_density_matrix






   !####################################################################
   !                           I/O ROUTINES
   !####################################################################

   !+-------------------------------------------------------------------+
   !PURPOSE  : write density matrices to file
   !+-------------------------------------------------------------------+
   subroutine print_density_matrix(dm,N)
      integer                  ,intent(in)            :: N
      complex(8),dimension(:,:),intent(in)            :: dm
      integer                                         :: unit
      character(len=64)                               :: suffix
      integer                                         :: io,jo
      !
      if(size(dm,1)/=N.OR.size(dm,2)/=N)then
         stop "ERROR: actual dm argument has incogruent size wrt explicitly passed N"
      endif
      !
      if(N<4**Ns)then
         suffix = "reduced_density_matrix"//"_rank"//reg(str(N))//".dat"
         ! TODO: better to write the number of sites, i.e. nint(log4(N))
         !       probably something like str(nint(log4(N)))//"sites"
         !       For now let's stay this way to be consistent with CDMFT
         !       code current version, then we'll change both.
      else
         suffix = "lattice_density_matrix.dat"
      endif
      !
      unit = free_unit()
      open(unit,file=suffix,action="write",position="rewind",status='unknown')
      !
      do io=1,N
         write(unit,"(90(F15.9,1X))") (dreal(dm(io,jo)),jo=1,N)
      enddo
      write(unit,*)
      !
      if(any(dimag(dm)/=0d0))then
         do io=1,N
            write(unit,"(90(F15.9,1X))") (dimag(dm(io,jo)),jo=1,N)
         enddo
         write(unit,*)
      endif
      !
      close(unit)
      !
   end subroutine print_density_matrix


end MODULE ED_DENSITY_MATRIX
















