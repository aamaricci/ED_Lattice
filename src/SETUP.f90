MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SECTOR
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private



  public :: init_ed_structure
  public :: setup_global

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure()
    logical                          :: control
    integer                          :: i,iud,iorb,jorb,ispin,jspin
    integer,dimension(:),allocatable :: DimUps,DimDws
    !
    Jhflag=.FALSE.
    if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0.OR.Jk/=0d0))Jhflag=.TRUE.
    !
    call ed_checks_global
    !
    !>Setup Dimensions of the problem
    Ns = sum(Nsites(1:Norb))
    !
    Ns_Orb = Ns
    Ns_Ud  = 1
    Nsectors = ((Ns_Orb+1)*(Ns_Orb+1))**Ns_Ud
    !
    !
    allocate(DimUps(Ns_Ud))
    allocate(DimDws(Ns_Ud))
    do iud=1,Ns_Ud
       DimUps(iud) = get_sector_dimension(Ns_Orb,Ns_Orb/2)
       DimDws(iud) = get_sector_dimension(Ns_Orb,Ns_Orb-Ns_Orb/2)
    enddo
    if(MpiMaster)then
       write(LOGfile,"(A)")"Summary:"
       write(LOGfile,"(A)")"--------------------------------------------"
       write(LOGfile,"(A,I15)")'# of levels           = ',Ns
       write(LOGfile,"(A,I15)")'# of spins            = ',2*Ns
       write(LOGfile,"(A,I15)")'# of orbitals         = ',Norb
       write(LOGfile,"(A,5I4)")'# Nsites              = ',Nsites
       write(LOGfile,"(A,I15)")'# of  sectors         = ',Nsectors
       write(LOGfile,"(A,I15)")'Ns_Orb                = ',Ns_Orb
       write(LOGfile,"(A,I15)")'Ns_Ud                 = ',Ns_Ud
       write(LOGfile,"(A,"//str(Ns_Ud)//"I8,2X,"//str(Ns_Ud)//"I8,I20)")&
            'Largest Sector(s)     = ',DimUps,DimDws,product(DimUps)*product(DimDws)
       write(LOGfile,"(A)")"--------------------------------------------"
    endif
    !
    !
    allocate(spH0ups(Ns_Ud))
    allocate(spH0dws(Ns_Ud))
    !
    !Allocate indexing arrays
    allocate(getCsector(Ns_Ud,2,Nsectors))  ;getCsector  =0
    allocate(getCDGsector(Ns_Ud,2,Nsectors));getCDGsector=0
    !
    allocate(getDim(Nsectors));getDim=0
    !
    allocate(twin_mask(Nsectors))
    allocate(sectors_mask(Nsectors))
    allocate(neigen_sector(Nsectors))
    !
    !
    finiteT = ed_finite_temp
    !
    select case(ed_method)
    case default
       if(finiteT)then
          if(mod(lanc_nstates_sector,2)/=0)then
             lanc_nstates_sector=lanc_nstates_sector+1
             write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
          endif
          !
          lanc_nstates_total=lanc_nstates_sector*Nsectors+10
          if(mod(lanc_nstates_total,2)/=0)then
             lanc_nstates_total=lanc_nstates_total+1
             write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
          endif
          write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
          write(LOGfile,"(A,I6)")"Nstates   Total  = ", lanc_nstates_total
          write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
       else
          write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
       endif
    case('lapack','full')
       !
       if(finiteT)then
          write(LOGfile,"(A)")"Full ED finite T calculation"
       else
          ed_method='lanczos'
          lanc_nstates_total=1
          lanc_dim_threshold=product(DimUps)*product(DimDws)
          write(LOGfile,"(A)")"Full ED T=0 calculation. Set LANC_DIM_THRESHOLD to "//str(lanc_dim_threshold)
          if(lanc_dim_threshold>2**13)stop "Full ED T=0: LANC_DIM_THRESHOLD > 2**13=8192!"
       endif
    end select
    !
    !
    !allocate functions
    allocate(impSmats(Nspin,Ns,Ns,Lmats))
    allocate(impSreal(Nspin,Ns,Ns,Lreal))
    impSmats=zero
    impSreal=zero
    !
    allocate(impGmats(Nspin,Ns,Ns,Lmats))
    allocate(impGreal(Nspin,Ns,Ns,Lreal))
    impGmats=zero
    impGreal=zero
    !
    allocate(impG0mats(Nspin,Ns,Ns,Lmats))
    allocate(impG0real(Nspin,Ns,Ns,Lreal))
    impG0mats=zero
    impG0real=zero
    !
    chi_flag=.false.
    if(any([chispin_flag,chidens_flag,chipair_flag,chiexct_flag]))chi_flag=.true.
    !
    allocate(spinChi_tau(Ns,Ns,0:Ltau))
    allocate(spinChi_w(Ns,Ns,Lreal))
    allocate(spinChi_iv(Ns,Ns,0:Lmats))
    !
    ! allocate(densChi_tau(Ns,Ns,0:Ltau))
    ! allocate(densChi_w(Ns,Ns,Lreal))
    ! allocate(densChi_iv(Ns,Ns,0:Lmats))
    ! !
    ! allocate(pairChi_tau(Ns,Ns,0:Ltau))
    ! allocate(pairChi_w(Ns,Ns,Lreal))
    ! allocate(pairChi_iv(Ns,Ns,0:Lmats))
    ! !
    ! allocate(exctChi_tau(0:2,Ns,Ns,0:Ltau))
    ! allocate(exctChi_w(0:2,Ns,Ns,Lreal))
    ! allocate(exctChi_iv(0:2,Ns,Ns,0:Lmats))
    !
    !allocate observables
    allocate(ed_dens(Ns),ed_docc(Ns),ed_dens_up(Ns),ed_dens_dw(Ns),ed_mag(Ns))
    ed_dens=0d0
    ed_docc=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    ed_mag =0d0
    !
    allocate(Drude_weight(Norb))
    allocate(OptCond_w(Norb,Lreal))
    Drude_weight = 0d0
    OptCond_w    = zero
    !
  end subroutine init_ed_structure


  subroutine ed_checks_global
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if(Norb>5)stop "ED ERROR: Norb > 5 is currently not supported"
    !
    if(Nspin>1.AND.(ed_twin))then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
    end if
    !
    select case(ed_method)
    case default
       if(lanc_method=="lanczos")then
          if(ed_finite_temp)stop "ED ERROR: lanc_method==lanczos available only for T=0"
          if(lanc_nstates_sector>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_sector==1, T=0"
       endif
       !
       lanc_nstates_total=1
    case('lapack','full')
       if(MpiStatus.AND.mpiSIZE>1)stop "ED ERROR: ed_diag_type=FULL + MPIsize>1: not possible at the moment"
    end select
  end subroutine ed_checks_global




  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  subroutine setup_global
    integer                          :: DimUp,DimDw
    integer                          :: DimUps(Ns_Ud),DimDws(Ns_Ud)
    integer                          :: Indices(2*Ns_Ud),Jndices(2*Ns_Ud)
    integer                          :: Nups(Ns_ud),Ndws(Ns_ud)
    integer                          :: Jups(Ns_ud),Jdws(Ns_ud)
    integer                          :: i,iud,iorb
    integer                          :: isector,jsector,gsector,ksector,lsector
    integer                          :: unit,status,istate,ishift,isign
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    type(sector) :: sectorI,sectorJ,sectorK,sectorG,sectorL
    !
    !Store full dimension of the sectors:
    do isector=1,Nsectors
       call get_DimUp(isector,DimUps)
       call get_DimDw(isector,DimDws)
       DimUp = product(DimUps)
       DimDw = product(DimDws)  
       getDim(isector)  = DimUp*DimDw
    enddo
    !
    !
    do isector=1,Nsectors
       neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector) !init every sector to required eigenstates
    enddo
    !
    !
    twin_mask=.true.
    if(ed_twin)then
       do isector=1,Nsectors
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          if(any(Nups < Ndws))twin_mask(isector)=.false.
       enddo
       write(LOGfile,"(A,I6,A,I9)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    !
    !
    getCsector  = 0
    getCDGsector= 0
    do isector=1,Nsectors
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       !
       !UPs:
       !DEL:
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jups(iud)=Jups(iud)-1; if(Jups(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,1,isector)=jsector
       enddo
       !ADD
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws
          Jups(iud)=Jups(iud)+1; if(Jups(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,1,isector)=jsector
       enddo
       !
       !DWs:
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)-1; if(Jdws(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,2,isector)=jsector
       enddo
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)+1; if(Jdws(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,2,isector)=jsector
       enddo
    enddo
    return
  end subroutine setup_global



  !##################################################################
  !##################################################################
  !SECTOR PROCEDURES - Sectors,Nup,Ndw,DimUp,DimDw,...
  !##################################################################
  !##################################################################
  elemental function get_sector_dimension(n,np) result(dim)
    integer,intent(in) :: n,np
    integer            :: dim
    dim = binomial(n,np)
  end function get_sector_dimension






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



end MODULE ED_SETUP
