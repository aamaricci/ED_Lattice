MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SECTOR
  USE SF_TIMER
  USE SF_MISC, only: sort
  USE SF_IOTOOLS, only:free_unit,reg,file_length,save_array
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
    integer                          :: i,iud,iorb,jorb,ispin,jspin,unit
    integer,dimension(:),allocatable :: DimUps,DimDws
    integer                          :: Nup,Ndw
    integer                          :: Tstep,Dim
    integer,allocatable              :: Tord(:)
    logical                          :: Tbool
    !
    Jhflag=.FALSE.
    if(Norb>1.AND.any([Jx,Jp,Jk_z,Jk_xy]/=0d0))Jhflag=.TRUE.
    !
    KondoFlag=.FALSE.
    if(Nimp>0)KondoFlag=.TRUE.
    !
    call ed_checks_global
    !
    !>Setup Dimensions of the problem
    Ns = sum(Nsites(1:Norb))
    !
    Ns_Orb = Ns
    Ns_Ud  = 1
    Ns_Imp = Ns+Nimp
    !
    Nsectors = ((Ns_Orb+1)*(Ns_Orb+1))**Ns_Ud
    if(KondoFlag)Nsectors = (Ns_imp+1)*(Ns_Imp+1)-2
    !
    !
    if(MpiMaster)then
       write(LOGfile,"(A)")"Summary:"
       write(LOGfile,"(A)")"--------------------------------------------"
       write(LOGfile,"(A,I15)")'# of levels           = ',Ns_Imp
       write(LOGfile,"(A,I15)")'# of spins            = ',2*(Ns_Imp)
       write(LOGfile,"(A,I15)")'# of orbitals         = ',Norb
       write(LOGfile,"(A,I15)")'# of impurities       = ',Nimp
       write(LOGfile,"(A,5I4)")'# Nsites              = ',Nsites
       write(LOGfile,"(A,I15)")'# of  sectors         = ',Nsectors
    endif
    !

    Nup = (Ns_Imp)/2
    Ndw = (Ns_Imp)-Nup
    if(KondoFlag)then
       Dim = get_sector_dimension(Nup,Ndw)
       if(MpiMaster)then
          write(LOGfile,"(A,I20)")'Largest Sector(s)     = ',Dim
          write(LOGfile,"(A)")"--------------------------------------------"
       endif
    else
       allocate(DimUps(Ns_Ud))
       allocate(DimDws(Ns_Ud))
       do iud=1,Ns_Ud
          DimUps(iud) = get_sector_dimension(Nup)
          DimDws(iud) = get_sector_dimension(Ndw)
       enddo
       if(MpiMaster)then
          write(LOGfile,"(A,"//str(Ns_Ud)//"I8,2X,"//str(Ns_Ud)//"I8,I20)")&
               'Largest Sector(s)     = ',DimUps,DimDws,product(DimUps)*product(DimDws)
          write(LOGfile,"(A)")"--------------------------------------------"
       endif
    endif
    !
    !
    allocate(spH0ups(Ns_Ud))
    allocate(spH0dws(Ns_Ud))
    !
    !Allocate indexing arrays
    allocate(getCsector(Ns_Ud,2,Nsectors))  ;getCsector  =0
    allocate(getCDGsector(Ns_Ud,2,Nsectors));getCDGsector=0
    allocate(getDim(Nsectors));getDim=0
    allocate(getNup(Nsectors),getNdw(Nsectors));getNup=0;getNdw=0
    allocate(getSector(0:Ns+Nimp,0:Ns+Nimp));getSector=0
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
    Tstep = 1
    allocate(temperature_list(Tstep))
    temperature_list = temp
    if(finiteT)then
       inquire(file=trim(Tfile)//".restart",exist=Tbool)
       if(Tbool)then
          deallocate(temperature_list)
          write(LOGfile,"(A)")'Reading temperature list from file '//trim(Tfile)//".restart"
          Tstep = file_length(trim(Tfile)//".restart")
          open(free_unit(unit),file=trim(Tfile)//".restart")
          allocate(temperature_list(Tstep),Tord(Tstep))
          do i=1,Tstep
             read(unit,*)temperature_list(i)
          enddo
          close(unit)
          call sort(temperature_list,Tord)                !sort from smallest to largest
          temperature_list = temperature_list(Tstep:1:-1) !invert order
          call save_array(trim(Tfile)//".used",temperature_list)
          temp             = temperature_list(1)          !set actual Temp to largest
       endif
    endif
    !
    !allocate functions
    allocate(impGmats(Nspin,Ns_Imp,Ns_Imp,Lmats))
    allocate(impGreal(Nspin,Ns_Imp,Ns_Imp,Lreal))
    impGmats=zero
    impGreal=zero
    !
    allocate(impGMatrix(Nspin,Ns_Imp,Ns_Imp))
    allocate(SpinChiMatrix(Ns_Imp,Ns_Imp))
    allocate(OcMatrix(Ns_Imp))
    !
    global_chi_flag=.false.
    if(any([chispin_flag]))global_chi_flag=.true.
    global_gf_flag=.false.
    if(any([gf_flag]))global_gf_flag=.true.
    global_oc_flag=.false.
    if(any([oc_flag]))global_oc_flag=.true.
    !
    offdiag_gf_flag=offdiag_gf_flag.AND.Norb>1
    !
    allocate(spinChi_tau(Ns_Imp,Ns_Imp,0:Ltau))
    allocate(spinChi_w(Ns_Imp,Ns_Imp,Lreal))
    allocate(spinChi_iv(Ns_Imp,Ns_Imp,0:Lmats))
    !
    !allocate observables
    allocate(ed_dens(Ns_Imp))
    allocate(ed_docc(Ns_Imp))
    allocate(ed_dens_up(Ns_Imp))
    allocate(ed_dens_dw(Ns_Imp))
    allocate(ed_mag(Ns_Imp))
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
    if(KondoFlag.AND.ed_twin)then
       write(LOGfile,"(A)")"WARNING: can not yet use twin_sector with KondoFlag. Set to false."
       ed_twin=.false.
    endif
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
    integer                          :: i,iud,iorb,Nup,Ndw
    integer                          :: isector,jsector,gsector,ksector,lsector
    integer                          :: unit,status,istate,ishift,isign
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    type(sector) :: sectorI,sectorJ,sectorK,sectorG,sectorL
    !
    !Store full dimension of the sectors:
    if(KondoFlag)then
       isector=0
       do Nup=0,Ns_Imp
          do Ndw=0,Ns_Imp
             if(Nup==0.AND.Ndw==0)cycle
             if(Nup==Ns_Imp.AND.Ndw==Ns_Imp)cycle
             isector=isector+1
             getSector(Nup,Ndw)=isector
             getNup(isector)=Nup
             getNdw(isector)=Ndw
             getDim(isector)=get_sector_dimension(nup,ndw)
          enddo
       enddo
    else
       do isector=1,Nsectors
          call get_DimUp(isector,DimUps)
          call get_DimDw(isector,DimDws)
          DimUp = product(DimUps)
          DimDw = product(DimDws)  
          getDim(isector)  = DimUp*DimDw
       enddo
    endif
    !
    !
    do isector=1,Nsectors
       neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector)
    enddo
    !
    !
    twin_mask=.true.
    if(ed_twin)then             !not with KondoFlag
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
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jups(iud)=Jups(iud)-1;
          if(Jups(iud) < 0)cycle
          if(KondoFlag.AND.(Jups(iud) <= 0) .AND. (Jdws(iud)<=0) )cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,1,isector)=jsector
       enddo
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws
          Jups(iud)=Jups(iud)+1;
          if(Jups(iud) > Ns_Orb)cycle
          if( KondoFlag .AND. (Jups(iud) >=Ns_Orb) .AND. (Jdws(iud) >=Ns_Orb) )cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,1,isector)=jsector
       enddo
       !
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)-1
          if(Jdws(iud) < 0)cycle
          if( KondoFlag .AND. (Jups(iud) <= 0) .AND. (Jdws(iud)<=0) )cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,2,isector)=jsector
       enddo
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)+1;
          if(KondoFlag .AND. (Jups(iud) >=Ns_Orb) .AND. (Jdws(iud) >=Ns_Orb) )cycle
          if(Jdws(iud) > Ns_Orb)cycle
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
  elemental function get_sector_dimension(Nup,Ndw) result(dim)
    integer,intent(in)          :: Nup
    integer,optional,intent(in) :: Ndw
    integer                     :: i,ip,dim
    if(present(Ndw))then
       dim = 0
       do i=0,Nimp
          ip = Nimp-i
          dim = dim + binomial(Nimp,i)*binomial(Ns,Nup-i)*binomial(Ns,Ndw-ip)
       enddo
    else
       dim = binomial(Ns,Nup)
    endif
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
    if(n2>n1) then
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
