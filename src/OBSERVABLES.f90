MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
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
  public :: observables_lattice


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
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_lattice()
    if(MPIMASTER)then
       write(LOGfile,"(A)")"Get Observables:"
       call start_timer()
    endif
    select case(ed_method)
    case default
       call lanc_observables()
    case ("lapack","full")
       call full_observables()
    end select
    if(MPIMASTER)call stop_timer(unit=LOGfile)
  end subroutine observables_lattice






  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine lanc_observables()
    integer                             :: iprob,istate,Nud(2,Ns),iud(2),jud(2),val,iimp
    integer,dimension(2)                :: Indices,Jndices
    integer,dimension(1,Ns)             :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: IbUp,IbDw  ![Ns]
    real(8),dimension(Ns)               :: nup,ndw,Sz,nt
    complex(8),dimension(:),allocatable :: state_cvec
    real(8)                             :: boltzman_weight
    real(8)                             :: state_weight
    !
    allocate(dens(Ns),dens_up(Ns),dens_dw(Ns))
    allocate(docc(Ns))
    allocate(magz(Ns),sz2(Ns,Ns),n2(Ns,Ns))
    allocate(dens_ImpUp(2,eNs))
    allocate(dens_ImpDw(2,eNs))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    dens_impup = 0.d0
    dens_impdw = 0.d0
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
             state_weight=abs(state_cvec(i))**2
             call build_op_Ns(i,IbUp,Ibdw,sectorI)
             nup= dble(IbUp)
             ndw= dble(IbDw)
             sz = (nup-ndw)/2d0
             nt =  nup+ndw
             !
             !Evaluate averages of observables:
             do is=1,Ns
                dens(is)     = dens(is)      +  nt(is)*state_weight*boltzman_weight
                dens_up(is)  = dens_up(is)   +  nup(is)*state_weight*boltzman_weight
                dens_dw(is)  = dens_dw(is)   +  ndw(is)*state_weight*boltzman_weight
                docc(is)     = docc(is)      +  nup(is)*ndw(is)*state_weight*boltzman_weight
                magz(is)     = magz(is)      +  (nup(is)-ndw(is))*state_weight*boltzman_weight
                sz2(is,is) = sz2(is,is)  +  (sz(is)*sz(is))*state_weight*boltzman_weight
                n2(is,is)  = n2(is,is)   +  (nt(is)*nt(is))*state_weight*boltzman_weight
                do js=is+1,Ns
                   sz2(is,js) = sz2(is,js)  +  (sz(is)*sz(js))*state_weight*boltzman_weight
                   sz2(js,is) = sz2(js,is)  +  (sz(js)*sz(is))*state_weight*boltzman_weight
                   n2(is,js)  = n2(is,js)   +  (nt(is)*nt(js))*state_weight*boltzman_weight
                   n2(js,is)  = n2(js,is)   +  (nt(js)*nt(is))*state_weight*boltzman_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*state_weight*boltzman_weight

             if(any([Jk_z,Jk_xy]/=0d0))then
                do iimp=1,iNs
                   if(nup(eNs+iimp)==1)then
                      do is=1,eNs
                         dens_ImpUp(1,is)  = dens_ImpUp(1,is)   +  nup(is)*state_weight*boltzman_weight
                         dens_ImpUp(2,is)  = dens_ImpUp(2,is)   +  ndw(is)*state_weight*boltzman_weight
                      enddo
                   endif
                   if(ndw(eNs+iimp)==1)then
                      do is=1,eNs
                         dens_ImpDw(1,is)  = dens_ImpDw(1,is)   +  nup(is)*state_weight*boltzman_weight
                         dens_ImpDw(2,is)  = dens_ImpDw(2,is)   +  ndw(is)*state_weight*boltzman_weight
                      enddo
                   endif
                enddo
             endif
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
    !
    !
    if(MPIMASTER)then
       call write_observables()
       if(KondoFlag)then
          if(iNs==eNs)then
             write(LOGfile,"(A,"//str(iNs)//"f15.8)")&
                  "dens=",(dens(eNs+is),is=1,iNs)
          else
             write(LOGfile,"(A,"//str(max(1,Jkindx(1)-1))//"A15,"//str(iNs)//"f15.8)")&
                  "dens=",("",is=1,max(1,Jkindx(1)-1)),(dens(eNs+is),is=1,iNs)
          endif
          write(LOGfile,"(A,100f15.8,f15.8)")&
               "dens=",(dens(is),is=1,eNs),sum(dens)
       else
          write(LOGfile,"(A,100f15.8,f15.8)")&
               "dens=",(dens(is),is=1,Ns),sum(dens)
       endif
       if(Nspin==2)then
          write(LOGfile,"(A,10f18.12,A)")&
               "magZ=",(magz(is),is=1,Ns)
       endif
       !
       ed_dens_up = dens_up
       ed_dens_dw = dens_dw
       ed_dens    = dens
       ed_docc    = docc
       ed_mag     = dens_up-dens_dw
    endif
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
       call Bcast_MPI(MpiComm,ed_mag)
       call Bcast_MPI(MpiComm,dens_ImpUp)
       call Bcast_MPI(MpiComm,dens_Impdw)
    endif
#endif
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2,dens_impup,dens_impdw)
  end subroutine lanc_observables





  subroutine full_observables()
    integer                         :: iprob,i,j
    integer                         :: izero,istate
    integer                         :: isector,jsector
    integer                         :: idim,jdim
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite
    integer                         :: r,m,k,val
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: boltzman_weight
    real(8)                         :: state_weight
    real(8)                         :: weight
    real(8)                         :: Ei
    real(8)                         :: temp_
    integer                         :: Nud(2,Ns),iud(2),jud(2)
    integer,dimension(2)      :: Indices,Jndices
    integer,dimension(1)        :: Nups,Ndws
    integer,dimension(Ns)       :: IbUp,IbDw  ![Ns]
    real(8),dimension(Ns)       :: nup,ndw,Sz,nt
    complex(8),dimension(:),pointer :: evec
    !
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Ns),dens_up(Ns),dens_dw(Ns))
    allocate(docc(Ns))
    allocate(magz(Ns),sz2(Ns,Ns),n2(Ns,Ns))

    !
    egs     = gs_energy
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    temp_ = temp
    if(.not.finiteT)temp_=0.001d0
    !
    do isector=1,Nsectors
       call get_Nup(isector,nups)
       call get_Ndw(isector,ndws)
       if(ed_filling/=0 .AND. (sum(Nups)+sum(Ndws)/=ed_filling) )cycle

       call build_sector(isector,sectorI)
       !
       do istate=1,sectorI%Dim
          Ei=espace(isector)%e(istate)
          boltzman_weight=exp(-Ei/temp_)/zeta_function
          if(boltzman_weight < cutoff)cycle
          !
          evec => espace(isector)%M(:,istate)
          !
          do i=1,sectorI%Dim
             state_weight=conjg(evec(i))*evec(i)
             call build_op_Ns(i,IbUp,Ibdw,sectorI)
             nup= dble(IbUp)
             ndw= dble(IbDw)
             sz = (nup-ndw)/2d0
             nt =  nup+ndw
             !
             !Evaluate averages of observables:
             do io=1,Ns
                dens(io)     = dens(io)      +  nt(io)*boltzman_weight*state_weight
                dens_up(io)  = dens_up(io)   +  nup(io)*boltzman_weight*state_weight
                dens_dw(io)  = dens_dw(io)   +  ndw(io)*boltzman_weight*state_weight
                docc(io)     = docc(io)      +  nup(io)*ndw(io)*boltzman_weight*state_weight
                magz(io)     = magz(io)      +  (nup(io)-ndw(io))*boltzman_weight*state_weight
                sz2(io,io) = sz2(io,io)  +  (sz(io)*sz(io))*boltzman_weight*state_weight
                n2(io,io)  = n2(io,io)   +  (nt(io)*nt(io))*boltzman_weight*state_weight
                do jo=io+1,Ns
                   sz2(io,jo) = sz2(io,jo)  +  (sz(io)*sz(jo))*boltzman_weight*state_weight
                   sz2(jo,io) = sz2(jo,io)  +  (sz(jo)*sz(io))*boltzman_weight*state_weight
                   n2(io,jo)  = n2(io,jo)   +  (nt(io)*nt(jo))*boltzman_weight*state_weight
                   n2(jo,io)  = n2(jo,io)   +  (nt(jo)*nt(io))*boltzman_weight*state_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*boltzman_weight*state_weight
          enddo
          !
       enddo
       call delete_sector(sectorI)
       if(associated(evec))nullify(evec)
    enddo
    !
    call write_observables()
    write(LOGfile,"(A,100f18.12,f18.12)")&
         "dens=",(dens(is),is=1,Ns),sum(dens)
    if(Nspin==2)then
       write(LOGfile,"(A,10f18.12,A)")&
            "magZ=",(magz(is),is=1,Ns)
    endif
    !
    !
    ed_dens_up = dens_up
    ed_dens_dw = dens_dw
    ed_dens    = dens
    ed_docc    = docc
    ed_mag     = dens_up-dens_dw
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2)
  end subroutine full_observables










  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: io,jo,iorb
    unit = fopen("parameters_last.ed",.true.)
    write(unit,"(90F15.9)")xmu,temp,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp,Jk_z,Jk_xy
    close(unit)

    unit = fopen("dens.ed",.true.)
    write(unit,"(90(F20.12,1X))")(dens(io),io=1,Ns)
    close(unit)

    unit = fopen("docc.ed",.true.)
    write(unit,"(90(F20.12,1X))")(docc(io),io=1,Ns)
    close(unit)

    unit = fopen("magz.ed",.true.)
    write(unit,"(90(F20.12,1X))")(magz(io),io=1,Ns)
    close(unit)         

    unit = fopen("egs.ed",.true.)
    write(unit,*)egs
    close(unit)

    unit = fopen("Sz_corr.ed",.true.)
    do io=1,Ns
       write(unit,"(90(F15.9,1X))")(sz2(io,jo),jo=1,io)
    enddo
    close(unit)         

    unit = fopen("N_corr.ed",.true.)
    do io=1,Ns
       write(unit,"(90(F15.9,1X))")(n2(io,jo),jo=1,io)
    enddo
    close(unit)         


    if(any([Jk_z,Jk_xy]/=0d0))then
       unit = fopen("dens_up_w_impup.ed",.true.)
       write(unit,"(90(F20.12,1X))")(dens_ImpUp(1,io),io=1,eNs)
       close(unit)
       unit = fopen("dens_dw_w_impup.ed",.true.)
       write(unit,"(90(F20.12,1X))")(dens_ImpUp(2,io),io=1,eNs)
       close(unit)         

       unit = fopen("dens_up_w_impdw.ed",.true.)
       write(unit,"(90(F20.12,1X))")(dens_ImpDw(1,io),io=1,eNs)
       close(unit)
       unit = fopen("dens_dw_w_impdw.ed",.true.)
       write(unit,"(90(F20.12,1X))")(dens_ImpDw(2,io),io=1,eNs)
       close(unit)         
    endif
  end subroutine write_observables


end MODULE ED_OBSERVABLES
















