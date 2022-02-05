MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_TIMER,  only: start_timer,stop_timer
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
  public :: energy_lattice



  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:,:),allocatable :: sz2,n2
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
  real(8)                            :: gs_weight
  !
  real(8)                            :: peso
  real(8)                            :: norm
  !
  integer                            :: i,j,ii
  integer                            :: isector,jsector
  !
  complex(8),dimension(:),pointer       :: state_cvec
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



  subroutine energy_lattice()
    if(MPIMASTER)then
       write(LOGfile,"(A)")"Get Energy:"
       call start_timer()
    endif
    select case(ed_method)
    case default
       call lanc_energy()
    case ("lapack","full")
       call full_energy()
    end select
    if(MPIMASTER)call stop_timer(unit=LOGfile)
  end subroutine energy_lattice



  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine lanc_observables()
    integer                           :: iprob,istate,Nud(2,Ns),iud(2),jud(2),val
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)             :: IbUp,IbDw  ![Ns]
    real(8),dimension(Ns)             :: nup,ndw,Sz,nt
    !
    allocate(dens(Ns),dens_up(Ns),dens_dw(Ns))
    allocate(docc(Ns))
    allocate(magz(Ns),sz2(Ns,Ns),n2(Ns,Ns))
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
    !
    call es_trim_size(state_list,temp,cutoff)
    do istate=1,state_list%trimd_size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-(Ei-Egs)/temp)
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i = 1,sectorI%Dim
             gs_weight=peso*abs(state_cvec(i))**2
             call build_op_Ns(i,IbUp,Ibdw,sectorI)
             nup= dble(IbUp)
             ndw= dble(IbDw)
             sz = (nup-ndw)/2d0
             nt =  nup+ndw
             !
             !Evaluate averages of observables:
             do is=1,Ns
                dens(is)     = dens(is)      +  nt(is)*gs_weight
                dens_up(is)  = dens_up(is)   +  nup(is)*gs_weight
                dens_dw(is)  = dens_dw(is)   +  ndw(is)*gs_weight
                docc(is)     = docc(is)      +  nup(is)*ndw(is)*gs_weight
                magz(is)     = magz(is)      +  (nup(is)-ndw(is))*gs_weight
                sz2(is,is) = sz2(is,is)  +  (sz(is)*sz(is))*gs_weight
                n2(is,is)  = n2(is,is)   +  (nt(is)*nt(is))*gs_weight
                do js=is+1,Ns
                   sz2(is,js) = sz2(is,js)  +  (sz(is)*sz(js))*gs_weight
                   sz2(js,is) = sz2(js,is)  +  (sz(js)*sz(is))*gs_weight
                   n2(is,js)  = n2(is,js)   +  (nt(is)*nt(js))*gs_weight
                   n2(js,is)  = n2(js,is)   +  (nt(js)*nt(is))*gs_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*gs_weight
          enddo
          !
          call delete_sector(sectorI)
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    !
    !
    !
    !
    if(MPIMASTER)then
       call write_observables()
       write(LOGfile,"(A,100f18.12,f18.12)")&
            "dens=",(dens(is),is=1,Ns),sum(dens)
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
    endif
#endif
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2)
  end subroutine lanc_observables





  subroutine full_observables()
    integer                         :: iprob,i,j
    integer                         :: izero,istate
    integer                         :: isector,jsector
    integer                         :: idim,jdim
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite
    integer                         :: numstates
    integer                         :: r,m,k,val
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: boltzman_weight
    real(8)                         :: state_weight
    real(8)                         :: weight
    real(8)                         :: Ei
    real(8)                         :: norm,temp_
    integer                         :: Nud(2,Ns),iud(2),jud(2)
    integer,dimension(2*Ns_Ud)      :: Indices,Jndices
    integer,dimension(Ns_Ud)        :: Nups,Ndws
    integer,dimension(Ns)           :: IbUp,IbDw  ![Ns]
    real(8),dimension(Ns)           :: nup,ndw,Sz,nt
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
             weight = boltzman_weight*state_weight
             call build_op_Ns(i,IbUp,Ibdw,sectorI)
             nup= dble(IbUp)
             ndw= dble(IbDw)
             sz = (nup-ndw)/2d0
             nt =  nup+ndw
             !
             !Evaluate averages of observables:
             do io=1,Ns
                dens(io)     = dens(io)      +  nt(io)*weight
                dens_up(io)  = dens_up(io)   +  nup(io)*weight
                dens_dw(io)  = dens_dw(io)   +  ndw(io)*weight
                docc(io)     = docc(io)      +  nup(io)*ndw(io)*weight
                magz(io)     = magz(io)      +  (nup(io)-ndw(io))*weight
                sz2(io,io) = sz2(io,io)  +  (sz(io)*sz(io))*weight
                n2(io,io)  = n2(io,io)   +  (nt(io)*nt(io))*weight
                do jo=io+1,Ns
                   sz2(io,jo) = sz2(io,jo)  +  (sz(io)*sz(jo))*weight
                   sz2(jo,io) = sz2(jo,io)  +  (sz(jo)*sz(io))*weight
                   n2(io,jo)  = n2(io,jo)   +  (nt(io)*nt(jo))*weight
                   n2(jo,io)  = n2(jo,io)   +  (nt(jo)*nt(io))*weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*weight
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









  !+-------------------------------------------------------------------+
  !PURPOSE  : Get energy from the lattice problem.
  !+-------------------------------------------------------------------+
  subroutine lanc_energy()
    integer                           :: istate,iud(2),jud(2)
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)             :: Nup,Ndw,Sz
    complex(8),dimension(Nspin,Ns,Ns) :: Hij,Hloc
    complex(8),dimension(Nspin,Ns)    :: Hdiag
    !
    Egs     = state_list%emin
    ed_Ekin    = 0.d0
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    ed_Dk      = 0.d0
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = diagonal(Hloc(ispin,:,:))
    enddo
    !
    !
    call es_trim_size(state_list,temp,cutoff)
    do istate=1,state_list%trimd_size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-(Ei-Egs)/temp)
       peso = peso/zeta_function
       !
       !Master:
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i=1,sectorI%Dim
             !
             call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
             do ii=1,Ns_Ud
                mup = sectorI%H(ii)%map(Indices(ii))
                mdw = sectorI%H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb)
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nup =  Nups(1,:)
             Ndw =  Ndws(1,:)
             Sz  = 0.5d0*(Nup-Ndw)
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !> H_Imp: Diagonal Elements, i.e. local part
             do io=1,Ns
                ed_Eknot = ed_Eknot + Hdiag(1,io)*Nup(io)*gs_weight
                ed_Eknot = ed_Eknot + Hdiag(Nspin,io)*Ndw(io)*gs_weight
             enddo
             !
             !> H_imp: Off-diagonal elements, i.e. non-local part.
             !UP - hopping
             iup = Indices(1)
             idw = Indices(2)
             mup = sectorI%H(1)%map(iup)
             do io=1,Ns
                do jo=1,Ns
                   Jcondition = &
                        (Hij(1,io,jo)/=zero) .AND. &
                        (Nup(jo)==1) .AND. (Nup(io)==0)
                   if (Jcondition) then
                      call c(jo,mup,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jup = binary_search(sectorI%H(1)%map,k2)
                      j   = jup + (idw-1)*sectorI%DimUp
                      ed_Ekin = ed_Ekin + &
                           Hij(1,io,jo)*sg1*sg2*state_cvec(i)*conjg(state_cvec(j))*peso
                   endif
                enddo
             enddo
             !
             !DW - hopping
             iup = Indices(1)
             idw = Indices(2)
             mdw = sectorI%H(2)%map(idw)
             do io=1,Ns
                do jo=1,Ns
                   Jcondition = &
                        (Hij(Nspin,io,jo)/=zero) .AND. &
                        (Ndw(jo)==1) .AND. (Ndw(io)==0)
                   if (Jcondition) then
                      call c(jo,mdw,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jdw = binary_search(sectorI%H(2)%map,k2)
                      j   = iup + (jdw-1)*sectorI%DimUp
                      ed_Ekin = ed_Ekin + &
                           Hij(Nspin,io,jo)*sg1*sg2*state_cvec(i)*conjg(state_cvec(j))*peso
                   endif
                enddo
             enddo
             !
             !
             !SPIN-EXCHANGE Jx
             if(Jhflag.AND.Jx/=0d0)then
                iup = Indices(1)
                idw = Indices(2)
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)
                            jo = pack_indices(isite,jorb)
                            Jcondition=(&
                                 (io/=jo).AND.&
                                 (nup(jo)==1).AND.&
                                 (ndw(io)==1).AND.&
                                 (ndw(jo)==0).AND.&
                                 (nup(io)==0))
                            if(Jcondition)then
                               call c(io,mdw,k1,sg1)  !DW
                               call cdg(jo,k1,k2,sg2) !DW
                               jdw=binary_search(sectorI%H(2)%map,k2)
                               call c(jo,mup,k3,sg3)  !UP
                               call cdg(io,k3,k4,sg4) !UP
                               jup=binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                               ed_Dse = ed_Dse + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                               !
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             ! PAIR-HOPPING Jp
             if(Jhflag.AND.Jp/=0d0)then
                iup = Indices(1)
                idw = Indices(2)
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)
                            jo = pack_indices(isite,jorb)
                            Jcondition=(&
                                 (nup(jo)==1).AND.&
                                 (ndw(jo)==1).AND.&
                                 (ndw(io)==0).AND.&
                                 (nup(io)==0))
                            if(Jcondition)then
                               call c(jo,mdw,k1,sg1)       !c_jo_dw
                               call cdg(io,k1,k2,sg2)      !c^+_io_dw
                               jdw = binary_search(sectorI%H(2)%map,k2)
                               call c(jo,mup,k3,sg3)       !c_jo_up
                               call cdg(io,k3,k4,sg4)      !c^+_io_up
                               jup = binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                               ed_Dph = ed_Dph + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                               !
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             if(Jhflag.AND.Jk/=0d0)then
                iup = Indices(1)
                idw = Indices(2)
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)!a
                            jo = pack_indices(isite,jorb)!b
                            !
                            ![cdg_io c_jo]_up [cdg_jo c_io]_dw
                            Jcondition=(&
                                 (ndw(io)==1).AND.&
                                 (ndw(jo)==0).AND.&
                                 (nup(jo)==1).AND.&
                                 (nup(io)==0))
                            if(Jcondition)then
                               call c(io,mdw,k1,sg1)       !c_io.dw
                               call cdg(jo,k1,k2,sg2)      !c^+_jo.dw
                               jdw = binary_search(sectorI%H(2)%map,k2)
                               call c(jo,mup,k3,sg3)       !c_jo.up
                               call cdg(io,k3,k4,sg4)      !c^+_io.up
                               jup = binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jk/2d0*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                               ed_Dk = ed_Dk + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                            endif
                            !
                            ![cdg_jo c_io]_up [cdg_io c_jo]_dw
                            Jcondition=(&
                                 (ndw(jo)==1).AND.&
                                 (ndw(io)==0).AND.&
                                 (nup(io)==1).AND.&
                                 (nup(jo)==0))
                            if(Jcondition)then
                               call c(jo,mdw,k1,sg1)       !c_jo.dw
                               call cdg(io,k1,k2,sg2)      !c^+_io.dw
                               jdw = binary_search(sectorI%H(2)%map,k2)
                               call c(io,mup,k3,sg3)       !c_io.up
                               call cdg(jo,k3,k4,sg4)      !c^+_jo.up
                               jup = binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jk/2d0*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                               ed_Dk = ed_Dk + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                            endif
                            !
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             do iorb=1,Norb
                do isite=1,Nsites(iorb)          
                   io = pack_indices(isite,iorb)
                   ed_Epot = ed_Epot + Uloc(iorb)*nup(io)*ndw(io)*gs_weight
                enddo
             enddo
             !
             !Eust = \sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)
                            jo = pack_indices(isite,jorb)
                            ed_Epot = ed_Epot + Ust*(nup(io)*ndw(jo) + nup(jo)*ndw(io))*gs_weight
                            ed_Dust = ed_Dust + (nup(io)*ndw(jo) + nup(jo)*ndw(io))*gs_weight
                            !
                            ed_Epot = ed_Epot + (Ust-Jh)*(nup(io)*nup(jo) + ndw(io)*ndw(jo))*gs_weight
                            ed_Dund = ed_Dund + (nup(io)*nup(jo) + ndw(io)*ndw(jo))*gs_weight
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !
             if(Jhflag.AND.Jk/=0d0)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            ed_Epot = ed_Epot - 2*Jk*Sz(io)*Sz(jo)*gs_weight
                            ed_Dk = ed_Dk + Sz(io)*Sz(jo)*gs_weight
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                do iorb=1,Norb
                   do isite=1,Nsites(iorb)          
                      io = pack_indices(isite,iorb)
                      ed_Ehartree=ed_Ehartree - 0.5d0*Uloc(iorb)*(nup(io)+ndw(io))*gs_weight
                   enddo
                enddo
                !
                if(Norb>1)then
                   do iorb=1,Norb
                      do jorb=iorb+1,Norb
                         do isite=1,Nsites(iorb)
                            do jsite=1,Nsites(jorb)
                               if(isite/=jsite)cycle !local interaction only:
                               io = pack_indices(isite,iorb)
                               jo = pack_indices(isite,jorb)
                               ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(io)+ndw(io)+nup(jo)+ndw(jo))*gs_weight
                               ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(io)+ndw(io)+nup(jo)+ndw(jo))*gs_weight
                            enddo
                         enddo
                      enddo
                   enddo
                   !
                endif
             endif
             !
             !
          enddo
          call delete_sector(sectorI)         
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_Ekin)
       call Bcast_MPI(MpiComm,ed_Epot)
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
       call Bcast_MPI(MpiComm,ed_Dph)
       call Bcast_MPI(MpiComm,ed_Dse)
       call Bcast_MPI(MpiComm,ed_Dk)
    endif
#endif
    !
    !ed_Ekin = ed_Ekin/Ns        !Rescale to avoid linear increasing with Ns
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(MPIMASTER)then
       if(ed_verbose>2)then
          write(LOGfile,"(A,10f18.12)")"<K>     =",ed_Ekin
          write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
          write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
          write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
          write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
          write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
          write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
          write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
          write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
          write(LOGfile,"(A,10f18.12)")"Dk      =",ed_Dk
       endif
       !
       call write_energy()
    endif
    !
    !
  end subroutine lanc_energy



  subroutine full_energy()
    integer                           :: i,j
    integer                           :: izero,istate
    integer                           :: isector
    integer                           :: iorb,jorb,ispin
    integer                           :: numstates
    integer                           :: m,k1,k2,k3,k4
    real(8)                           :: sg1,sg2,sg3,sg4
    real(8)                           :: Ei
    real(8)                           :: boltzman_weight
    real(8)                           :: state_weight
    real(8)                           :: weight
    real(8)                           :: norm,temp_
    real(8),dimension(Nspin,Norb)     :: eloc
    complex(8),dimension(:),pointer      :: evec
    integer                           :: iud(2),jud(2)
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)             :: Nup,Ndw,Sz
    logical                           :: Jcondition
    complex(8),dimension(Nspin,Ns,Ns) :: Hij,Hloc
    complex(8),dimension(Nspin,Ns)    :: Hdiag
    integer            :: Iups(Ns_Ud)
    integer            :: Idws(Ns_Ud)
    !
    !
    ed_Ehartree= 0.d0
    ed_Ekin    = 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    ed_Dk      = 0.d0
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = diagonal(Hloc(ispin,:,:))
    enddo
    !
    temp_ = temp
    if(.not.finiteT)temp_=0.001d0
    !
    do isector=1,Nsectors
       call get_Nup(isector,Iups)
       call get_Ndw(isector,Idws)
       if(ed_filling/=0 .AND. (sum(Iups)+sum(Idws)/=ed_filling) )cycle
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
          do i=1,sectorI%Dim
             !
             call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
             do ii=1,Ns_Ud
                mup = sectorI%H(ii)%map(Indices(ii))
                mdw = sectorI%H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb)
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nup =  Nups(1,:)
             Ndw =  Ndws(1,:)
             Sz  = 0.5d0*(Nup-Ndw)
             !
             state_weight = conjg(evec(i))*evec(i)
             weight = boltzman_weight*state_weight
             !
             !
             !start evaluating the Tr(H_loc) to estimate potential energy
             !
             !LOCAL ENERGY
             !> H_Imp: Diagonal Elements, i.e. local part
             do io=1,Ns
                ed_Eknot = ed_Eknot + Hdiag(1,io)*Nup(io)*boltzman_weight*state_weight
                ed_Eknot = ed_Eknot + Hdiag(Nspin,io)*Ndw(io)*boltzman_weight*state_weight
             enddo
             !
             !> H_imp: Off-diagonal elements, i.e. non-local part.
             !UP - hopping
             iup = Indices(1)
             idw = Indices(2)
             mup = sectorI%H(1)%map(iup)
             do io=1,Ns
                do jo=1,Ns
                   Jcondition = &
                        (Hij(1,io,jo)/=zero) .AND. &
                        (Nup(jo)==1) .AND. (Nup(io)==0)
                   if (Jcondition) then
                      call c(jo,mup,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jup = binary_search(sectorI%H(1)%map,k2)
                      j   = jup + (idw-1)*sectorI%DimUp
                      ed_Ekin = ed_Ekin + &
                           Hij(1,io,jo)*sg1*sg2*evec(i)*conjg(evec(j))*boltzman_weight
                   endif
                enddo
             enddo
             !
             !DW - hopping
             iup = Indices(1)
             idw = Indices(2)
             mdw = sectorI%H(2)%map(idw)
             do io=1,Ns
                do jo=1,Ns
                   Jcondition = &
                        (Hij(Nspin,io,jo)/=zero) .AND. &
                        (Ndw(jo)==1) .AND. (Ndw(io)==0)
                   if (Jcondition) then
                      call c(jo,mdw,k1,sg1)
                      call cdg(io,k1,k2,sg2)
                      jdw = binary_search(sectorI%H(2)%map,k2)
                      j   = iup + (jdw-1)*sectorI%DimUp
                      ed_Ekin = ed_Ekin + &
                           Hij(Nspin,io,jo)*sg1*sg2*evec(i)*conjg(evec(j))*boltzman_weight
                   endif
                enddo
             enddo
             !
             !SPIN-EXCHANGE Jx
             if(Jhflag.AND.Jx/=0d0)then
                iup = Indices(1)
                idw = Indices(2)
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)
                            jo = pack_indices(isite,jorb)
                            Jcondition=(&
                                 (io/=jo).AND.&
                                 (nup(jo)==1).AND.&
                                 (ndw(io)==1).AND.&
                                 (ndw(jo)==0).AND.&
                                 (nup(io)==0))
                            if(Jcondition)then
                               call c(io,mdw,k1,sg1)  !DW
                               call cdg(jo,k1,k2,sg2) !DW
                               jdw=binary_search(sectorI%H(2)%map,k2)
                               call c(jo,mup,k3,sg3)  !UP
                               call cdg(io,k3,k4,sg4) !UP
                               jup=binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                               ed_Dse = ed_Dse + sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                               !
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             ! PAIR-HOPPING Jp
             if(Jhflag.AND.Jp/=0d0)then
                iup = Indices(1)
                idw = Indices(2)
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)
                            jo = pack_indices(isite,jorb)
                            Jcondition=(&
                                 (nup(jo)==1).AND.&
                                 (ndw(jo)==1).AND.&
                                 (ndw(io)==0).AND.&
                                 (nup(io)==0))
                            if(Jcondition)then
                               call c(jo,mdw,k1,sg1)       !c_jo_dw
                               call cdg(io,k1,k2,sg2)      !c^+_io_dw
                               jdw = binary_search(sectorI%H(2)%map,k2)
                               call c(jo,mup,k3,sg3)       !c_jo_up
                               call cdg(io,k3,k4,sg4)      !c^+_io_up
                               jup = binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                               ed_Dph = ed_Dph + sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                               !
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             if(Jhflag.AND.Jk/=0d0)then
                iup = Indices(1)
                idw = Indices(2)
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)!a
                            jo = pack_indices(isite,jorb)!b
                            !
                            ![cdg_io c_jo]_up [cdg_jo c_io]_dw
                            Jcondition=(&
                                 (ndw(io)==1).AND.&
                                 (ndw(jo)==0).AND.&
                                 (nup(jo)==1).AND.&
                                 (nup(io)==0))
                            if(Jcondition)then
                               call c(io,mdw,k1,sg1)       !c_io.dw
                               call cdg(jo,k1,k2,sg2)      !c^+_jo.dw
                               jdw = binary_search(sectorI%H(2)%map,k2)
                               call c(jo,mup,k3,sg3)       !c_jo.up
                               call cdg(io,k3,k4,sg4)      !c^+_io.up
                               jup = binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jk/2d0*sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                               ed_Dk = ed_Dk + sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                            endif
                            !
                            ![cdg_jo c_io]_up [cdg_io c_jo]_dw
                            Jcondition=(&
                                 (ndw(jo)==1).AND.&
                                 (ndw(io)==0).AND.&
                                 (nup(io)==1).AND.&
                                 (nup(jo)==0))
                            if(Jcondition)then
                               call c(jo,mdw,k1,sg1)       !c_jo.dw
                               call cdg(io,k1,k2,sg2)      !c^+_io.dw
                               jdw = binary_search(sectorI%H(2)%map,k2)
                               call c(io,mup,k3,sg3)       !c_io.up
                               call cdg(jo,k3,k4,sg4)      !c^+_jo.up
                               jup = binary_search(sectorI%H(1)%map,k4)
                               j = jup + (jdw-1)*sectorI%DimUp
                               !
                               ed_Epot = ed_Epot + Jk/2d0*sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                               ed_Dk = ed_Dk + sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                            endif
                            !
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             do iorb=1,Norb
                do isite=1,Nsites(iorb)          
                   io = pack_indices(isite,iorb)
                   ed_Epot = ed_Epot + Uloc(iorb)*nup(io)*ndw(io)*state_weight
                enddo
             enddo
             !
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)
                            jo = pack_indices(isite,jorb)
                            ed_Epot = ed_Epot + Ust*(nup(io)*ndw(jo) + nup(jo)*ndw(io))*state_weight*boltzman_weight
                            ed_Dust = ed_Dust + (nup(io)*ndw(jo) + nup(jo)*ndw(io))*state_weight*boltzman_weight
                            ed_Epot = ed_Epot + (Ust-Jh)*(nup(io)*nup(jo) + ndw(io)*ndw(jo))*state_weight*boltzman_weight
                            ed_Dund = ed_Dund + (nup(io)*nup(jo) + ndw(io)*ndw(jo))*state_weight*boltzman_weight
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             if(Jhflag.AND.Jk/=0d0)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            ed_Epot = ed_Epot - 2*Jk*Sz(io)*Sz(jo)*state_weight*boltzman_weight
                            ed_Dk = ed_Dk + Sz(io)*Sz(jo)*state_weight*boltzman_weight
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                do iorb=1,Norb
                   do isite=1,Nsites(iorb)          
                      io = pack_indices(isite,iorb)
                      ed_Ehartree=ed_Ehartree - 0.5d0*Uloc(iorb)*(nup(io)+ndw(io))*state_weight*boltzman_weight
                   enddo
                enddo
                !
                if(Norb>1)then
                   do iorb=1,Norb
                      do jorb=iorb+1,Norb
                         do isite=1,Nsites(iorb)
                            do jsite=1,Nsites(jorb)
                               if(isite/=jsite)cycle !local interaction only:
                               io = pack_indices(isite,iorb)
                               jo = pack_indices(isite,jorb)
                               ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(io)+ndw(io)+nup(jo)+ndw(jo))*state_weight*boltzman_weight
                               ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(io)+ndw(io)+nup(jo)+ndw(jo))*state_weight*boltzman_weight
                            enddo
                         enddo
                      enddo
                   enddo
                   !
                endif
                !
                !
             endif
          enddo
       enddo
       call delete_sector(sectorI)
       if(associated(evec))nullify(evec)
    enddo
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose>2)then
       write(LOGfile,"(A,10f18.12)")"<K>     =",ed_Ekin
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
       write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
       write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
       write(LOGfile,"(A,10f18.12)")"Dk      =",ed_Dk
    endif
    call write_energy()
    !
    !
  end subroutine full_energy





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
    write(unit,"(90F15.9)")xmu,temp,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp,Jk
    close(unit)

    unit = fopen("dens.ed",.true.)
    write(unit,"(90(F20.12,1X))")(dens(io),io=1,Ns)
    close(unit)         

    unit = fopen("dens_up.ed",.true.)
    write(unit,"(90(F20.12,1X))")(dens_up(io),io=1,Ns)
    close(unit)         

    unit = fopen("dens_dw.ed",.true.)
    write(unit,"(90(F20.12,1X))")(dens_dw(io),io=1,Ns)
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
       write(unit,"(90(F15.9,1X))")(sz2(io,jo),jo=1,Ns)
    enddo
    close(unit)         

    unit = fopen("N_corr.ed",.true.)
    do io=1,Ns
       write(unit,"(90(F15.9,1X))")(n2(io,jo),jo=1,Ns)
    enddo
    close(unit)         
    !
  end subroutine write_observables

  subroutine write_energy()
    integer :: unit

    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<K>",&
         reg(txtfy(2))//"<Hi>",&
         reg(txtfy(3))//"<V>=<Hi-Ehf>",&
         reg(txtfy(4))//"<Eloc>",&
         reg(txtfy(5))//"<Ehf>",&
         reg(txtfy(6))//"<Dst>",&
         reg(txtfy(7))//"<Dnd>",&
         reg(txtfy(8))//"<Dse>",&
         reg(txtfy(9))//"<Dph>",&
         reg(txtfy(10))//"<Dk>"
    close(unit)

    unit = fopen("energy_last.ed",.true.)
    write(unit,"(90F15.9)")&
         ed_Ekin,ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,&
         ed_Dust,ed_Dund,ed_Dse,ed_Dph,ed_Dk
    close(unit)
  end subroutine write_energy


end MODULE ED_OBSERVABLES
















