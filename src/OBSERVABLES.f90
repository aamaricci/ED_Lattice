MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
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
  public :: observables_impurity
  public :: local_energy_impurity



  logical,save                       :: iolegend=.true.
  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:,:),allocatable :: sz2,n2
  real(8),dimension(:,:),allocatable :: zimp,simp
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
  real(8),dimension(:),allocatable   :: vvinit
  real(8),dimension(:),pointer       :: state_cvec
  logical                            :: Jcondition
  !
  type(sector)                       :: sectorI,sectorJ


contains 


  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    integer                           :: iprob,istate,Nud(2,Ns),iud(2),jud(2),val
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)             :: IbUp,IbDw  ![Ns]
    real(8),dimension(Ns)             :: nup,ndw,Sz,nt
    !
    allocate(dens(Ns),dens_up(Ns),dens_dw(Ns))
    allocate(docc(Ns))
    allocate(magz(Ns),sz2(Ns,Ns),n2(Ns,Ns))
    allocate(simp(Nspin,Ns),zimp(Nspin,Ns))
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
    do istate=1,state_list%size
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
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
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
             ! !Configuration probability
             ! iprob=1
             ! do iorb=1,Norb
             !    iprob=iprob+nint(nt(iorb))*3**(iorb-1)
             ! end do
             ! Prob(iprob) = Prob(iprob) + gs_weight
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
    !IMPURITY DENSITY MATRIX
    if(allocated(imp_density_matrix)) deallocate(imp_density_matrix)
    allocate(imp_density_matrix(Nspin,Ns,Ns));imp_density_matrix=zero
    do istate=1,state_list%size
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
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i=1,sectorI%Dim
             call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
             !
             call build_op_Ns(i,Nud(1,:),Nud(2,:),sectorI)
             !
             !Diagonal densities
             do ispin=1,Nspin
                do is=1,Ns
                   imp_density_matrix(ispin,is,is) = &
                        imp_density_matrix(ispin,is,is) + &
                        peso*Nud(ispin,is)*(state_cvec(i))*state_cvec(i)
                enddo
             enddo
             !
             !Off-diagonal
             do ispin=1,Nspin
                do is=1,Ns
                   do js=1,Ns
                      if((Nud(ispin,js)==1).and.(Nud(ispin,is)==0))then
                         iud(1) = sectorI%H(1)%map(Indices(1))
                         iud(2) = sectorI%H(2)%map(Indices(2))
                         call c(js,iud(ispin),r,sgn1)
                         call cdg(is,r,k,sgn2)
                         Jndices = Indices
                         Jndices(1+(ispin-1)*Ns_Ud) = &
                              binary_search(sectorI%H(1+(ispin-1)*Ns_Ud)%map,k)
                         call indices2state(Jndices,[sectorI%DimUps,sectorI%DimDws],j)
                         !
                         imp_density_matrix(ispin,is,js) = imp_density_matrix(ispin,is,js) + &
                              peso*sgn1*state_cvec(i)*sgn2*(state_cvec(j))
                      endif
                   enddo
                enddo
             enddo
          enddo
          call delete_sector(sectorI)         
       endif
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
       call get_szr
       if(iolegend)call write_legend
       call write_observables()
       write(LOGfile,"(A,10f18.12,f18.12)")&
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
    deallocate(simp,zimp)
  end subroutine observables_impurity





  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_impurity()
    integer                           :: istate,iud(2),jud(2)
    integer,dimension(2*Ns_Ud)        :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb)   :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)             :: Nup,Ndw
    complex(8),dimension(Nspin,Ns,Ns) :: Hij,Hloc
    complex(8),dimension(Nspin,Ns)    :: Hdiag
    !
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = diagonal(Hloc(ispin,:,:))
    enddo
    !
    !
    do istate=1,state_list%size
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
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
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
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nup =  Nups(1,:)!Breorder(Nups)
             Ndw =  Ndws(1,:)!Breorder(Ndws)
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !> H_Imp: Diagonal Elements, i.e. local part
             do io=1,Ns
                ed_Eknot = ed_Eknot + Hdiag(1,io)*Nup(io)*gs_weight
                ed_Eknot = ed_Eknot + Hdiag(Nspin,io)*Ndw(io)*gs_weight
             enddo
             ! !> H_imp: Off-diagonal elements, i.e. non-local part.
             !
             !UP - hopping
             iup = Indices(1)
             mup = sectorI%H(1)%map(iup)
             do iorb=1,Norb
                do jorb=1,Norb
                   do isite=1,Nsites(iorb)
                      do jsite=1,Nsites(jorb)
                         if(isite/=jsite)cycle !local terms only:
                         io = pack_indices(isite,iorb)
                         jo = pack_indices(isite,jorb)
                         Jcondition = &
                              (Hij(1,io,jo)/=zero) .AND. &
                              (Nup(jo)==1) .AND. (Nup(io)==0)
                         if (Jcondition) then
                            call c(jo,mup,k1,sg1)
                            call cdg(io,k1,k2,sg2)
                            jup = binary_search(sectorI%H(1)%map,k2)
                            j   = jup + (idw-1)*sectorI%DimUp
                            ed_Eknot = ed_Eknot + &
                                 Hij(1,io,jo)*sg1*sg2*state_cvec(i)*(state_cvec(j))*peso
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !
             !DW - hopping
             idw = Indices(2)
             mdw = sectorI%H(2)%map(idw)
             do iorb=1,Norb
                do jorb=1,Norb
                   do isite=1,Nsites(iorb)
                      do jsite=1,Nsites(jorb)
                         if(isite/=jsite)cycle !local terms only:
                         io = pack_indices(isite,iorb)
                         jo = pack_indices(isite,jorb)
                         Jcondition = &
                              (Hij(Nspin,io,jo)/=zero) .AND. &
                              (Ndw(jo)==1) .AND. (Ndw(io)==0)
                         if (Jcondition) then
                            call c(jo,mdw,k1,sg1)
                            call cdg(io,k1,k2,sg2)
                            jdw = binary_search(sectorI%H(2)%map,k2)
                            j   = iup + (jdw-1)*sectorI%DimUp
                            ed_Eknot = ed_Eknot + &
                                 Hij(Nspin,io,jo)*sg1*sg2*state_cvec(i)*(state_cvec(j))*peso
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !
             !
             !SPIN-EXCHANGE Jx
             if(Jhflag.AND.Jx/=0d0)then
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
                               ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
                               ed_Dse = ed_Dse + sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
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
                               ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
                               ed_Dph = ed_Dph + sg1*sg2*sg3*sg4*state_cvec(i)*state_cvec(j)*peso
                               !
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             endif

             !
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
             do iorb=1,Norb
                do isite=1,Nsites(iorb)          
                   io = pack_indices(isite,iorb)
                   ed_Epot = ed_Epot + Uloc(iorb)*nup(io)*ndw(io)*gs_weight
                enddo
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
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
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      do isite=1,Nsites(iorb)
                         do jsite=1,Nsites(jorb)
                            if(isite/=jsite)cycle !local interaction only:
                            io = pack_indices(isite,iorb)
                            jo = pack_indices(isite,jorb)
                            ed_Epot = ed_Epot + (Ust-Jh)*(nup(io)*nup(jo) + ndw(io)*ndw(jo))*gs_weight
                            ed_Dund = ed_Dund + (nup(io)*nup(jo) + ndw(io)*ndw(jo))*gs_weight
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
                do iorb=1,Norb
                   do isite=1,Nsites(iorb)          
                      io = pack_indices(isite,iorb)
                      ed_Ehartree=ed_Ehartree - 0.5d0*Uloc(iorb)*(nup(io)+ndw(io))*gs_weight + 0.25d0*uloc(iorb)*gs_weight
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
                               ed_Ehartree=ed_Ehartree + 0.25d0*Ust*gs_weight
                               ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(io)+ndw(io)+nup(jo)+ndw(jo))*gs_weight
                               ed_Ehartree=ed_Ehartree + 0.25d0*(Ust-Jh)*gs_weight
                            enddo
                         enddo
                      enddo
                   enddo
                   !
                endif
             endif
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
       call Bcast_MPI(MpiComm,ed_Epot)
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
    endif
#endif
    !
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose==3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
    endif
    !
    if(MPIMASTER)then
       call write_energy_info()
       call write_energy()
    endif
    !
    !
  end subroutine local_energy_impurity









  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ispin,is
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3d0*pi/beta
    do ispin=1,Nspin
       do is=1,Ns
          simp(ispin,is) = dimag(impSmats(ispin,is,is,1)) - &
               wm1*(dimag(impSmats(ispin,is,is,2))-dimag(impSmats(ispin,is,is,1)))/(wm2-wm1)
          zimp(ispin,is)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,is,is,1))/wm1 ))
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin,stride
    unit = free_unit()
    open(unit,file="observables_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+1))//"s2",&
         reg(txtfy(5*Norb+2))//"egs",&
         ((reg(txtfy(5*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((5+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((5+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         ((reg(txtfy((5+2*Norb)*Norb+2+Norb*Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)

    close(unit)
    !

    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend

  subroutine write_energy_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>"
    close(unit)
  end subroutine write_energy_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_all.ed",position='append')
    write(unit,"(90(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs,&
         ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)    
    !
    unit = free_unit()
    open(unit,file="parameters_last.ed")
    write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
    close(unit)
    !
    unit = free_unit()
    open(unit,file="observables_last.ed")
    write(unit,"(90(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs,&
         ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         
    !
    !
  end subroutine write_observables

  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last.ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy


end MODULE ED_OBSERVABLES

















