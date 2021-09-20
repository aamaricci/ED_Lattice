  do i=1,Nloc
     !
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     !
     !> HxV_imp: Diagonal Elements, i.e. local part
     htmp = zero
     do io=1,Ns
        htmp = htmp + Hdiag(1,io)*nup(io) + Hdiag(Nspin,io)*ndw(io)
        htmp = htmp - xmu*(Nup(io)+Ndw(io))
     enddo
     !
     !
     !> H_Int: Kanamori interaction part. non-local S-E and P-H terms commented below.
     !
     ! density-density interaction: same orbital, opposite spins:
     !  = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     do iorb=1,Norb
        do isite=1,Nsites(iorb)          
           io = pack_indices(iorb,isite)
           htmp = htmp + Uloc(iorb)*Nup(io)*Ndw(io)
        enddo
     enddo
     !
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins:
        ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              do isite=1,Nsites(iorb)
                 do jsite=1,Nsites(jorb)
                    if(isite/=jsite)cycle !local interaction only:
                    io = pack_indices(iorb,isite)
                    jo = pack_indices(jorb,isite)
                    htmp = htmp + Ust*(Nup(io)*Ndw(jo) + Nup(jo)*Ndw(io))
                 enddo
              enddo
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              do isite=1,Nsites(iorb)
                 do jsite=1,Nsites(jorb)
                    if(isite/=jsite)cycle !local interaction only:
                    io = pack_indices(iorb,isite)
                    jo = pack_indices(jorb,isite)
                    htmp = htmp + (Ust-Jh)*(Nup(io)*Nup(jo) + Ndw(io)*Ndw(jo))
                 enddo
              enddo
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        do iorb=1,Norb
           do isite=1,Nsites(iorb)          
              io = pack_indices(iorb,isite)
              htmp = htmp - 0.5d0*Uloc(iorb)*(nup(io)+ndw(io)) + 0.25d0*Uloc(iorb)
           enddo
        enddo
        if(Norb>1)then
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 do isite=1,Nsites(iorb)
                    do jsite=1,Nsites(jorb)
                       if(isite/=jsite)cycle !local interaction only:
                       io = pack_indices(iorb,isite)
                       jo = pack_indices(jorb,isite)
                       htmp=htmp-0.5d0*Ust*(nup(io)+ndw(io)+nup(jo)+ndw(jo))+0.25d0*Ust
                       htmp=htmp-0.5d0*(Ust-Jh)*(nup(io)+ndw(io)+nup(jo)+ndw(jo))+0.25d0*(Ust-Jh)
                    enddo
                 enddo
              enddo
           enddo
        endif
     endif
     !
     !
     hv(i) = hv(i) + htmp*vin(i)
     !
  enddo



