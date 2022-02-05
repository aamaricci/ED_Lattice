  do i=1,Nloc
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     Nup = bdecomp(mup,Ns)
     Ndw = bdecomp(mdw,Ns)
     Sz  = 0.5d0*(Nup-Ndw)
     !
     !
     !> HxV_imp: Diagonal Elements, i.e. local part
     htmp = zero
     do io=1,Ns
        htmp = htmp + Hdiag(1,io)*Nup(io) + Hdiag(Nspin,io)*Ndw(io)
     enddo
     !
     !
     !> H_Int: Kanamori interaction part.
     !  = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     if(any(Uloc/=0d0))then
        do iorb=1,Norb
           do isite=1,Nsites(iorb)          
              io = pack_indices(isite,iorb)
              htmp = htmp + Uloc(iorb)*Nup(io)*Ndw(io)
           enddo
        enddo
     endif
     !
     !density-density interaction: different orbitals, opposite spins:
     ! =  (Uprime==Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           do isite=1,Nsites(iorb)
              do jsite=1,Nsites(jorb)
                 if(isite/=jsite)cycle
                 io = pack_indices(isite,iorb)
                 jo = pack_indices(isite,jorb)
                 htmp = htmp + Ust*(Nup(io)*Ndw(jo) + Nup(jo)*Ndw(io))
              enddo
           enddo
        enddo
     enddo
     !
     !density-density interaction: different orbitals, parallel spins
     ! = \sum_{i<j} (Usecond==Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           do isite=1,Nsites(iorb)
              do jsite=1,Nsites(jorb)
                 if(isite/=jsite)cycle
                 io = pack_indices(isite,iorb)
                 jo = pack_indices(isite,jorb)
                 htmp = htmp + (Ust-Jh)*(Nup(io)*Nup(jo) + Ndw(io)*Ndw(jo))
              enddo
           enddo
        enddo
     enddo
     !
     !Sz_a.Sz_b part of the Kondo coupling:
     if(Jhflag.AND.Jk/=0d0)then
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              do isite=1,Nsites(iorb)
                 do jsite=1,Nsites(jorb)
                    ! if(isite/=jsite)cycle !local interaction only:
                    ! io = pack_indices(isite,iorb)
                    ! jo = pack_indices(isite,jorb)
                    if(jsite/=Jkindx(isite))cycle
                    io = pack_indices(isite,iorb)
                    jo = pack_indices(jsite,jorb)
                    htmp = htmp - 2*Jk*Sz(io)*Sz(jo)
                 enddo
              enddo
           enddo
        enddo
     endif
     !
     !
     hv(i) = hv(i) + htmp*vin(i)
     !
  enddo



