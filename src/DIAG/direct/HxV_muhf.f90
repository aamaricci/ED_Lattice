  do i=1,Nloc
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     Nup = bdecomp(mup,Ns)
     Ndw = bdecomp(mdw,Ns)
     !
     !
     htmp = zero
     do io=1,Ns
        htmp = htmp - xmu*(Nup(io)+Ndw(io))
     enddo
     !
     if(hfmode)then
        if(any(Uloc/=0d0))then
           do iorb=1,Norb
              do isite=1,Nsites(iorb) 
                 io = pack_indices(isite,iorb)
                 htmp = htmp-0.5d0*Uloc(iorb)*(Nup(io)+Ndw(io)) !+ 0.25d0*Uloc(iorb)
              enddo
           enddo
        endif
        !
        if(Norb>1)then
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 do isite=1,Nsites(iorb)
                    do jsite=1,Nsites(jorb)
                       if(isite/=jsite)cycle !local interaction only:
                       io = pack_indices(isite,iorb)
                       jo = pack_indices(isite,jorb)
                       htmp=htmp - 0.5d0*Ust*(Nup(io)+Ndw(io)+Nup(jo)+Ndw(jo))
                       htmp=htmp - 0.5d0*(Ust-Jh)*(Nup(io)+Ndw(io)+Nup(jo)+Ndw(jo))
                       ! htmp=htmp + 0.25d0*Ust
                       ! htmp=htmp + 0.25d0*(Ust-Jh)
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



