  do j=1,Nloc
     jup = iup_index(j,DimUp)
     jdw = idw_index(j,DimUp)
     !
     mup = Hsector%H(1)%map(jup)
     mdw = Hsector%H(2)%map(jdw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
     !    S-E: J c^+_io_up c^+_jorb_dw c_io_dw c_jorb_up  (i.ne.j) 
     !    S-E: J c^+_{io} c^+_{jorb+Ns} c_{io+Ns} c_{jorb}
     if(Jhflag.AND.Jx/=0d0)then
        do iorb=1,Norb
           do jorb=1,Norb
              do isite=1,Nsites(iorb)
                 do jsite=1,Nsites(jorb)
                    if(isite/=jsite)cycle !local interaction only:
                    io = pack_indices(iorb,isite)
                    jo = pack_indices(jorb,isite)
                    Jcondition=(&
                         (io/=jo).AND.&
                         (nup(jo)==1).AND.& 
                         (ndw(io)==1).AND.&
                         (ndw(jo)==0).AND.&
                         (nup(io)==0))
                    if(Jcondition)then
                       call c(io,mdw,k1,sg1)  !DW
                       call cdg(jo,k1,k2,sg2) !DW
                       idw=binary_search(Hsector%H(2)%map,k2)
                       call c(jo,mup,k3,sg3)  !UP
                       call cdg(io,k3,k4,sg4) !UP
                       iup=binary_search(Hsector%H(1)%map,k4)
                       htmp = Jx*sg1*sg2*sg3*sg4
                       i = iup + (idw-1)*dimup
                       !
                       Hv(i) = Hv(i) + htmp*vin(j)
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif
     !
     ! PAIR-HOPPING (P-H) TERMS
     !    P-H: J c^+_io_up c^+_io_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
     !    P-H: J c^+_{io}  c^+_{io+Ns} c_{jorb+Ns} c_{jorb}
     if(Jhflag.AND.Jp/=0d0)then
        do iorb=1,Norb
           do jorb=1,Norb
              do isite=1,Nsites(iorb)
                 do jsite=1,Nsites(jorb)
                    if(isite/=jsite)cycle !local interaction only:
                    io = pack_indices(iorb,isite)
                    jo = pack_indices(jorb,isite)
                    Jcondition=(&
                         (nup(jo)==1).AND.&
                         (ndw(jo)==1).AND.&
                         (ndw(io)==0).AND.&
                         (nup(io)==0))
                    if(Jcondition)then
                       call c(jo,mdw,k1,sg1)       !c_jorb_dw
                       call cdg(io,k1,k2,sg2)      !c^+_iorb_dw
                       idw = binary_search(Hsector%H(2)%map,k2)
                       call c(jo,mup,k3,sg3)       !c_jorb_up
                       call cdg(io,k3,k4,sg4)      !c^+_iorb_up
                       iup = binary_search(Hsector%H(1)%map,k4)
                       htmp = Jp*sg1*sg2*sg3*sg4
                       i = iup + (idw-1)*DimUp
                       !
                       Hv(i) = Hv(i) + htmp*vin(j)
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif
     !
  enddo
