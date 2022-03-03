  do j=1,Nloc
     jup = iup_index(j+mpiIshift,DimUp)
     jdw = idw_index(j+mpiIshift,DimUp)
     !
     mup = Hsector%H(1)%map(jup)
     mdw = Hsector%H(2)%map(jdw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     ! SPIN-EXCHANGE (S-E) TERMS
     !    S-E: J Cdg_a.up Cdg_b.dw C_a.dw C_b.up
     !    S-E: J (-1)*2 [Cdg_a C_b]_up [Cdg_b C_a]_dw
     if(Jx/=0d0)then
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
                       idw=binary_search(Hsector%H(2)%map,k2)
                       call c(jo,mup,k3,sg3)  !UP
                       call cdg(io,k3,k4,sg4) !UP
                       iup=binary_search(Hsector%H(1)%map,k4)
                       htmp = one*Jx*sg1*sg2*sg3*sg4
                       i = iup + (idw-1)*DimUp
                       !
                       Hv(j) = Hv(j) + htmp*vt(i)
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif
     !
     ! PAIR-HOPPING (P-H) TERMS
     !    P-H: J Cdg_a.up Cdg_a.dw   C_b.dw   C_b.up
     !    P-H: J (-1)**2 [Cdg_a C_b]_up [Cdg_a C_b]_dw  
     if(Jp/=0d0)then
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
                       call c(jo,mdw,k1,sg1)       !c_jorb_dw
                       call cdg(io,k1,k2,sg2)      !c^+_iorb_dw
                       idw = binary_search(Hsector%H(2)%map,k2)
                       call c(jo,mup,k3,sg3)       !c_jorb_up
                       call cdg(io,k3,k4,sg4)      !c^+_iorb_up
                       iup = binary_search(Hsector%H(1)%map,k4)
                       htmp = one*Jp*sg1*sg2*sg3*sg4
                       i = iup + (idw-1)*DimUp
                       !
                       Hv(j) = Hv(j) + htmp*vt(i)
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif
     !
  enddo
