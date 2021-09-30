  !We build the transposed H_non_local here (symmetric)
  !to comply with the MPI decomposition of the matrix.
  !A better MPI handling might be necessary here...
  do i=MpiIstart,MpiIend
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
     ! SPIN-EXCHANGE (S-E) TERMS
     !    S-E: J Cdg_a.up Cdg_b.dw C_a.dw C_b.up
     !    S-E: J (-1)*2 [Cdg_a C_b]_up [Cdg_b C_a]_dw
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
                       jdw=binary_search(Hsector%H(2)%map,k2)
                       call c(jo,mup,k3,sg3)  !UP
                       call cdg(io,k3,k4,sg4) !UP
                       jup=binary_search(Hsector%H(1)%map,k4)
                       htmp = Jx*sg1*sg2*sg3*sg4
                       j = jup + (jdw-1)*DimUp
                       !
                       select case(MpiStatus)
                       case (.true.)
                          call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                       case (.false.)
                          call sp_insert_element(spH0nd,htmp,i,j)
                       end select
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
                       jdw = binary_search(Hsector%H(2)%map,k2)
                       call c(jo,mup,k3,sg3)       !c_jo_up
                       call cdg(io,k3,k4,sg4)      !c^+_io_up
                       jup = binary_search(Hsector%H(1)%map,k4)
                       htmp = Jp*sg1*sg2*sg3*sg4
                       j = jup + (jdw-1)*DimUp
                       !
                       select case(MpiStatus)
                       case (.true.)
                          call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                       case (.false.)
                          call sp_insert_element(spH0nd,htmp,i,j)
                       end select
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif

     ! KONDO COUPLING - non-local (in memory) part of the Kondo coupling
     ! Sa.Sb--> 1/2[Sa^+ Sb^- + Sa^-Sb^+] + {Sa^z Sb^z--> local}
     ! Jk: cdg_a.up c_a.dw cdg_b.dw c_b.up + cdg_a.dw c_a.up cdg_b.up c_b.dw
     ! Jk: -(-1)**[2+1]  [cdg_a.up c_b.up] [cdg_b.dw c_a.dw] +
     !   : -(-1)**[2+1] [cdg_a.dw c_b.dw] [cdg_b.up c_a.up ]
     !   : Jk*{[cdg_a c_b]_up [cdg_b c_a]_dw + [cdg_b c_a]_up [cdg_a c_b]_dw }
     if(Jhflag.AND.Jk/=0d0)then
        do iorb=1,Norb
           do jorb=1,Norb
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
                       jdw = binary_search(Hsector%H(2)%map,k2)
                       call c(jo,mup,k3,sg3)       !c_jo.up
                       call cdg(io,k3,k4,sg4)      !c^+_io.up
                       jup = binary_search(Hsector%H(1)%map,k4)
                       htmp = Jk*sg1*sg2*sg3*sg4
                       j = jup + (jdw-1)*DimUp
                       !
                       select case(MpiStatus)
                       case (.true.)
                          call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                       case (.false.)
                          call sp_insert_element(spH0nd,htmp,i,j)
                       end select
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
                       jdw = binary_search(Hsector%H(2)%map,k2)
                       call c(io,mup,k3,sg3)       !c_io.up
                       call cdg(jo,k3,k4,sg4)      !c^+_jo.up
                       jup = binary_search(Hsector%H(1)%map,k4)
                       htmp = Jk*sg1*sg2*sg3*sg4
                       j = jup + (jdw-1)*DimUp
                       !
                       select case(MpiStatus)
                       case (.true.)
                          call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                       case (.false.)
                          call sp_insert_element(spH0nd,htmp,i,j)
                       end select
                    endif
                    !
                 enddo
              enddo
           enddo
        enddo
     endif

  enddo
