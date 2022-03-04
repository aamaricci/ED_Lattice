  if(Jk/=0d0)then
     do i=MpiIstart,MpiIend
        m  = Hsector%H(1)%map(i)
        ib  = bdecomp(m,2*Ns_imp)
        !
        Nup = ib(1:Ns)
        Ndw = ib(Ns+1:2*Ns)
        NpUp= ib(2*Ns+1:2*Ns+Nimp)
        NpDw= ib(2*Ns+Nimp+1:2*Ns+2*Nimp)
        Sz  = 0.5d0*(Nup-Ndw)
        Szp = 0.5d0*(NpUp-NpDw)
        !
        ! KONDO COUPLING
        ! -Jk(s.S) = -Jk*[Sx.Simpx + Sy.Simpy]  -2Jk[Sz.Simpz]
        !            = -Jk*1/2[S^+.Simp^- + S^-.Simp^+] - 2Jk*{S^z.Simp^z}
        !
        !- 2Jk*{S^z.Simp^z} part of the Kondo coupling:
        do iimp=1,Nimp
           do iorb=1,Norb
              do isite=1,Nsites(iorb)
                 if(isite/=Jkindx(iimp))cycle
                 io = pack_indices(isite,iorb)
                 htmp = htmp - 2d0*Jk*Sz(io)*Szp(iimp)
              enddo
           enddo
        enddo
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0d,htmp,i,i)
        case (.false.)
           call sp_insert_element(spH0d,htmp,i,i)
        end select
        !
        !
        !  -Jk[S^+ Simp^- + S^-Simp^+]
        !  -Jk[c^+_up c_dw d^+_dw d_up + c^+_dw c_up d^+_up d_dw]
        !  -Jk*(-1)**[2+1] {[c^+_up.d_up] [d^+_dw.c_dw] + [c^+_dw.d_dw] [d^+_up.c_up]}
        !   Jk*{[c^+.d]_up [d^+.c]_dw + [d^+.c]_up [c^+.d]_dw }
        do iimp=1,Nimp
           do iorb=1,Norb
              do isite=1,Nsites(iorb)
                 if(isite/=Jkindx(iimp))cycle
                 io    = pack_indices(isite,iorb)
                 io_up = io
                 io_dw = io + Ns
                 imp_up= 2*Ns + iimp
                 imp_dw= 2*Ns + iimp + Nimp
                 ![c^+.d]_up [d^+.c]_dw
                 Jcondition=(&
                      (ndw(io)==1).AND.(npdw(iimp)==0).AND.&
                      (npup(iimp)==1).AND.(nup(io)==0) )
                 if(Jcondition)then
                    call c(io_dw,m,k1,sg1)     !c_dw
                    call cdg(imp_dw,k1,k2,sg2) !d^+_dw
                    call c(imp_up,k2,k3,sg3)   !d_up
                    call cdg(io_up,k3,k4,sg4)  !c^+_up
                    j=binary_search(Hsector%H(1)%map,k4)
                    htmp = one*Jk*sg1*sg2*sg3*sg4
                    !
                    select case(MpiStatus)
                    case (.true.)
                       call sp_insert_element(MpiComm,spH0d,htmp,i,j)
                    case (.false.)
                       call sp_insert_element(spH0d,htmp,i,j)
                    end select
                 endif
                 !
                 ![d^+.c]_up [c^+.d]_dw 
                 io    = pack_indices(isite,iorb)
                 io_up = io
                 io_dw = io + Ns
                 imp_up= 2*Ns + iimp
                 imp_dw= 2*Ns + iimp + Nimp
                 Jcondition=(&
                      (npdw(iimp)==1).AND.(ndw(io)==0).AND.&
                      (nup(io)==1).AND.(npup(iimp)==0) )
                 if(Jcondition)then
                    call c(imp_dw,m,k1,sg1)    !d_dw
                    call cdg(io_dw,k1,k2,sg2)  !c^+_dw
                    call c(io_up,k2,k3,sg3)    !c_up
                    call cdg(imp_up,k3,k4,sg4) !d^+_up
                    j=binary_search(Hsector%H(1)%map,k4)
                    htmp = one*Jk*sg1*sg2*sg3*sg4
                    !
                    select case(MpiStatus)
                    case (.true.)
                       call sp_insert_element(MpiComm,spH0d,htmp,i,j)
                    case (.false.)
                       call sp_insert_element(spH0d,htmp,i,j)
                    end select
                 endif
                 !
              enddo
           enddo
        enddo
        !
     enddo
  endif
