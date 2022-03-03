  if(Jk/=0d0)then
     do i=MpiIstart,MpiIend
        m  = Hsector%H(1)%map(i)
        ib  = bdecomp(m,2*Ns+Nimp)
        !
        nup = ib(1:Ns)
        ndw = ib(Ns+1:2*Ns)
        sz  = 0.5d0*(Nup-Ndw)
        np  = ib(2*Ns+1:)        !0=DW, 1=UP
        szp = np-0.5d0          !-1/2:DW, 1/2:UP
        !
        ! KONDO COUPLING - non-local (in memory) part of the Kondo coupling
        ! -Jk(s.S) = -Jk*[sx.Sx + sy.Sy]  -2Jk[sz.Sz]
        !            = -Jk*1/2[s^+.S^- + s^-.S^+] - 2Jk*{s^z.S^z}
        !
        !- 2Jk*{s^z.S^z} part of the Kondo coupling:
        do iimp=1,Nimp
           do iorb=1,Norb
              do isite=1,Nsites(iorb)
                 if(isite/=Jkindx(iimp))cycle
                 io = pack_indices(isite,iorb)
                 htmp = htmp - 2d0*Jk*Szp(iimp)*Sz(io)
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
        ! -Jk/2[s^+.S^- + s^-.S^+] - 2Jk*{s^z.S^z}
        ! -Jk/2[Cdg_up.C_dw.S^- + Cdg_dw.C_up.S^+]
        do iimp=1,Nimp
           do iorb=1,Norb
              do isite=1,Nsites(iorb)
                 if(isite/=Jkindx(iimp))cycle
                 io = pack_indices(isite,iorb)
                 !
                 !Cdg_up.C_dw.S^- 
                 Jcondition=((nup(io)==0).AND.(ndw(io)==1).AND.(np(iimp)==1))
                 if(Jcondition)then
                    call Sminus(2*Ns+iimp,m,k1)
                    call c(io+Ns,k1,k2,sg2)
                    call cdg(io,k2,k3,sg3)
                    j=binary_search(Hsector%H(1)%map,k3)
                    htmp = -one*Jk*sg2*sg3
                    !
                    select case(MpiStatus)
                    case (.true.)
                       call sp_insert_element(MpiComm,spH0d,htmp,i,j)
                    case (.false.)
                       call sp_insert_element(spH0d,htmp,i,j)
                    end select
                 endif
                 !
                 !Cdg_dw.C_up.S^+
                 Jcondition=((ndw(io)==0).AND.(nup(io)==1).AND.(np(iimp)==0))
                 if(Jcondition)then
                    call Splus(2*Ns+iimp,m,k1)
                    call c(io,k1,k2,sg2)
                    call cdg(io+Ns,k2,k3,sg3)
                    j=binary_search(Hsector%H(1)%map,k3)
                    htmp = -one*Jk*sg2*sg3
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
