  ! KONDO COUPLING 
  ! -Jk(s.S) = -Jk*[sx.Sx + sy.Sy]  -2Jk[sz.Sz]
  !            = -Jk*1/2[s^+.S^- + s^-.S^+] - 2Jk*{s^z.S^z}
  if(Jk/=0d0)then     
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
     !
     i = j
     hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
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
                 i=binary_search(Hsector%H(1)%map,k3)
                 htmp = -one*Jk*sg2*sg3
                 !
                 if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
              endif
              !
              !Cdg_dw.C_up.S^+
              Jcondition=((ndw(io)==0).AND.(nup(io)==1).AND.(np(iimp)==0))
              if(Jcondition)then
                 call Splus(2*Ns+iimp,m,k1)
                 call c(io,k1,k2,sg2)
                 call cdg(io+Ns,k2,k3,sg3)
                 i=binary_search(Hsector%H(1)%map,k3)
                 htmp = -one*Jk*sg2*sg3
                 !
                 if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
              endif
              !
           enddo
        enddo
     enddo
     !
  endif
