  ! KONDO COUPLING
  if(any([Jk_z,Jk_xy]/=0d0))then     
     !
     htmp = zero
     !
     do iimp=1,Nimp
        do iorb=1,Norb
           do isite=1,Nsites(iorb)
              if(isite/=Jkindx(iimp))cycle
              io = pack_indices(isite,iorb)
              htmp = -2d0*Jk_z*Sz(io)*Szp(iimp)
           enddo
        enddo
     enddo
     i = j
     hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
     !
     !
     do iimp=1,Nimp
        do iorb=1,Norb
           do isite=1,Nsites(iorb)
              if(isite/=Jkindx(iimp))cycle
              io = pack_indices(isite,iorb)
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
                 i=binary_search(Hsector%H(1)%map,k4)
                 htmp = one*Jk_xy*sg1*sg2*sg3*sg4
                 !
                 if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
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
                   (nup(io)==1).AND.(npup(io)==0) )
              if(Jcondition)then
                 call c(imp_dw,m,k1,sg1)    !d_dw
                 call cdg(io_dw,k1,k2,sg2)  !c^+_dw
                 call c(io_up,k2,k3,sg3)    !c_up
                 call cdg(imp_up,k3,k4,sg4) !d^+_up
                 i=binary_search(Hsector%H(1)%map,k4)
                 htmp = one*Jk_xy*sg1*sg2*sg3*sg4
                 !
                 if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
              endif
              !
           enddo
        enddo
     enddo
     !
     !Kondo exchange (U'-Jk/4): use Jk_z as Jk here
     if(Ust/=0d0)then
        htmp=zero
        do iimp=1,Nimp
           do iorb=1,Norb
              do isite=1,Nsites(iorb)
                 if(isite/=Jkindx(iimp))cycle
                 io = pack_indices(isite,iorb)                 
                 htmp = htmp + Ust*(Nup(io)*Nup(iimp) + Nup(io)*Ndw(iimp) + Ndw(io)*Nup(iimp) + Ndw(io)*Ndw(iimp))
              enddo
           enddo
        enddo
        i = j
        hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
     endif
  endif
