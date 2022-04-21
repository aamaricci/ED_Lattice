  ! KONDO COUPLING
  if(any([Jk_z,Jk_xy]/=0d0))then     
     !
     htmp = zero
     !
     do iimp=1,iNs
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
     do iimp=1,iNs
        do iorb=1,Norb
           do isite=1,Nsites(iorb)
              if(isite/=Jkindx(iimp))cycle
              io = pack_indices(isite,iorb)
              io_up  = ket_site_index(io,1)
              io_dw  = ket_site_index(io,2)
              imp_up = ket_imp_index(iimp,1)
              imp_dw = ket_imp_index(iimp,2)
              ! io_up = io
              ! io_dw = io + eNs
              ! imp_up= 2*Ns + iimp
              ! imp_dw= 2*Ns + iimp + iNs
              !c^+_up d^+_dw c_dw  d_up
              Jcondition=(&
                   (npup(iimp)==1).AND.&
                   (ndw(io)   ==1).AND.&
                   (npdw(iimp)==0).AND.&
                   (nup(io)   ==0) )
              if(Jcondition)then
                 call c(imp_up,  m, k1,sg1)  !d_up
                 call c(io_dw,   k1,k2,sg2)  !c_dw
                 call cdg(imp_dw,k2,k3,sg3)  !d^+_dw
                 call cdg(io_up, k3,k4,sg4)  !c^+_up
                 i=binary_search(Hsector%H(1)%map,k4)
                 htmp = one*Jk_xy*sg1*sg2*sg3*sg4
                 !
                 if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
              endif
              !
              ! c^+_dw d^+_up c_up  d_dw                 
              Jcondition=(&
                   (npdw(iimp)==1).AND.&
                   (nup(io)   ==1).AND.&
                   (npup(iimp)==0).AND.& 
                   (ndw(io)   ==0) )
              if(Jcondition)then
                 call c(imp_dw,m,k1,sg1)    !d_dw
                 call c(io_up,k1,k2,sg2)    !c_up
                 call cdg(imp_up,k2,k3,sg3) !d^+_up
                 call cdg(io_dw,k3,k4,sg4)  !c^+_dw
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
        do iimp=1,iNs
           do iorb=1,Norb
              do isite=1,Nsites(iorb)
                 if(isite/=Jkindx(iimp))cycle
                 io = pack_indices(isite,iorb)                 
                 htmp = htmp + Ust*(Nup(io)*NpUp(iimp) + Nup(io)*NpDw(iimp) + Ndw(io)*NpUp(iimp) + Ndw(io)*NpDw(iimp))
              enddo
           enddo
        enddo
        i = j
        hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
     endif
  endif
