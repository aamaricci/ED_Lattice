  htmp = zero
  !
  !UP electrons
  do io=1,Ns
     do jo=1,Ns
        Jcondition = &
             (Hij(1,io,jo)/=zero) .AND. &
             (nup(jo)==1) .AND. (nup(io)==0)
        if (Jcondition) then
           call c(jo,m,k1,sg1)
           call cdg(io,k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(1,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo

  !DW electrons
  do io=1,Ns
     do jo=1,Ns
        Jcondition = &
             (Hij(Nspin,io,jo)/=zero) .AND. &
             (ndw(jo)==1) .AND. (ndw(io)==0)
        if (Jcondition) then
           call c(jo+Ns,m,k1,sg1)
           call cdg(io+Ns,k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(Nspin,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo




  !UP impurity
  do io=1,iNs
     do jo=1,iNs
        Jcondition = (npup(jo)==1) .AND. (npup(io)==0)
        if (Jcondition) then
           call c(2*Ns + jo,m,k1,sg1)
           call cdg(2*Ns + io,k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo
  !
  !
  !DW impurity
  do io=1,iNs
     do jo=1,iNs
        Jcondition = (npdw(jo)==1) .AND. (npdw(io)==0)
        if (Jcondition) then
           call c(2*Ns + iNs + jo,m,k1,sg1)
           call cdg(2*Ns + iNs + io,k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo

