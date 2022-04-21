  htmp = zero
  !

  do io=1,eNs
     do jo=1,eNs
        !UP electrons
        Jcondition = &
             (Hij(1,io,jo)/=zero) .AND. &
             (nup(jo)==1) .AND. (nup(io)==0)
        if (Jcondition) then
           call c(ket_site_index(jo,1),m,k1,sg1)
           call cdg(ket_site_index(io,1),k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(1,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        !DW electrons
        Jcondition = &
             (Hij(Nspin,io,jo)/=zero) .AND. &
             (ndw(jo)==1) .AND. (ndw(io)==0)
        if (Jcondition) then
           call c(ket_site_index(jo,2),m,k1,sg1)
           call cdg(ket_site_index(io,2),k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(Nspin,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo


  do io=1,iNs
     do jo=1,iNs
        !UP impurity
        Jcondition = (npup(jo)==1) .AND. (npup(io)==0)
        if (Jcondition) then
           call c(ket_imp_index(jo,1),m,k1,sg1)
           call cdg(ket_imp_index(io,1),k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        !DW impurity
        Jcondition = (npdw(jo)==1) .AND. (npdw(io)==0)
        if (Jcondition) then
           call c(ket_imp_index(jo,2),m,k1,sg1)
           call cdg(ket_imp_index(io,2),k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo

