  !We build the transposed H here, so we take conjg part
  !to comply with the MPI decomposition of the matrix.
  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
     ib  = bdecomp(m,2*Ns)
     Nup = ib(1:eNs)
     Ndw = ib(eNs+1:2*eNs)
     NpUp= ib(2*eNs+1:2*eNs+iNs)
     NpDw= ib(2*eNs+iNs+1:2*eNs+2*iNs)
     !
     htmp = zero
     !
     do io=1,eNs
        do jo=1,eNs
           !UP electrons
           Jcondition = (Hij(1,io,jo)/=zero) .AND. (nup(jo)==1) .AND. (nup(io)==0)
           if (Jcondition) then
              call c(ket_site_index(jo,1),m,k1,sg1)
              call cdg(ket_site_index(io,1),k1,k2,sg2)
              j    = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(Hij(1,io,jo))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0d,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0d,htmp,i,j)
              end select
           endif
           !DW electrons
           Jcondition = (Hij(Nspin,io,jo)/=zero) .AND. (ndw(jo)==1) .AND. (ndw(io)==0)
           if (Jcondition) then
              call c(ket_site_index(jo,2),m,k1,sg1)
              call cdg(ket_site_index(io,2),k1,k2,sg2)
              j    = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(Hij(Nspin,io,jo))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0d,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0d,htmp,i,j)
              end select
           endif
        enddo
     enddo


     do io=1,iNs-1
        !UP impurity: C^+_i C_(i+1) i-->i+1
        Jcondition =  (t_imp/=0d0) .AND. (npup(io+1)==1) .AND. (npup(io)==0)
        if (Jcondition) then
           call c(ket_imp_index(io+1,1),m,k1,sg1)
           call cdg(ket_imp_index(io,1),k1,k2,sg2)
           j    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           select case(MpiStatus)
           case (.true.)
              call sp_insert_element(MpiComm,spH0d,htmp,i,j)
           case (.false.)
              call sp_insert_element(spH0d,htmp,i,j)
           end select
        endif
        !UP impurity: C^+_(i+1) C_i i<--i+1
        Jcondition =  (t_imp/=0d0) .AND. (npup(io)==1) .AND. (npup(io+1)==0)
        if (Jcondition) then
           call c(ket_imp_index(io,1),m,k1,sg1)
           call cdg(ket_imp_index(io+1,1),k1,k2,sg2)
           j    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           select case(MpiStatus)
           case (.true.)
              call sp_insert_element(MpiComm,spH0d,htmp,i,j)
           case (.false.)
              call sp_insert_element(spH0d,htmp,i,j)
           end select
        endif
        !
        !DW impurity C^+_i C_(i+1) i-->i+1
        Jcondition = (t_imp/=0d0) .AND. (npdw(io+1)==1) .AND. (npdw(io)==0)
        if (Jcondition) then
           call c(ket_imp_index(io+1,2),m,k1,sg1)
           call cdg(ket_imp_index(io,2),k1,k2,sg2)
           j    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           select case(MpiStatus)
           case (.true.)
              call sp_insert_element(MpiComm,spH0d,htmp,i,j)
           case (.false.)
              call sp_insert_element(spH0d,htmp,i,j)
           end select
        endif
        !DW impurity C^+_(i+1) C_i i<--i+1
        Jcondition = (t_imp/=0d0) .AND. (npdw(io)==1) .AND. (npdw(io+1)==0)
        if (Jcondition) then
           call c(ket_imp_index(io,2),m,k1,sg1)
           call cdg(ket_imp_index(io+1,2),k1,k2,sg2)
           j    = binary_search(Hsector%H(1)%map,k2)
           htmp = t_imp*sg1*sg2
           !
           select case(MpiStatus)
           case (.true.)
              call sp_insert_element(MpiComm,spH0d,htmp,i,j)
           case (.false.)
              call sp_insert_element(spH0d,htmp,i,j)
           end select
        endif
     enddo

  enddo


