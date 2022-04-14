  !We build the transposed H here, so we take conjg part
  !to comply with the MPI decomposition of the matrix.
  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
     ib  = bdecomp(m,2*Ns_imp)
     Nup = ib(1:Ns)
     Ndw = ib(Ns+1:2*Ns)
     NpUp= ib(2*Ns+1:2*Ns+iNs)
     NpDw= ib(2*Ns+iNs+1:2*Ns+2*iNs)
     !
     htmp = zero


     do io=1,Ns
        do jo=1,Ns
           !UP electrons
           Jcondition = (Hij(1,io,jo)/=zero) .AND. (nup(jo)==1) .AND. (nup(io)==0)
           if (Jcondition) then
              call c(jo,m,k1,sg1)
              call cdg(io,k1,k2,sg2)
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
              call c(jo+Ns,m,k1,sg1)
              call cdg(io+Ns,k1,k2,sg2)
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


     do io=1,iNs
        do jo=1,iNs
           !UP impurity
           Jcondition =  (t_imp/=0d0) .AND. (npup(jo)==1) .AND. (npup(io)==0)
           if (Jcondition) then
              call c(2*Ns + jo,m,k1,sg1)
              call cdg(2*Ns + io,k1,k2,sg2)
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
           !DW impurity
           Jcondition = (t_imp/=0d0) .AND. (npdw(jo)==1) .AND. (npdw(io)==0)
           if (Jcondition) then
              call c(2*Ns + iNs + jo,m,k1,sg1)
              call cdg(2*Ns + iNs + io,k1,k2,sg2)
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

  enddo


