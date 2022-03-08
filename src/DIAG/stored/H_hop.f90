  !We build the transposed H here, so we take conjg part
  !to comply with the MPI decomposition of the matrix.
  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
     ib  = bdecomp(m,2*Ns_imp)
     Nup = ib(1:Ns)
     Ndw = ib(Ns+1:2*Ns)
     !
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
        enddo
     enddo
     !
     !
     !DW electrons
     do io=1,Ns
        do jo=1,Ns
           Jcondition = &
                (Hij(Nspin,io,jo)/=zero) .AND. &
                (ndw(jo)==1) .AND. (ndw(io)==0)
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
     !
  enddo


