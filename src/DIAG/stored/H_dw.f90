  do jdw=1,DimDw
     mdw  = Hsector%H(2)%map(jdw)
     Ndw  = bdecomp(mdw,Ns)
     !
     htmp = zero
     !
     !> Hdw: Off-diagonal elements, i.e. non-local part. 
     !remark: io=jo cant have simultaneously n=0 and n=1 (Jcondition)
     !        so diagonal element (in H_local) are neglected
     do io=1,Ns
        do jo=1,Ns
           Jcondition = &
                (Hij(Nspin,io,jo)/=zero) .AND. &
                (Ndw(jo)==1) .AND. (Ndw(io)==0)
           if (Jcondition) then
              call c(jo,mdw,k1,sg1)
              call cdg(io,k1,k2,sg2)
              idw = binary_search(Hsector%H(2)%map,k2)
              htmp = conjg(Hij(Nspin,io,jo))*sg1*sg2
              !
              call sp_insert_element(spH0dws(1),htmp,idw,jdw)
              !
           endif
        enddo
     enddo
     !
  enddo

