  do jup=1,DimUp
     do jdw=1,DimDw
        mdw  = Hsector%H(2)%map(jdw)
        ndw  = bdecomp(mdw,Ns)
        j    = jup + (jdw-1)*dimUp
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
                 iup = jup
                 i   = iup + (idw-1)*DimUp
                 htmp = Hij(Nspin,io,jo)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
     enddo
  enddo
