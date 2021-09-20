  do jdw=1,MpiQup
     do jup=1,DimDw
        mdw  = Hsector%H(2)%map(jup)
        ndw  = bdecomp(mdw,Ns)
        j    = jup + (jdw-1)*DimDw
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
                 iup = binary_search(Hsector%H(2)%map,k2)
                 idw = jdw
                 i   = iup + (idw-1)*DimDw
                 htmp = Hij(Nspin,io,jo)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*vt(j)
                 !
              endif
           enddo
        enddo
        !
     end do
  enddo
