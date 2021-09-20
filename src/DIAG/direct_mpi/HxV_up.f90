  do iph=1,DimPh
     do jdw=1,MpiQdw
        do jup=1,DimUp
           mup  = Hsector%H(1)%map(jup)
           nup  = bdecomp(mup,Ns)
           j    = jup + (jdw-1)*dimUp + (iph-1)*DimUp*MpiQdw
           !
           !> Hup: Off-diagonal elements, i.e. non-local part. 
           !remark: io=jo cant have simultaneously n=0 and n=1 (Jcondition)
           !        so diagonal element (in H_local) are neglected
           do io=1,Ns
              do jo=1,Ns
                 Jcondition = &
                      (Hij(1,io,jo)/=zero) .AND. &
                      (Nup(jo)==1) .AND. (Nup(io)==0)
                 if (Jcondition) then
                    call c(jo,mup,k1,sg1)
                    call cdg(io,k1,k2,sg2)
                    iup  = binary_search(Hsector%H(1)%map,k2)
                    idw  = jdw
                    i    = iup + (idw-1)*dimUp + (iph-1)*DimUp*MpiQdw
                    htmp = Hij(1,io,jo)*sg1*sg2
                    !
                    Hv(i) = Hv(i) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
           !
           !
        enddo
     enddo
  enddo


