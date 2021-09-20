  !phonon coupling to electron density operators
  do i=1,Nloc
     i_el = mod(i-1,DimUp*MpiQdw) + 1
     iph = (i-1)/(DimUp*MpiQdw) + 1
     !
     iup = iup_index(i_el+mpiIshift,DimUp)
     idw = idw_index(i_el+mpiIshift,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     htmp=zero
     do io=1,Ns
        htmp = htmp + g_ph(io)*(nup(io)+ndw(io) - 1.d0)
     enddo
     Hv(i) = Hv(i) + htmp*vin(i)
     !
     do jj = 1,DimPh
        if(jj .eq. iph+1) then
           j = i_el + (jj-1)*DimUp*MpiQdw 
           Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
        endif
        !
        if(jj .eq. iph-1) then
           j = i_el + (jj-1)*DimUp*MpiQdw
           Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
        endif
     enddo
  enddo

