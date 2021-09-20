  !We build the electronic part of the electron-phonon interaction: Sum_iorb g_iorb n_iorb
  !Coupling to combinations of electron density operators
  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     htmp = zero
     do io=1,Ns
        htmp = htmp + g_ph(io)*(nup(io)+ndw(io) - 1.d0)
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0e_eph,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0e_eph,htmp,i,i)
     end select
     !
  enddo
  !


  !Here we build the phononc part of the electron-phonon interaction: (b^+ + b)
  htmp = zero
  do iph=1,DimPh
     i = iph + 1
     if(i <= DimPh) then
        htmp = sqrt(dble(iph))
        call sp_insert_element(spH0ph_eph,htmp,iph,i)
     end if
     i = iph - 1
     if(i>0) then
        htmp = sqrt(dble(iph - 1))
        call sp_insert_element(spH0ph_eph,htmp,iph,i)
     end if
  end do













