MODULE EDLAT
  USE ED_INPUT_VARS, only: &
       ed_read_input , &
       Nsites        , &
       Norb          , &
       Nspin         , &
       Uloc          , &
       Ust           , &
       Jh            , &
       Jx            , &
       Jp            , &
       xmu           , &
       temp          , &
       eps           , &
       wini          , &
       wfin          , &
       sb_field      , &
       Lmats         , &
       Lreal         , &
       LOGfile   

  USE ED_GRAPH_MATRIX, only: &
       ed_Hij_init     => Hij_init     , &
       ed_Hij_add_link => Hij_add_link , &
       ed_Hij_read     => Hij_read     , &
       ed_Hij_get      => Hij_get      , &
       ed_Hij_local    => Hij_local    , &
       ed_Hij_info     => Hij_info     , &
       ed_Hij_write    => Hij_write

  USE ED_DENSITY_MATRIX, only: &
       ed_get_density_matrix

  USE ED_IO, only: &
       ed_get_gimp_matsubara  , &
       ed_get_gimp_realaxis   , &
       ed_get_dens            , &
       ed_get_mag             , &
       ed_get_docc            , &
       ed_get_energy          , &
       ed_get_doubles


  USE ED_MAIN, only: &
       ed_init_solver , &
       ed_solve


END MODULE EDLAT

