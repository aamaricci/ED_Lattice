module ED_MAIN
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace,delete_eigenspace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_CHI_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  implicit none
  private
  !
  public :: ed_init_solver
  public :: ed_solve


  !Boolean to select MPI mode for the inequivalent sites solver:
  !T: solve each site SERIALLY, using MPI lanczos, 
  !F: solve each site PARALLEL, using Serial Lanczos   
  logical :: mpi_lanc_=.true.



contains





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate Memory and Initialize ED -+!
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver()
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    if(mpi_lanc_)call ed_set_MpiComm()
    !
    if(MpiMaster)write(LOGfile,"(A)")"INIT ED SOLVER"
    !
    !Init ED Structure & memory
    call init_ed_structure()
    !
    !Check Hmatrix is allocated:
    if(.not.Hmatrix%status)stop "ED_INIT_SOLVER ERROR: Hmatrix is not allocated"
    !
    call setup_global
    !
    call ed_del_MpiComm()
    !
  end subroutine ed_init_solver



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve()
    !
    !SET THE LOCAL MPI COMMUNICATOR
    !Only if mpi_lanc_ is TRUE: solve using MPI Lanczos  
    if(mpi_lanc_)call ed_set_MpiComm()
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    if(.not.Hmatrix%status)stop "ED_INIT_SOLVER ERROR: Hmatrix is not allocated"
    call Hij_write(unit=LOGfile)
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()
    call observables_impurity()
    call internal_energy_impurity()
    if(gf_flag)call buildgf_impurity()
    if(chi_flag)call buildchi_impurity()
    !
    select case(ed_method)
    case default
       call es_delete_espace(state_list)
    case("lapack","full")
       call delete_eigenspace()
    end select
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    nullify(spHtimesV_p)
    if(MpiMaster)write(Logfile,"(A)")""
  end subroutine ed_solve

end module ED_MAIN







