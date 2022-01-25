module ED_MAIN
  USE SF_IOTOOLS, only: str,reg,free_unit,file_length
  USE SF_TIMER,only: start_timer,stop_timer
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace,delete_eigenspace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_CHI_FUNCTIONS
  USE ED_OC_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  implicit none
  private
  !
  public :: ed_init_solver
  public :: ed_solve



  integer             :: Bstep
  logical             :: Bbool
  real(8),allocatable :: Blist(:)
contains





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate Memory and Initialize ED -+!
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver()
    integer :: unit,ibeta
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm()
    !
    if(MpiMaster)write(LOGfile,"(A)")"INIT ED SOLVER"
    !
    !Init ED Structure & memory
    call init_ed_structure()
    !
    inquire(file=trim(Bfile)//".restart",exist=Bbool)
    if(Bbool)then
       write(LOGfile,"(A)")'Reading temperature list from file'//trim(Bfile)//".restart"
       Bstep = file_length(trim(Bfile)//".restart")
       open(free_unit(unit),file=trim(Bfile)//".restart")
       allocate(Blist(Bstep))
       do ibeta=1,Bstep
          read(unit,*)Blist(ibeta)
       enddo
       close(unit)
    else
       Bstep = 1
       allocate(Blist(Bstep))
       Blist = beta
    endif
    !
    !Check Hmatrix is allocated:
    if(.not.Hmatrix%status)stop "ED_INIT_SOLVER ERROR: Hmatrix is not allocated"
    !
    call setup_global
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
  end subroutine ed_init_solver



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve()
    integer :: ibeta
    !
    !SET THE LOCAL MPI COMMUNICATOR
    call ed_set_MpiComm()
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    if(.not.Hmatrix%status)stop "ED_INIT_SOLVER ERROR: Hmatrix is not allocated"
    call Hij_write(unit=LOGfile)
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_lattice        !-> get state_list, independent of TEMP
    !
    do ibeta=1,Bstep
       beta = Blist(ibeta)
       if(Bstep>1)ed_file_suffix=str(beta)
       call partition_function_lattice !-> get trimmed state_list
       call observables_lattice        !-> get static observables
       call energy_lattice             !-> get energies 
       !
       if(gf_flag)call build_gf_lattice    !-> get gmatrix w/p
       if(chi_flag)call build_chi_lattice
       if(oc_flag)call build_oc_lattice
    enddo
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







