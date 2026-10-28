module deriv_det_mod
  use constants_mod, only: dp
  implicit none
  private

  public :: orb_pmap, orb_active_basis, orb_term_params, orb_psign
  public :: orb_wup, orb_wdn, orb_oup, orb_odn, orb_osum, orb_nactive
  public :: det_orb_scratch_t, det_orb_scratch

  ! Storage shared by sequential derivative evaluations in one process.
  ! common_allocations_deriv_det rebuilds the parameter map for each new setup.
  integer, allocatable, save :: orb_pmap(:,:), orb_active_basis(:), orb_term_params(:)
  real(dp), allocatable, save :: orb_psign(:,:), orb_wup(:), orb_wdn(:)
  real(dp), allocatable, save :: orb_oup(:,:), orb_odn(:,:), orb_osum(:)
  integer, save :: orb_nactive = 0

  type :: det_orb_scratch_t
    integer :: ns = 0
    logical, allocatable :: use_rank_up(:), use_rank_dn(:)
    integer, allocatable :: bases(:), mb(:), mt(:), mp(:), groups(:), where(:)
    real(dp), allocatable :: s(:,:), r(:,:), u(:,:), op(:,:), work(:,:)
    real(dp), allocatable :: ms(:), omu(:), proj(:,:), y(:,:,:)
  end type det_orb_scratch_t
  type(det_orb_scratch_t), save :: det_orb_scratch

end module deriv_det_mod
