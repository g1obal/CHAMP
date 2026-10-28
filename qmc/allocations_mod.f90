module allocations_mod

  use all_tools_mod

  contains

! ===================================================================================
  subroutine common_allocations
! -----------------------------------------------------------------------------------
! Description   : allocate arrays used in FIT, VMC and DMC
! To improved   : allocate only when necessary!
!
! Created       : J. Toulouse, 09 May 2009
! -----------------------------------------------------------------------------------
  use orbitals_mod, only: orb_tot_nb
  use all_modules_mod
  implicit none

  call object_provide ('nelec')
  call object_provide ('nelec_pair')
  call object_provide ('ncent')
  call object_provide ('nbasis')
  call object_provide ('orb_tot_nb')
  call object_provide ('nctype')
  call object_provide ('nparm')
  call object_provide ('nparmjs')

  call alloc ('phin', phin, nbasis, nelec)
  call alloc ('dphin', dphin, 3, nbasis, nelec)
  call alloc ('d2phin', d2phin, nbasis, nelec)
  
  call alloc ('d2phin_xx', d2phin_xx, nbasis, nelec)
  call alloc ('d2phin_yy', d2phin_yy, nbasis, nelec)
  call alloc ('d2phin_xy', d2phin_xy, nbasis, nelec)

  call alloc ('orb', orb, nelec, orb_tot_nb)
  call alloc ('dorb', dorb, 3, nelec, orb_tot_nb)
  call alloc ('ddorb', ddorb, nelec, orb_tot_nb)
  call alloc ('orbe', orbe, orb_tot_nb)
  call alloc ('dorbe', dorbe, 3, orb_tot_nb)
  call alloc ('ddorbe', ddorbe, orb_tot_nb)
  
  call alloc ('d2orb_xx', d2orb_xx, nelec, orb_tot_nb)
  call alloc ('d2orb_yy', d2orb_yy, nelec, orb_tot_nb)
  call alloc ('d2orb_xy', d2orb_xy, nelec, orb_tot_nb)

  call object_provide ('nupdn_square')
  call object_provide ('ndetupdn')
  call alloc ('slmui', slmui, nupdn_square, ndetupdn)
  call alloc ('slmdi', slmdi, nupdn_square, ndetupdn)
  call alloc ('fpu', fpu, 3, nupdn_square, ndetupdn)
  call alloc ('fpd', fpd, 3, nupdn_square, ndetupdn)
  call alloc ('fppu', fppu, nupdn_square, ndetupdn)
  call alloc ('fppd', fppd, nupdn_square, ndetupdn)
  call alloc ('detu', detu, ndetupdn)
  call alloc ('detd', detd, ndetupdn)
  call alloc ('ddeti_deti', ddeti_deti, 3, nelec, ndetupdn)
  call alloc ('d2edeti_deti', d2edeti_deti, nelec, ndetupdn)

  call alloc ('slmin', slmin, nupdn_square, ndetupdn)
  call alloc ('detn', detn, ndetupdn)
  call alloc ('ddeti_detin', ddeti_detin, 3, nelec, ndetupdn)
  call alloc ('d2edeti_detin', d2edeti_detin, nelec, ndetupdn)

  call object_provide ('nparmd')
  call alloc ('deti_det', deti_det, nparmd)
  call alloc ('ddeti_det', ddeti_det, 3, nelec, nparmd)
  call alloc ('d2deti_det', d2deti_det, nparmd)
  call alloc ('detij_det', detij_det, nparmd, nparmd)
! only for pseudo?
  call alloc ('deti_new', deti_new, nparmd)

  call alloc ('rshift', rshift, 3, nelec, ncent)
  call alloc ('rvec_en', rvec_en, 3, nelec, ncent)
  call alloc ('r_en', r_en, nelec, ncent)
  call alloc ('rvec_ee', rvec_ee, 3, nelec_pair)
  call alloc ('r_ee', r_ee, nelec_pair)
  call alloc ('pot_ee', pot_ee, nelec)

  call alloc ('rshift_sav', rshift_sav, 3, ncent)
  call alloc ('rvec_en_sav', rvec_en_sav, 3, ncent)
  call alloc ('r_en_sav', r_en_sav, ncent)
  call alloc ('rvec_ee_sav', rvec_ee_sav, 3, nelec)
  call alloc ('r_ee_sav', r_ee_sav, nelec)

  call alloc ('d1d2a', d1d2a, nctype)
  call alloc ('d2d2a', d2d2a, nctype)
  call alloc ('d1d2b', d1d2b, 2)
  call alloc ('d2d2b', d2d2b, 2)
  call alloc ('didk', didk, nparmjs)

  call alloc ('vd', vd, 3, nelec)
  call alloc ('vj', vj, 3, nelec)
  call alloc ('div_vo', div_vo, nelec)

  call alloc ('ekineo', ekineo, nelec)
  call alloc ('ekinen', ekinen, nelec)

! jastrow
  call alloc ('fso', fso, nelec, nelec)
  call alloc ('fijo', fijo, 3, nelec, nelec)
  call alloc ('d2ijo', d2ijo, nelec, nelec)
  call alloc ('fjo', fjo, 3, nelec)
  call alloc ('fsn', fsn, nelec, nelec)
  call alloc ('fijn', fijn, 3, nelec, nelec)
  call alloc ('d2ijn', d2ijn, nelec, nelec)
  call alloc ('fjn', fjn, 3, nelec)

  call alloc ('gvalue', gvalue, nparmjs)
  call alloc ('g', g, 3, nelec, nparmjs)
  call alloc ('d2g', d2g, nparmjs)
  call alloc ('go', go, nelec, nelec, nparmjs)

! optimization
  call alloc ('denergy', denergy, nparm)

  call alloc ('try', try, NRAD)
  call alloc ('suc', suc, NRAD)
  call alloc ('trunfb', trunfb, NRAD)
  call alloc ('rprob', rprob, NRAD)
  call alloc ('rprobup', rprobup, NRAD)
  call alloc ('rprobdn', rprobdn, NRAD)
  call alloc ('ekin', ekin, NRAD)
  call alloc ('ekin2', ekin2, NRAD)

  if(ibasis.eq.3) then
   call alloc ('cdorb', cdorb, 3, orb_tot_nb)
   call alloc ('cddorb', cddorb, orb_tot_nb)
  endif

! for pseudo
  if(nloc > 0) then
    call alloc ('vps', vps, nelec, ncent, MPS_L)
  endif

  if (nparmot > 0) call common_allocations_deriv_det
  if (ndim == 2 .and. (ibasis == 4 .or. ibasis == 5 .or. ibasis == 8)) &
    call initialize_gaussian_basis

 end subroutine common_allocations

! ===================================================================================
  subroutine common_allocations_deriv_det
! -----------------------------------------------------------------------------------
! Validate the fixed setup and prepare storage for determinant parameter derivatives.
! Call again if the parameter configuration, occupations or optimization mode changes.
! -----------------------------------------------------------------------------------
  use all_modules_mod
  implicit none
  character(len=*), parameter :: ddw_here = 'common_allocations_deriv_det'
  integer :: ddw_it, ddw_ip, ddw_ib, ddw_ic, ddw_p, ddw_slave, ddw_master
  integer :: ddw_ns, ddw_cap, ddw_mcap, ddw_spin, ddw_n, ddw_nd, ddw_d, ddw_j
  integer :: ddw_row, ddw_a, ddw_nb, ddw_nm, ddw_ng
  integer :: ddw_proj_nb, ddw_proj_nm, ddw_y_ns, ddw_y_ng
  real(dp) :: ddw_cost_rank, ddw_cost_group

  if (nparmot == 0) return
  if (nparmot < 0) call die(ddw_here,'negative total orbital-parameter count')

  ! Validate setup before indexing the parameter arrays or allocating large buffers.
  if (ndim /= 2 .or. nup < 0 .or. ndn < 0 .or. nelec /= nup+ndn) &
    call die(ddw_here,'invalid electron or spatial dimensions')
  if (notype < 1 .or. nbasis < 1 .or. norb < 1) &
    call die(ddw_here,'invalid basis, orbital or parameter-type count')
  if (ibasis /= 4 .and. ibasis /= 5 .and. ibasis /= 8) &
    call die(ddw_here,'supported ibasis values are 4,5,8')
  if (inum_orb /= 0 .or. (iperiodic /= 0 .and. iperiodic /= 1)) &
    call die(ddw_here,'analytic localized basis derivatives are required')
  if (coef_is_diag .and. norb > nbasis) &
    call die(ddw_here,'invalid diagonal coefficient dimensions')

  ! Keep allocation tests separate from SIZE/SHAPE: Fortran need not short-circuit.
  if (.not. allocated(nparmo) .or. .not. allocated(iwo)) &
    call die(ddw_here,'independent-parameter arrays are not allocated')
  if (size(nparmo) < notype .or. size(iwo,2) < notype) &
    call die(ddw_here,'independent-parameter array dimensions disagree')
  do ddw_it = 1, notype
    if (nparmo(ddw_it) < -nbasis .or. nparmo(ddw_it) > nbasis) &
      call die(ddw_here,'parameter count exceeds the number of bases')
    if (abs(nparmo(ddw_it)) > size(iwo,1)) &
      call die(ddw_here,'independent-parameter list is too short')
  enddo
  if (iconstrain_gauss_orbs /= 0 .and. iconstrain_gauss_orbs /= 1) &
    call die(ddw_here,'invalid global constraint flag')
  if (any(nparmo(1:notype) < 0)) then
    if (iconstrain_gauss_orbs /= 1) &
      call die(ddw_here,'constrained parameter type requires the global constraint flag')
    if (.not. allocated(norb_constraints)) &
      call die(ddw_here,'constraint-count array is not allocated')
    if (size(norb_constraints) < notype) &
      call die(ddw_here,'constraint-count array is too short')
  endif

  ! These occupation tables also determine the diagonal workspace capacities.
  if (ndetup < 0 .or. ndetdn < 0) &
    call die(ddw_here,'negative spin-determinant count')
  if (nup > 0) then
    if (ndetup == 0 .or. nup > norb) &
      call die(ddw_here,'invalid spin-up determinant dimensions')
    if (.not. allocated(iworbdup)) &
      call die(ddw_here,'spin-up occupations are not allocated')
    if (size(iworbdup,1) < nup .or. size(iworbdup,2) < ndetup) &
      call die(ddw_here,'spin-up occupation dimensions disagree')
    if (any(iworbdup(1:nup,1:ndetup) < 1) .or. any(iworbdup(1:nup,1:ndetup) > norb)) &
      call die(ddw_here,'spin-up orbital index out of range')
  endif
  if (ndn > 0) then
    if (ndetdn == 0 .or. ndn > norb) &
      call die(ddw_here,'invalid spin-down determinant dimensions')
    if (.not. allocated(iworbddn)) &
      call die(ddw_here,'spin-down occupations are not allocated')
    if (size(iworbddn,1) < ndn .or. size(iworbddn,2) < ndetdn) &
      call die(ddw_here,'spin-down occupation dimensions disagree')
    if (any(iworbddn(1:ndn,1:ndetdn) < 1) .or. any(iworbddn(1:ndn,1:ndetdn) > norb)) &
      call die(ddw_here,'spin-down orbital index out of range')
  endif

  call alloc ('orb_pmap', orb_pmap, notype, nbasis)
  call alloc ('orb_psign', orb_psign, notype, nbasis)
  orb_pmap = 0
  orb_psign = 0.0_dp
  ddw_p = 0
  do ddw_it = 1, notype
    do ddw_ip = 1, abs(nparmo(ddw_it))
      ddw_ib = iwo(ddw_ip,ddw_it)
      if (ddw_ib < 1 .or. ddw_ib > nbasis) &
        call die(ddw_here,'independent basis index out of range')
      if (orb_pmap(ddw_it,ddw_ib) /= 0) &
        call die(ddw_here,'duplicate independent parameter')
      ddw_p = ddw_p + 1
      if (ddw_p > nparmot) call die(ddw_here,'inconsistent total parameter count')
      orb_pmap(ddw_it,ddw_ib) = ddw_p
      orb_psign(ddw_it,ddw_ib) = 1.0_dp
    enddo
  enddo
  if (ddw_p /= nparmot) call die(ddw_here,'inconsistent total parameter count')

  if (iconstrain_gauss_orbs == 1) then
    do ddw_it = 1, notype
      if (nparmo(ddw_it) >= 0) cycle
      if (norb_constraints(ddw_it) < 0) &
        call die(ddw_here,'negative constraint count')
      if (norb_constraints(ddw_it) == 0) cycle
      if (.not. allocated(orb_constraints)) &
        call die(ddw_here,'constraint table is not allocated')
      if (size(orb_constraints,1) < notype .or. size(orb_constraints,3) < 2) &
        call die(ddw_here,'constraint table dimensions disagree')
      if (norb_constraints(ddw_it) > size(orb_constraints,2)) &
        call die(ddw_here,'constraint table is too short')
      do ddw_ic = 1, norb_constraints(ddw_it)
        ddw_slave = orb_constraints(ddw_it,ddw_ic,1)
        ddw_master = orb_constraints(ddw_it,ddw_ic,2)
        if (ddw_slave < 1 .or. ddw_slave > nbasis .or. ddw_master == 0 .or. &
            ddw_master < -nbasis .or. ddw_master > nbasis) &
          call die(ddw_here,'constraint basis index out of range')
        ddw_master = abs(ddw_master)
        ! A direct frozen master has map zero, so its slaves stay frozen too.
        if (any(orb_constraints(ddw_it,1:norb_constraints(ddw_it),1) == ddw_master)) &
          call die(ddw_here,'constraint chains are not supported')
        if (orb_pmap(ddw_it,ddw_slave) /= 0 .or. &
            any(orb_constraints(ddw_it,1:ddw_ic-1,1) == ddw_slave)) &
          call die(ddw_here,'duplicate or active constraint slave')
        orb_pmap(ddw_it,ddw_slave) = orb_pmap(ddw_it,ddw_master)
        orb_psign(ddw_it,ddw_slave) = sign(1.0_dp,real(orb_constraints(ddw_it,ddw_ic,2),dp))
      enddo
    enddo
  endif

  call alloc ('dparam', dparam, notype, nbasis, nelec)
  call alloc ('d2param', d2param, notype, notype, nbasis, nelec)
  call alloc ('ddparam', ddparam, 3, notype, nbasis, nelec)
  call alloc ('d2dparam', d2dparam, notype, nbasis, nelec)

  if (.not. allocated(dparam) .or. .not. allocated(d2param) .or. &
      .not. allocated(ddparam) .or. .not. allocated(d2dparam)) &
    call die(ddw_here,'raw basis derivatives are not allocated')
  if (any(shape(dparam) /= [notype,nbasis,nelec]) .or. &
      any(shape(d2param) /= [notype,notype,nbasis,nelec]) .or. &
      any(shape(ddparam) /= [3,notype,nbasis,nelec]) .or. &
      any(shape(d2dparam) /= [notype,nbasis,nelec])) &
    call die(ddw_here,'raw basis derivative dimensions disagree')

  orb_nactive = count(any(orb_pmap > 0,dim=1))
  call alloc ('orb_active_basis', orb_active_basis, orb_nactive)
  ddw_p = 0
  do ddw_ib = 1, nbasis
    if (.not. any(orb_pmap(:,ddw_ib) > 0)) cycle
    ddw_p = ddw_p + 1
    orb_active_basis(ddw_p) = ddw_ib
  enddo
  call alloc ('orb_wup', orb_wup, ndetup)
  call alloc ('orb_wdn', orb_wdn, ndetdn)
  call alloc ('orb_oup', orb_oup, nparmot, ndetup)
  call alloc ('orb_odn', orb_odn, nparmot, ndetdn)
  call alloc ('orb_osum', orb_osum, nparmot)
  call alloc ('orb_term_params', orb_term_params, nparmot)

  ! One scratch pool serves the two spin sectors sequentially. Every use
  ! in spin_part is restricted to the current n, nb, nm and ng extents.
  ddw_ns = max(nup,ndn)
  ddw_cap = orb_nactive
  if (coef_is_diag) ddw_cap = min(ddw_ns,orb_nactive)
  ddw_mcap = notype*ddw_cap
  det_orb_scratch%ns = ddw_ns
  call alloc ('orb_bases', det_orb_scratch%bases, ddw_cap)
  call alloc ('orb_mb', det_orb_scratch%mb, ddw_mcap)
  call alloc ('orb_mt', det_orb_scratch%mt, ddw_mcap)
  call alloc ('orb_mp', det_orb_scratch%mp, ddw_mcap)
  call alloc ('orb_groups', det_orb_scratch%groups, ddw_mcap)
  call alloc ('orb_where', det_orb_scratch%where, nparmot)
  call alloc ('orb_s', det_orb_scratch%s, ddw_ns, ddw_ns)
  call alloc ('orb_r', det_orb_scratch%r, ddw_cap, ddw_ns)
  call alloc ('orb_u', det_orb_scratch%u, ddw_ns, ddw_mcap)
  call alloc ('orb_op', det_orb_scratch%op, ddw_ns, ddw_ns)
  call alloc ('orb_work', det_orb_scratch%work, ddw_ns, ddw_mcap)
  call alloc ('orb_ms', det_orb_scratch%ms, ddw_mcap)
  call alloc ('orb_omu', det_orb_scratch%omu, ddw_mcap)

  call alloc ('orb_use_rank_up', det_orb_scratch%use_rank_up, ndetup)
  call alloc ('orb_use_rank_dn', det_orb_scratch%use_rank_dn, ndetdn)
  det_orb_scratch%use_rank_up = .false.
  det_orb_scratch%use_rank_dn = .false.
  ddw_proj_nb = 0
  ddw_proj_nm = 0
  ddw_y_ns = 0
  ddw_y_ng = 0

  ! Both contractions are exact. Select once from the parameter structure.
  ! Dense coefficients conservatively include every active physical basis,
  ! so changing coefficient values or their zero pattern cannot exceed these
  ! capacities. No unused worst-case alternative buffer is allocated.
  if (mod(iopt,10) == 2) then
    do ddw_spin = 1, 2
      if (ddw_spin == 1) then
        ddw_n = nup
        ddw_nd = ndetup
      else
        ddw_n = ndn
        ddw_nd = ndetdn
      endif
      if (ddw_n == 0) cycle

      do ddw_d = 1, ddw_nd
        ddw_nb = 0
        if (coef_is_diag) then
          do ddw_j = 1, ddw_n
            if (ddw_spin == 1) then
              ddw_ib = iworbdup(ddw_j,ddw_d)
            else
              ddw_ib = iworbddn(ddw_j,ddw_d)
            endif
            if (.not. any(orb_pmap(:,ddw_ib) > 0)) cycle
            ddw_nb = ddw_nb + 1
            det_orb_scratch%bases(ddw_nb) = ddw_ib
          enddo
        else
          ddw_nb = orb_nactive
          det_orb_scratch%bases(1:ddw_nb) = orb_active_basis
        endif
        if (ddw_nb == 0) cycle

        ddw_nm = 0
        ddw_ng = 0
        det_orb_scratch%where = 0
        do ddw_row = 1, ddw_nb
          ddw_ib = det_orb_scratch%bases(ddw_row)
          do ddw_it = 1, notype
            ddw_a = orb_pmap(ddw_it,ddw_ib)
            if (ddw_a == 0) cycle
            ddw_nm = ddw_nm + 1
            if (det_orb_scratch%where(ddw_a) /= 0) cycle
            ddw_ng = ddw_ng + 1
            det_orb_scratch%where(ddw_a) = ddw_ng
          enddo
        enddo

        ddw_cost_rank = real(ddw_nb,dp)*real(ddw_nm,dp)*real(ddw_n,dp) + real(ddw_nm,dp)**2
        ddw_cost_group = (real(ddw_nm,dp)+real(ddw_ng,dp)**2)*real(ddw_n,dp)**2
        if (ddw_cost_rank <= ddw_cost_group) then
          if (ddw_spin == 1) then
            det_orb_scratch%use_rank_up(ddw_d) = .true.
          else
            det_orb_scratch%use_rank_dn(ddw_d) = .true.
          endif
          ddw_proj_nb = max(ddw_proj_nb,ddw_nb)
          ddw_proj_nm = max(ddw_proj_nm,ddw_nm)
        else
          ddw_y_ns = max(ddw_y_ns,ddw_n)
          ddw_y_ng = max(ddw_y_ng,ddw_ng)
        endif
      enddo
    enddo
  endif

  call alloc ('orb_proj', det_orb_scratch%proj, ddw_proj_nb, ddw_proj_nm)
  call alloc ('orb_y', det_orb_scratch%y, ddw_y_ns, ddw_y_ns, ddw_y_ng)
  end subroutine common_allocations_deriv_det

end module allocations_mod
