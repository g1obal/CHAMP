      subroutine read_input
! Written by Cyrus Umrigar
      use all_tools_mod
      use constants_mod
      use variables_mod
      use basis_mod, only : which_analytical_basis, l_purely_analytical_basis
      use control_mod
      use montecarlo_mod
      use orbitals_mod
      use orbpar_mod
      use optimization_mod
      use fitdet_mod
      use orbital_grid_mod
      use atom_mod
      use dorb_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use basis1_mod
      use contrl_mod
      use const_mod
      use const2_mod
      use dim_mod
      use numbas_mod
      use basis2_mod
      use contr2_mod
      use gradhess_mod
      use contrl_opt2_mod
      use forcepar_mod
      use wfsec_mod
      use doefp_mod
      use pseudo_mod
      use contrl_per_mod
      use contrl_opt_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      use bparm_mod
      use jasparread_mod
      use pointer_mod
      use contr3_mod
      use periodic_mod
      use qua_mod
      use optimo_mod
      use jel_sph2_mod
      use contr_names_mod
      use contr_ylm_mod
      use pars_mod
      use jaspar1_mod
      use jaspar2_mod
      use ncusp_mod
      use contrldmc_mod
      use confg_mod
      use rlobxy_mod
      use pairden_mod
      use fourier_mod
      use zigzag_mod
      use branch_mod
      use periodic2_mod
      use periodic_1d_mod
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-4)

      character*80 fmt
      character*30 section
      character*10 eunit
      character*80000 input_line
      
      logical insidelist !GO

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /cyldot/ cyldot_v, cyldot_s, cyldot_rho !GO
      common /gndot/ gndot_v0, gndot_rho, gndot_s, gndot_k !GO
      common /dotcenter/ dot_bump_height, dot_bump_radius, dot_bump_radius_inv2
      common /wire/ wire_w,wire_length,wire_length2,wire_radius2, wire_potential_cutoff,wire_prefactor,wire_root1
      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      common /angularpert/ ang_perturb,amp_perturb,shrp_perturb,omg_perturb,iperturb
      common /compferm/ emagv,nv,idot

! These commons for reading fit input.  We should separate these into another
! subroutine that is called both from fit and read_input.
!     common /fit/ nsig,ncalls,iopt,ipr_opt

!     namelist /opt_list/ igradhess
      namelist /opt_list/ iring_coulomb, iantiferromagnetic, iper_gaussian_type, xmax,xfix,fmax1,fmax2,rring,ifixe,nv,idot,ifourier &
     &,iperturb,ang_perturb,amp_perturb,shrp_perturb,omg_perturb,rmin,rmax,nmeshr,nmesht,icoosys,dot_bump_height,dot_bump_radius &
     &,nmeshk1,izigzag,zzdelyr,gndot_k,gauss_width_max

      common /jel_sph1/ dn_background,rs_jel,radius_b ! RM

      dimension cent_tmp(3)
      integer, allocatable :: iflag(:)

      character*25 lhere

! Inputs not described in mainvmc:
! The first line of input is fixed-format, all the rest are free.
! title      title
! irand_seed random number seeds (four 4-digit integers)
! ijas       form of Jastrow. (between 1 and 6, mostly we use 4)
! isc        form of scaling function for ri,rj,rij in Jastrow (between 1 and 10, mostly use 2,4,6,7,16,17)
! iperiodic  0  finite system,
!            >0: numer of dimensions in which system is periodic
!            =1 system periodic in one dimension
!            =3 system periodic in three dimensions

! ibasis     =1 localized Slater or gaussian or gauss-Slater or numerical basis
!            =2 planewave basis, also for extended orbitals on grid
!            =3 complex basis for 2D quantum dots / composite fermions (if numr=0 then Fock-Darwin basis)
!            =4 2d localized floating gaussians in cartesian coord. (wigner crystal).
!            =5 2d localized floating gaussians in circular coord. (wigner crystal)
!            =6 2d localized non-circular floating gaussians in cartesian coord. (wigner crystal)
!            =7 2d localized periodic floating gaussians in cartesian coord (periodic 1d wigner crystal)
!            Warning I would like to be able to use for dots a gaussian radial basis with complex
!            spherical harmonics, but at present we cannot do that because the radial and angular bases are tied together.
! which_analytical_basis  If ibasis==1, then this is used to select between slater, gaussian and gauss-slater
! hb         hbar=0.5 for Hartree units
! etrial     guess for energy
! eunit      'Hartree'
! nstep      number of steps per block
! nblk       number of blocks
! nblkeq     number of equilibration blocks
! nconf      target number of MC configurations per processor in all dmc modes
! nconf_global = nconf in dmc, dmc_mov1 and dmc_mov1_mpi1 modes
! nconf_new  number of new MC configs. saved per processor.
! idump      = 1 dump restart file
! irstar     = 1 restart from restart file
! isite      if le 0, read starting MC config. in vmc from mc_configs_start
! isite      if ge 1, call sites to generate starting MC config. in vmc
! ipr        print level
! ipr_opt    print level in fit
! imetro     form of Metropolis (6 is most efficient choice for most systems)
!            1 simple algorithm with force-bias
!            6 accelerated Metropolis algorithm from Cyrus' 1993 PRL
! delta      step-size for simple algorithm
! deltar     radial step-size for accelerated algorithm
! deltat     angular step-size for accelerated algorithm
! fbias      force-bias.  (Use 1 always).
! idmc       form of dmc algorithm
!            1  simple dmc
!            2  improved dmc from Umrigar, Nightingale, Runge 1993 JCP
!            < 0, same as |idmc| but turn off branching to do vmc
! ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v
! nfprod     number of products to undo for estimating population control bias in dmc
! tau        time-step in dmc
! nloc       external potential (a positive value => nonlocal pseudopotential)
!            -9 numerical dot potential read in from potential_num (not yet implemented)
!            -7 gaussian dot potential, gndot_v0*dexp(-(r**2 / gndot_rho**2)**gndot_s)  !GO 
!            -6 cylindrical dot potential, 0.5d0*cyldot_v*( tanh(cyldot_s*(r+cyldot_rho)) - tanh(cyldot_s*(r-cyldot_rho)) ) !GO
!            -5 quadratic dot potential 0.5*w0^2*r^2 with barrier at center
!                (dot_bump_height)*exp(1 - 1/(1-(x/dot_bump_radius)^2))
!            -4 quantum wire, if finite:   Vwire(x) + 0.5 w0 * y^2
!                   if iperiodic.eq.1, then no V(x) - periodic in x direction
!            -3 Jellium sphere with nucleus at center, Ryo Maezono(RM) and Masayoshi Shimomoto(MS)
!            -2 quartic dot potential p1*x^4 + p2*y^4-2*p3*(xy)^2 + p4*(x-y)*x*y*r
!            -1 quadratic dot potential .5*w0*(r-rring)^2
!               If rring=0 then it is a dot, if it is not 0 it is a ring.
!            0  local, -Z/r
!            1  in Fahy format
!            2  in Troullier-Martins format (unformatted)
!            3  in Troullier-Martins format (formatted)
!            4  in CHAMP format (formatted)
!            5  in fhi format (formatted)
!            6  chemistry pseudopotentials in GAMESS-like format with 1 extra line (formatted)
! numr        Select either analytic radial basis functions (Slater, asymptotic, gaussian, gauss-slater) or numerical ones.
!             Can use mixed analytical and numerical radial functions?
!   3-dim systems:
!              =1 numerical radial functions read in from file basis.<ictype>.
!                 Use purely numerical radial basis functions (read in by read_bas_num).
!                 The nature of the basis functions is specified by just the l value, not nl,
!                 i.e., s,p,d... rather than 1s,2s,2p...
!             <=0 Read in 1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d,5f,5g,sa,pa,da
!                 Use purely analytical radial basis functions or mixed analytical and numerical
!              =0 The LCAO coefs. are read in the order 1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d,5f,5g,sa,pa,da
!                 For efficiency, within the code the coefs get multiplied by the normalization of the radial and the angular part of the basis fn.
!            <=-1 The order in which n1s etc. is read in is similar for numr=0 and -1,-2,-3
!                 However, the LCAO coefs are read in the order: all s's, all px's etc.
!                 For efficiency, Within the code the coefs get multiplied by the normalization of just the radial part of the basis fn., i.e., similar to numerical radial functions
!             =-1 Basis fns. upto 4p are read in, i.e.,  1s,2s,2p,3s,3p,3d,4s,4p,sa,pa,da
!                 Used for input from GAMESS
!             =-2 Basis fns. upto 4f are read in, i.e.,  1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,sa,pa,da
!             =-3 Basis fns. upto 4g are read in, i.e.,  1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d,5f,5g,sa,pa,da
!               The analytic radial basis functions can be slater, gaussian or gauss-slater.
!               In addition the ones at the end (sa,pa,da) are asymptotic functions that I have not used in a long time.
!               Whether one is using Slater or gaussian basis fns. used to be inputted by having n1s,n2s etc. be either > 0 or < 0
!               but now use which_analytical_basis to select between slater, gaussian and gauss-slater.
!               The order I read in the p functions is x,y,z, which is m=1,-1,0 and
!               the order that I read in the d functions is
!               3z^-r^2, x^2-y^2, xy, xz, yz, which is m=0,2,-2,1,-1, in order to be able to use old inputs.
!               All others (f, g, etc) are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
!               i.e. the order in which we were reading in the p functions.
!   2-dim systems:
!             = 0 The basis functions read in are: n1s,n2p(1)2p(-1),n3d(2),n3d(-2),etc since |m|=l
!             =-1 Read in (m_bas(ib),ib=1,nbasis) instead of n1s etc.
! nforce     number of geometries (i.e. # of forces +1)
! nefp       effective fluctuation potential for optimizing wf.
! w0         spring constant for quantum dot
! bext       external magnetic field in a.u. (only for quantum dots)
!            1 a.u. = (meff/epsilon_rel)^2 epsilon_0 /(2 mu_bohr) = 6.86219 Tesla for GaAs
! we         effective spring constant = sqrt(w0*w0+0.25*bext*bext)
! wire_w       "omega" in transverse harmonic potential 0.5 r^2 wire_w^2 for finite quantum wire
! wire_length  Length of finite quantum wire
! wire_radius2 "radius" of 3D cylindrical leads used to calculate confining potential at edges of finite quantum wire ( = 1/ (2 wire_w)
! wire_potential_cutoff  Number of wire-lengths of the 3D cylindrical leads on each side to consider when calculating quantum wire confining potential
! wire_prefactor  constant used in calculation of confining potential of quantum wire, saved to speed up calculation
! wire_root1  constant used in calculation of confining potential of quantum wire, saved to speed up calculation.
! emaglz     "magnetic energy" due to B-Lz coupling=-0.5*B*Lz (favor positive Lz)
! emagsz     "magnetic energy" due to B-Sz coupling=-0.5*B*glande*Sz (favor spin up)
! glande     effective lande factor (only for quantum dots) considered positive
! nquad      number of angular quadrature points for nonlocal psp.
! cyldot_v   height of the cylindrical quantum dot !GO
! cyldot_s   stiffness parameter of the cylindrical quantum dot !GO
! cyldot_rho radius of the cylindrical quantum dot !GO
! gndot_v0   height of the gaussian quantum dot !GO
! gndot_rho  radius of the gaussian quantum dot !GO
! gndot_s    stiffness parameter of the gaussian quantum dot !GO
! gndot_k    strength of optional quadratic gate potential centered on the origin, default is 0, k*r^2
! nelec      number of electrons
! nup        number of up-spin electrons
! npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big
! cutg       max g-vector in primitive cell.  Value determined by 2 factors:
!            a) must be large enough that all plane-wave components in wf. are covered
!            b) controls quality of Ewald fit for -Z/r and pseudopot. in primitive cell
! cutg_sim   max g-vector in simulation cell.
!            Controls quality of Ewald fit for 1/r in simulation cell
! alattice   lattice constant to multiply rlatt
!            (i.e., length of system if iperiodic=1)
! rlatt      lattice vectors of primitive cell
! rlatt_sim  lattice vectors of simulation cell
! rkvec_shift_latt k-shift for generating k-vector lattice, in reciprocal simulation cell units
! rkvec_shift k-shift for generating k-vector lattice, in cartesian units (computed, not input)
! nctype     number of atom/center types
! ncent      number of atoms/centers
! iwctype    specify atom-type for each atom
! znuc       nuclear charge
! cent       atom positions
! ndet       number of determinants in wavefunction
! nbasis     number of basis functions
! norb       number of orbitals
! cdet       coefficients of determinants (obsolete)
! ncsf       number of configuration state functions (CSFs)
! csf_coef   coefficients of each CSF
! ndet_in_csf  number of determinants in each CSF
! iwdet_in_csf which determinants enter in each CSF
! cdet_in_csf  coef. of determinants in each CSF
! iworbd     which orbitals enter in which determinants
! inum_orb   numerical orbitals on grid or not
!         =0 analytic orbs
!        !=0 use Lagrange-interpolated orbitals for periodic systems (4-pt in each direction)
!            and cubic spline-interpolated orbitals for finite systems
!         =4 numerical orbitals using 4-pt interpolation in each direction
!            if file orbitals_num exists, read from it, else write to it
!        =-4 numerical orbitals using 4-pt interpolation in each direction
!            if file orbitals_num exists, read from it, but do not write to it
!         =+-5 numerical orbitals using 4-pt interpolating pp-spline in each direction
!WP       =+-6 numerical orbitals using 4-pt approximating B-spline in each direction
!         =+-8 numerical orbitals using 4-pt interpolating B-spline in each direction
! iorb_used =1 compute only occupied orbitals
!           =0 compute all occupied and virtual orbitals if necessary (for orbital optimization)
! iorb_format 'tm' for orbitals file orbitals_pw_tm generated from Jose-Luis Martins' program
!             'pwscf' for orbitals file orbitals_pw_tm generated from PWSCF.
! ianalyt_lap analytic laplacian or not
! ijas       form of Jastrow. (between 1 and 6, mostly we use 4)
! isc        form of scaling function for ri,rj,rij in Jastrow (between 1 and 7, mostly use 2,4,6,7)
!          2  [1-exp(scalek*r)]/scalek
!          3  [1-exp{-scalek*r-(scalek*r)**2/2}]/scalek
!          4  r/(1+scalek*r)
!          5  r/{1+(scalek*r)**2}**.5
!          6  Short-range version of 2 (range given by cutjas)
!          7  Short-range version of 4 (range given by cutjas)
!          8  [1-exp(scalek*r)]
!          10 scalek*r/(1+scalek*r)

! ifock    0  no Fock's terms
!          1  phi20-like + phi21
!          2  phi20-like + phi21 + phi31-like terms
!          3  phi20-like + phi21 + phi31 + scale3
!          4  phi20-like + phi21 + phi31 + scale3 + phi20 + scale20

! nspin2   1,2,3,-1,-2 -> nspin2b=abs(nspin2)
! nspin2   > 0  nspin2 sets of a, c parms, nspin2b sets of b parms
!             nocuspb=0  parallel e-e cusp conditions satisfied (b=1/2,1/4)
! nspin2   < 0  -> nspin2=1
!               nspin2=1 sets of a and c parms, nspin2b sets of b parms
!               -1 nocuspb=1 parallel e-e cusp conditions not satisfied (1/2,1/2)
!               -2 nocuspb=0 parallel e-e cusp conditions satisfied (1/2,1/4)
! nord     order of the polynmial
! norda    order of the e-n polynmial in Jastrow4
! nordb    order of the e-e polynmial in Jastrow4
! nordc    order of the e-e-n polynmial in Jastrow4
! cjas1    simple jastrow1 (0.5 to satisfy cusps, parallel-spins automatically take half this value)
! cjas2    simple jastrow1 parameter
! scalek   scale factor for Jastrow
! a1,a2    Jastrow parameters for Jastrow2
! a,b,c    Jastrow parameters for Jastrow3
! a4,b,c   Jastrow parameters for Jastrow4,5,6
! cutjas   cutoff for Jastrow4,5,6 if isc=6,7
! rlobx(y) Lobachevsky parameters for Fock expansion
! iantiferromagnetic  Whether to constrain optimization so that
!           floating gaussians have antiferromagnetic order
!           Not implemented for odd number of electrons, or if nup \= ndn
!           = 1 by default for nup=ndn, = 0 otherwise
! iper_gaussian_type  Type of floating gaussian to use for periodic wires
!           = 1 for standard gaussian, periodic by summing images in each cell
!           = 2 (default) for form similar to that used in rings, made periodic with cos term
! iring_coulomb  whether to use a modified form of the coulomb potential
!                in rings
!           = 0  (default), V ~ 1/|rvec_i - rvec_j|
!           = 1: V ~ 1/r, where "r" is given by r^2 = "x"^2 + "y"^2,
!                   "x" = r_1 - r_2
!                   "y" = 2 * rring * sin((theta_1 - theta_2) / 2)
!            (useful for making "ringlike" geometry less important when
!               studying qpc's in rings.)
! ifixe   if > 0: which electron is fixed at a given position. 0 means none.
!          -1: calculate 2d density (not pair density)
!          -2: calculate full 2d pair density (no fixed electron)
!          -3: calculate full 2d pair density  AND 2d density.
! xfix(3)  can represent 2 different things:
!          if ifixe>0   :  coordinates of fixed electron
!                  -2/-3:  full pair-densities is written. 
!                          xfix(1) and xfix(2) are x and y coordinates of the fixed mesh point. !GO
!                          xfix(3) is lookup intervals in degrees. Rotation angle between fixed points. !GO
!                          For example; if xfix(3)=360 (no symmetry), only one mesh point is fixed.
!                                       if xfix(3)=60 (hexagon symmetry) 6 mesh points are fixed, the angle between the points is 60 degrees.   
!                                       if xfix(3)<0.6 (circular symmetry) 600 mesh points are fixed.
! icoosys  coordinate system for calculating 2d density/pairdensity
!          1: cartesian
!          2: circular
! rmin,rmax radius limits for circular coordinates
! nmeshr   2*nmeshr+1 is the number of radial mesh points for circular coordinates around rring
! nmesht   2*nmesht+1 is the number of angular mesh points for circular coordinates
! ifourier 1: "internal(?) fourier transform" of the 2 dimensional density is performed
!              (ie, power spectrum of density f(r,k_theta) for rings or f(y, k_x) for periodic wires)
!              The units of f are such that a density oscillation with period L/n
!              (where L is the length of the system; L=2*pi for rings and L=alattice for wires)
!              will cause the power spectrum have a peak at k = n.
!          0: ... is not performed (default value)
! nmeshk1   nmeshk1+1 is the number of mesh points in k space for the power spectrum
!           Note that it should be set so that delk1 = fmax1/nmeshk1 is an integer >=1, since the variable we are
!            fourier transforming is periodic over a length L (ie, n = 1.)
! izigzag  1: Write out zigzag amplitude (-1^i y_i) - quantity like staggered magnetization
!          2: Also write out quantities related to zigzag phase transition, namely reduced
!               pair density pden(r_i - r_j, theta_i - theta_j) if rings or
!               pden(y_i - y_j, x_i - x_j) if wires, as well as den(y_i - y_j, i-j)
!          0: do not write out these quantities (default value)
! zzdelyr  sets the size of delta r (or delta y) when plotting zigzag
!            pair density if izigzag = 2
! idot     0: pure complex quantum dots
!          1: dots with composite fermions
!          2: dots with laughlin wave functions
!          3: dots with projected composite fermions
! nv       2nv is the vorticity (or number of vortices per electron)
!          of composite fermions.nv is usually denoted as p in the litterature.
!          also used for laughlin wave functions (m = 2 nv + 1).
! emagv    vortices angular momentum magnetic energy due to the extra-momentum
!          carried by vortices of composite fermions (or laughlin wfs).
! rring    >0 quantum ring with potential given by 0.5 w0^2 (r-rring)^2
! iperturb 0: off (no angular perturbation for quantum rings)
!          1: smoothed-square-like perturbation for quantum rings
!     = amp_perturb*(tanh(shrp_perturb*(theta+ang_perturb))-tanh(shrp_perturb*(theta-ang_perturb)))
!             and for quantum wires: (ang_perturb becomes semi-length of constriction)
!     = amp_perturb*(tanh(shrp_perturb*(x(1,i)+ang_perturb))-tanh(shrp_perturb*(x(1,i)-ang_perturb)))
!          2: gaussian-like perturbation for quantum rings
!     = 2.d0*amp_perturb*dexp(-0.5d0*(theta/(ang_perturb))**2)
!              and for quantum wires:  (Note this is NOT periodic yet if iperiodic.eq.1)
!     = 2.d0*amp_perturb*dexp(-0.5d0*(x(1,i)/(ang_perturb))**2)
!          3: upside-down parabola, to give simple 1/2 w_y^2 y^2 - 1/2 w_x^2 x^2
!               saddle point for quantum wires.  Multiplied by a bump
!               function similar to that used in iperturb = 1 to make
!               a smooth connections to the leads.
!     = (2.d0*amp_perturb - 0.5d0*(omg_perturb**2)*((theta*rring)**2)
!       *0.5*(tanh(shrp_perturb*(theta+ang_perturb))-tanh(shrp_perturb*(theta-ang_perturb)))
!        where ang_perturb is set to be the value of theta where
!       (2.d0*amp_perturb - 0.5d0*(omg_perturb**2)*((theta*rring)**2) = 0
!       and shrp_perturb is hard-coded to make the bump function rise
!       much more quickly than its width.
!      and for quantum wires:
!     = (2.d0*amp_perturb - 0.5d0*(omg_perturb**2)*(x(1,i)**2)
!        *0.5*(tanh(shrp_perturb*(x(1,i)+ang_perturb))-tanh(x(1,i)*(theta-ang_perturb)))
!
! ang_perturb  angular perturbation range
! amp_perturb  amplitude of the angular perturbation
! shrp_perturb sharpness of the angular perturbation
! omg_perturb  omega for upside down parabola perturbation
!
! gauss_width_max  maximum width for 2d floating gaussian basis set (ibasis=4), optional parameter
!                  < 0: no limit (default)
!                  >= 0: limit with the value of gauss_width_max

! Optimization parameters:

! 1) 10 1000 1.d-6 0. 1.d-3     nopt_iter,nblk_max,add_diag(1),p_var,tol_energy
! This line is used only by vmc and dmc programs, not by fit
! nopt_iter   Number of wavefunction optimization steps to be done in vmc or dmc programs
!             = 0 one vmc run, no grad/hess
!             < 0 one vmc run, but calculate grad/hess, hamiltonian/overlap and save to file
!             > 0 do nopt_iter optimization steps, i.e., nopt_iter+1 vmc runs (unless convergence reached earlier)
! nblk_max    For the optimization the program starts with with nblk blocks and increases it
!             upto nblk_max if that is needed to get the desired to_energy
! add_diag(1) The starting value of add_diag.
!             If add_diag < 0, do not do correlated sampling runs to adjust add_diag automatically.
!             Instead use abs(add_diag) always.
! p_var       Specifies the linear combination of energy and variance to be optimized.
!          0. Optimize energy
!          1. Optimize variance
! tol_energy  The desired error bar on the energy.  This variable is also used to see if the
!             energy difference between the 3 wavefunctions (1 primary and 2 secondary) is small
!             enough to say that the wavefunction optimization is converged.

! 2) 4800 163 1 1 5 1000 21101 1 NDATA,NPARM,icusp,icusp2,NSIG,NCALLS,iopt,ipr_opt
! ndata       Number of MC configurations (data points) to optimize over in fit
! nparm       Number of parameters to be optimized
!             Both if we are doing optimization in fit or in vmc, this number can get changed within code.
!    icusp  = <= -1  zeroth order e-N cusp not imposed
!             >=  0  zeroth order e-N cusp imposed "exactly"
!    icusp2   is for imposing e-N and e-e cusps
!             >=  1 impose cusps exactly
!             <= -2 impose/check cusps via penalty (used mostly for ijas=2 for which exact imposition is difficult)
!                The penalty wt. cuspwt depends on whether icusp2 is divisible by 3 or not.
!                e.g. icusp2=-101 gives a wt of 100, icusp2=-102 gives a wt of .01
!                     icusp2=-1001 gives a wt of 1000, icusp2=-1002 gives a wt of .001
! nsig        Crude estimate of number if significant digits
! ncalls      Maximum number of calls to func permitted in fit
! iopt        In fit:
!               Choose between zxssq and quench
!             = 0,1 call zxssq from IMSL (obsolete)
!             = 2   call quench written by Peter Nightingale and Cyrus Umrigar
!             = 0 Do not check for strict downward descent in zxssq
!             = 1 Strict downward descent in zxssq
!               zxssq is obsolete so, if in fit mode, we reset iopt to 2 in read_input
!             In vmc:
!               Choose between linear, Newton and perturbation theory (if not using Julien's menu; one should use his)
!             = 1 linear method
!             = 2 modified Newton method
!             = 3 perturbation theory
!               At present other digits also have the foll. meaning for linear and Newton (but these will change)
!                             linear                          |           Newton
!             10 digit = 0   nonsym Ham (good choice)         |       nonsym Hess (not implemented)
!                        1      sym Ham (bad  choice)         |          sym Hess (not implemented)
!                        2      sym Ham (with covariances,    |
!                               not as bad as 1)              |
!             100      = 0   nonorthog basis                  |           rescale Hess
!                        1   semiorthog basis                 |       not rescale Hess
!             1000     = 0   use state with largest 0 coef    |
!                        1   use state with lowest eig        |
!             10000    = 0   c_i/c_0 for everything                           |
!                        1   c_i/(c_0-sum_1^{Nparm } S_i0*c_i) for everything |
!                        2   c_i/(c_0-sum_1^{NCSF-1} S_i0*c_i) for everything |
!                        3   c_i/(c_0-sum_1^{NCSF-1} S_i0*c_i) for CSF        |
!                            c_i/c_0                           for J          |
!                        4   c_i/(c_0+sum_1^{Nparmj} a_i*c_i)  for everything |
!                            where a_i=(S_0i*c_0 + sum_j^Nparmj S_ij*c_j)
!                                      ----------------------------------
!                                      (S_00*c_0 + sum_j^Nparmj S_0j*c_j)
!                            Note that the sum 4 lines up and the sume in the def. of a_i
!                            are over the nonlin parms only
!                            In practice c_0 is set to 1.
!                        5   c_i/(c_0+sum_1^{Nparmj} a_i*c_i)  for everything |
!                            where a_i=(S_0i*c_0 + sum_j^Nparm S_ij*c_j)
!                                      ----------------------------------
!                                      (S_00*c_0 + sum_j^Nparm S_0j*c_j)
!                            Note that the sum 4 lines up is over the nonlin parms only
!                            but the sums in the def. of a_i are over all parameters.
!                        In practice c_0 is set to 1.
!             So, good choices are:
!             = 21101 for linear
!             = 31101 for linear
!             = 41001 for linear (2nd semiorthog. -- too small changes)
!             = 51001 for linear (2nd semiorthog. -- too small changes)
!             =     2 for modified Newton with rescaling of J params.
!             =   102 for modified Newton without rescaling of J params (not as good when far from min).
!             =     3 perturbation theory (little used in this version; used more in Julien's version).
!             Bad choices useful for testing are
!             = 21111 for linear (symmetrized Nightingale)
!             = 21121 for linear (symmetric Hamiltonian different from symmetrized Nightingale (a bit better))
! ipr_opt     Print flag for optimization
!             <= -2  Minimal print out
!             >= -1  Print cusp monotonicity conditions
!             >=  0  Print configs and errors on 6
!             >=  2  Print out configs and wavefunction on 2

!    isc      is used for 3 purposes
!             a) for using spq rather than srt variables
!             b) for calling numerical (jastrow_num,psi) rather than analytical (jastrow)
!                gradient and Laplacian. Analytic if >= -5
!             c) for scaling ri,rj,rij in Pade.
!                -2,-7  (1-exp(scalek(1)*r)/scalek(1)
!                -3,-8  [1-exp{-scalek(1)*r-(scalek(1)*r)**2/2}]/scalek(1)
!                -4,-9  r/(1+scalek(1)*r)
!                -5,-10 r/{1+(scalek(1)*r)**2}**.5

! 3) 0 0 0 0 i3body,irewgt,iaver,istrech (rarely used)
! irewgt      = 0 do not reweight the local energy differences in fit
!             > 0 reweight by psi_currrent^2/psi_sampled^2, but limit the weights to wtmax=iabs(mod(irewgt,100))
!               This option allows the nodes of the wavefn. to go from one side of a sampled config to the other.
!               Without it, the divergence in the local energy at the node blocks it.
! 4) 0 0 0 0 0 0 0 0 0 0 ipos,idcds,idcdr,idcdt,id2cds,id2cdr,id2cdt,idbds,idbdr,idbdt (rarely used)
! 5) 0 0 0 0 0 1 1 1 1 0 1 1 0 0 1 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 0 1 1 0 0 0 0 0 1 1 0 1 1 0 0 0 (lo(iorb),iorb=1,norb) (calculated internally if necn < 0)

! 6) 102  4  5  15  0  37 0 0  nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg
! nparml      Number of LCAO coefs to be optimized
! nparma      Number of en  Jastrow coefs to be optimized (one entry for each atom type)
! nparmb      Number of ee  Jastrow coefs to be optimized
! nparmc      Number of een Jastrow coefs to be optimized (one entry for each atom type)
! Note: nparma and nparmc are specified for each center type
!       nparmb are specified for both spin types if nspin2=2
! nparmf      Number of Fock Jastrow coefs to be optimized (not yet implemented for Jastrow4)
! nparmd      Number of determinantal coefs to be optimized (obsolete)
! nparmcsf    Number of determinantal coefs to be optimized
! nparms      Number of Jastrow scale factor coefs to be optimized (0 or 1)
! nparmo(i)   Number of Orbital parameters of type i to be optimized
!             At present this is used for floating
!             gaussians and there are 3-4 types (x,y positions and width).
!             setting this to be negative indicates that there are constraints, the
!             absolute value is the number of free parameters
! nparmg      Do not use this.
! norb_constraints(i)  Number of constraints imposed on orbital params of type i to be optimized
! orb_constraints(type,constraint #, orb1:orb2) For each type and constraint,
!              the first orbital is set equal to the value (or negative of, if it has a '-' sign)
!              of the second orbital (orb2 must be one of the orbitals being optimized, ie.
!              it should be in iwo).
!    Eg, sample input if orbitals 1 and 3 have the same first coordinate, and
!         if the second coordinate of 2 had mirror symmetry with 30, and 3 with 29:
!  12 0 0    (norb_constraints(i),i=1,notype)  ! a total of 3 constraints
!  31             ((orb_constraints(1,i,j),j=1,2),i=1,norb_constrants(1))
!  30 -2   29 -3   ((orb_constraints(1,i,j),j=1,2),i=1,norb_constrants(1))
!                  ((orb_constraints(1,i,j),j=1,2),i=1,norb_constrants(1))
!                  ((orb_constraints(1,i,j),j=1,2),i=1,norb_constrants(1))

! 7) 1  1   1  3   1  4   1  5   1  6   1 15   ... (iworb(ibasis),iwbasi(ibasis),ibasis=1,nparml)
! 8) 1  2  3  4  5  6  7  8  9 10 19 (iwbase(ibasis),ibasis=1,nparme)
! 9)   2 3 4 5 6 7 8 9 10 ... (iwcsf(iparm),iparm=1,nparmcsf)

! 10)     3 4 5 6 (iwjasa(iparm),iparm=1,nparma) (one entry for each atom type)
! which e-n Jastrow parameters are varied

! 11) 2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)
! which e-e Jastrow parameters are varied

! 12)     3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 (iwjasc(iparm),iparm=1,nparmc) (one entry for each atom type)
! which e-e-n Jastrow parameters are varied (the above is for 5th-order, but code can handle any order)

! 13) 165 35    necn,nebase
! necn      Number of LCAO coefs set equal.  If it is -1, necn and nparml will be calculated within code
! nebase    Number of basis exponents .  If it is -1, nebase and nparme will be calculated within code

! 14) 1 24   1  1    1 25   1  2    ... ((ieorb(iorb,ibasis),iebasi(iorb,ibasis),iorb=1,2),ibasis=1,necn)
! 15) 11  7   12  8   13  9   14 10  ... ((iebase(iorb,ibasis),iorb=1,2),ibasis=1,nebase)
! 16) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (ipivot(j),j=1,norb) (not needed if necn < 0)
! 17) -109.51 eave (guess for energy, used in fit)
! 18) 1.d-6 5. 1 50 1 pmarquardt,tau,noutput,nstep,ibold (used in fit)
! 19) T F analytic,cholesky (used in fit)

! former ewald.h:
! NCOEFX     max number of coefficients in short-range polynomial for optimal Ewald (not used anymore)
! NPX        is the number of factors of (r/cutr-1) we multiply polynomial by (not used anymore)
! MKPTS      maximum number of k-pts.  This can be at most the ratio of the simulation to primitive cell,
!            volumes, V_sim/V.  However, since some k-pts have two indep. states (most if V_sim/V>>1)
!            the number of k-pts is between V_sim/(2V) and V_sim/V.  I have not tried to distinguish
!            between arrays that could be dimensioned to nkpts and those that need V_sim/V, but
!            have used the larger number for all.
! IVOL_RATIO is the ratio of the simulation to primitive cell volumes
!            However, since (cutg_sim/cutg) and (cutg_sim_big/cutg_big) can be chosen to
!            be (V_sim/V)^(-1/3), this ratio can just be set to 1 and I could remove this parameter.
! IBIG_RATIO is the ratio of the number of vectors or norms before the optimal Ewald separation
!            to the number after the separation. This ratio is (cutg_big/cutg)^3
!            and (cutg_sim_big/cutg_sim)^3 for the primitive and simulation cells resp.
!            This ratio is accurate for vectors and an overestimate for norms.
! NSYM       is the ratio of the number of vectors to the number of norms
!            and depends on the symmetry of the lattice.  Since many stars have less
!            vectors than the number of symmetry operations, this number can be set
!            somewhat smaller than the number of sym. ops.  e.g., for cubic one can
!            try 32 instead of 48.
!      parameter(NCOEFX=20, NPX=4, MKPTS=32, IVOL_RATIO=1, IBIG_RATIO=15, NSYM=32)
!      parameter(NGNORMX=700, NGVECX=NGNORMX*NSYM, NGVEC2X=2*NGVECX, NG1DX=60)
!      parameter(NGNORM_SIMX=NGNORMX*IVOL_RATIO, NGVEC_SIMX=NGVECX*IVOL_RATIO)
!      parameter(NGNORM_BIGX=IBIG_RATIO*NGNORMX, NGVEC_BIGX=IBIG_RATIO*NGVECX)
!      parameter(NGNORM_SIM_BIGX=IBIG_RATIO*NGNORM_SIMX, NGVEC_SIM_BIGX=IBIG_RATIO*NGVEC_SIMX)


      lhere = 'read_input' ! JT

      pi=four*datan(one)

!     correlated sampling index
      iwf = 1
      call object_modified ('iwf')

!     write(fmt,'(''(a,a'',i3,'')'')') len_trim(mode)
!     write(6,fmt) 'CHAMP version 3.08.0, mode=', mode
!     write(6,'(''Running in mode: '',a)') trim(mode)

      write(6,*)
      if(mode.eq.'fit') write(6,'(''Wavefn. optimization'')')
      if(mode.eq.'fit_mpi') write(6,'(''Wavefn. optimization mpi'')')
      if(mode.eq.'vmc') write(6,'(''Variational MC'')')
      if(mode.eq.'vmc_mov1') write(6,'(''Variational MC one-electron move'')')
      if(mode.eq.'vmc_mov1_mpi') write(6,'(''Variational MC one-electron move mpi'')')
      if(mode.eq.'dmc') write(6,'(''Diffusion MC'')')
      if(mode.eq.'dmc_mov1') write(6,'(''Diffusion MC 1-electron move'')')
      if(mode.eq.'dmc_mov1_mpi1') write(6,'(''Diffusion MC 1-electron move, mpi no global pop'')')
      if(mode.eq.'dmc_mov1_mpi2') write(6,'(''Diffusion MC 1-electron move, mpi global pop big comm'')')
      if(mode.eq.'dmc_mov1_mpi3') write(6,'(''Diffusion MC 1-electron move, mpi global pop small comm'')')

!     initializations
      nparm=0
      nparmd = 0
      call object_modified ('nparm')
      call object_modified ('nparmd')
      nparmj=0
      nparms=0
      nparmjs=nparmj+nparms
      call object_modified ('nparmjs')
      nwf = 3 ! for optimization runs
      call object_modified ('nwf')

!     read(5,'(a20,4x,4i4)') title,irand_seed
      read(5,*) title
      write(fmt,'(''(a'',i3,'')'')') len_trim(title)
      write(6,fmt) title

      read(5,'(4i4)') irand_seed
      read(5,*) iperiodic,ibasis,which_analytical_basis
      call object_modified ('which_analytical_basis')
      if(iperiodic.gt.0) then
        write (6,'(''System periodic in '',i1,'' dimensions'')') iperiodic
      else
        write(6,'(''Finite system'')')
      endif
      if(ibasis.eq.1) then
        write(6,'(''Localized Slater or gaussian or gauss-slater or numerical basis'')')
        write(6,'(2a)') 'If there are any analytical radial basis functions then they are ', trim(which_analytical_basis)
        if(trim(which_analytical_basis).ne.'slater' .and. trim(which_analytical_basis).ne.'gaussian' &
     &  .and. trim(which_analytical_basis).ne.'gauss-slater') then
          write(6,'(''which_analytical_basis must be either slater or gaussian or gauss-slater'')')
          stop 'which_analytical_basis must be either slater or gaussian or gauss-slater'
        endif
        notype=0
      elseif(ibasis.eq.2) then
        write(6,'(''Planewave basis, or extended orbitals on grid'')')
        notype=0
      elseif(ibasis.eq.3) then
        write(6,'(''Complex basis for 2D quantum dots / composite fermions'')')
        notype=0
      elseif(ibasis.eq.4) then
        write(6,'(''Floating Gaussian basis for 2D Wigner crystals'')')
        notype=3
      elseif(ibasis.eq.5) then
        write(6,'(''Floating Gaussian basis for 2D Wigner crystals in ring geom.'')')
        notype=4
      elseif(ibasis.eq.6) then
        write(6,'(''Floating, non-circular Gaussian basis for 2D Wigner crystals'')')
        notype=4
      elseif(ibasis.eq.7) then
        if(iperiodic.ne.1) then
           write(6, '(''read_input: iperiodic must =1 for ibasis=7'')')
           stop 'read_input: iperiodic must =1 for ibasis=7'
        endif
        write(6,'(''Floating, periodic Gaussian basis for periodic (quasi-)1D Wigner crystals'')')
        notype=4
      else
        write(6,'(''read_input: ibasis must be between 1 and 7'')')
        stop 'read_input: ibasis must be between 1 and 7'
      endif

!     if(index(mode,'vmc').ne.0 .and. iperiodic.gt.0) stop 'In order to do VMC calculation for periodic system run dmc or dmc.mov1 with idmc < 0'

!     if(index(mode,'vmc').ne.0 .or. index(mode,'dmc').ne.0) then
        write(6,'(/,''random number seeds'',t25,4i4)') irand_seed
        call setrn(irand_seed)
!     endif

      read(5,*) hb,etrial,eunit
      write(6,'(''hbar**2/(2.*m) ='',t31,f10.5)') hb
      write(6,'(''etrial'',t29,f12.6)') etrial
      write(6,'(''all energies are in'',t31,a10)') eunit

      dos_dele = etrial/NAX   !set spacing for single-particle DOS histogram

      if(index(mode,'vmc').ne.0 .or. index(mode,'dmc').ne.0) then
        read(5,*) nstep,nblk,nblkeq,nconf,nconf_new
        write(6,'(''no. of steps/block ='',t31,i10)') nstep
        write(6,'(''no. of blocks after eq.='',t31,i10)') nblk
        write(6,'(''no. of blocks before eq. ='',t31,i10)') nblkeq
! Warning: tmp +50 is ususally enough but increase to say +2000 for antisymmetric projector
        MWALK = nconf+50 ! default value of MWALK
!       MWALK = nconf+2000 ! default value of MWALK
        if(index(mode,'dmc').ne.0) then
          write(6,'(''target walker population/processor='',t36,i5)') nconf
          if(nconf.le.0) stop 'target population <= 0'
        endif
        write(6,'(''no. configurations saved ='',t31,i10)') nconf_new
        call object_modified ('nstep') !JT
        call object_modified ('nconf') !JT
      else
        read(5,*)
      endif
      if(mode.eq.'dmc' .or. mode.eq.'dmc_mov1' .or. mode.eq.'dmc_mov1_mpi1') then
        nconf_global=nconf
       elseif(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') then
        nconf_global=nconf*nproc
      endif

      read(5,*) idump,irstar,isite,ipr
      if(index(mode,'vmc').ne.0 .or. index(mode,'dmc').ne.0) then
        if(irstar.eq.1) nblkeq=0
        write(6,'(''no. of blocks before eq. ='',t31,i10)') nblkeq
        if(irstar.eq.1) write(6,'(''Job is starting from restart file'')')
        if(idump.eq.1) write(6,'(''Job will write restart file at end of run'')')
! Make sure that the printout is not huge
        if(nstep*(nblk+2*nblkeq).gt.104000 .and. ipr.gt.-1) then
          ipr=min(ipr,-1)
          write(6,'(''Warning: ipr set to'',i3,'' to avoid large output'')') ipr
        endif
      endif

        read(5,*) imetro,delta,deltar,deltat,fbias
        deltai=one/delta
        if(deltar.lt.one) then
          write(6,*) '**Warning value of deltar reset to 2.'
          deltar=two
        endif
        if(deltat.lt.zero .or. deltat.gt.two) then
          write(6,*) '**Warning value of deltat reset to 2.'
          deltat=two
        endif
! Truncate fbias so that it is never negative, and the quantity sampled is never negative
        fbias=dmin1(two,dmax1(zero,fbias))
        write(6,'(''Version of Metropolis ='',t31,i10)') imetro
        write(6,'(''step size ='',t31,f10.5)') delta
        write(6,'(''radial step multiplier ='',t31,f10.5)') deltar
        write(6,'(''cos(theta) step size ='',t31,f10.5)') deltat
        write(6,'(''force bias ='',t31,f10.5)') fbias
        if(imetro.ne.1 .and. imetro.ne.6) stop 'imetro must be 1 or 6 (accel. Metropolis)'
        if(imetro.ne.1 .and. iperiodic.gt.0) stop 'In order to do VMC calculation for periodic system run dmc or dmc.mov1 with &
       &idmc < 0 or run vmc with imetro=1'

! It has been updated now, but rather quickly
!     if(index(mode,'vmc_one').ne.0 .and. imetro.eq.1) stop 'metrop_mov1 has not been updated'

      if(index(mode,'dmc').ne.0) then
        read(5,*) idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
        write(6,'(/,''idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e='',9i4)') &
     &  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
        if(idmc.lt.0) write(6,'(''Running DMC program in VMC mode'')')
        if(iabs(idmc).lt.1 .or. iabs(idmc).gt.3) stop 'iabs(idmc) must be 1 or 2 or 3'
        read(5,*) nfprod,tau
        rttau=dsqrt(tau)
        write(6,'(''nfprod,tau'',t31,i5,f10.5)') nfprod,tau
      else
        read(5,*)
        read(5,*)
      endif
      call object_modified ('nfprod')

      read(5,*) nloc,numr,nforce,nefp
      write(6,'(''nloc,numr ='',t31,4i5)') nloc,numr
      write(6,'(''nforce,nefp ='',t31,4i5)') nforce,nefp
!     if(numr.gt.0) write(6,'(/,''numerical basis functions used'')')
      if(nloc.lt.-7 .or. nloc.gt.6) stop 'nloc must be between -7 and 6 inclusive' !GO
      if(nloc.ge.2) then
        read(5,*) nquad
        write(6,'(''nquad='',t31,i4)') nquad
      elseif(nloc.eq.-1) then
        read(5,*) w0,bext,glande
        we=dsqrt(w0*w0+0.25d0*bext*bext)
        write(6,'(''spring const of dot pot., w0='',t31,f10.5)') w0
        write(6,'(''applied magnetic field., bext='',t31,f10.5)') bext
        write(6,'(''effective spring const., we='',t31,f10.5)') we
        write(6,'(''Lande factor, glande='',t31,f10.5)') glande
      elseif(nloc.eq.-2) then
        read(5,*) p1,p2,p3,p4
        write(6,'(''quartic dot pot. p1,p2,p3,p4='',t31,9f9.6)') p1,p2,p3,p4
      elseif(nloc.eq.-4) then
        if(iperiodic.eq.0)then
          read(5,*) wire_w,wire_length,wire_potential_cutoff
          write(6,'(''wire_w,wire_length,wire_potential_cutoff='',t31,9f9.6)') wire_w,wire_length,wire_potential_cutoff
        else
         read(5,*) wire_w
         write(6,'(''wire_w='',t31,f10.5)') wire_w
        endif

        we=wire_w  !  this is a quick fix:  needed for the subroutine basis_fns_2dgauss

      elseif(nloc.eq.-5) then
        read(5,*) w0,bext,glande
        we=dsqrt(w0*w0+0.25d0*bext*bext)
        dot_bump_height = we*100.0d0 ! default value, this is set in optional features
        dot_bump_radius = 0.1d0 / dsqrt(we) ! default value, set in optional features
        write(6,'(''spring const of dot pot., w0='',t31,f10.5)') w0
        write(6,'(''applied magnetic field., bext='',t31,f10.5)') bext
        write(6,'(''effective spring const., we='',t31,f10.5)') we
        write(6,'(''Lande factor, glande='',t31,f10.5)') glande
       elseif(nloc.eq.-6) then
        read(5,*) cyldot_v,cyldot_s,cyldot_rho !GO
        w0 = 0.d0 ! in order to use gaussian orbitals
        bext = 0.d0 ! in order to use gaussian orbitals
        glande = 0.d0 ! in order to use gaussian orbitals
        we = 1.d0 ! in order to use gaussian orbitals
        write(6,'(''height cyldot pot., cyldot_v='',t31,f10.5)') cyldot_v
        write(6,'(''stiff. cyldot pot., cyldot_s='',t31,f10.5)') cyldot_s
        write(6,'(''rad. cyldot pot., cyldot_rho='',t31,f10.5)') cyldot_rho
       elseif(nloc.eq.-7) then 
        read(5,*) gndot_v0,gndot_rho,gndot_s !GO
        w0 = 0.d0 ! in order to use gaussian orbitals
        bext = 0.d0 ! in order to use gaussian orbitals
        glande = 0.d0 ! in order to use gaussian orbitals
        we = 1.d0 ! in order to use gaussian orbitals
        write(6,'(''height gndot pot., gndot_v0='',t31,f10.5)') gndot_v0
        write(6,'(''rad. gndot pot., gndot_rho='',t31,f10.5)') gndot_rho
        write(6,'(''stiff. gndot pot., gndot_s='',t31,f10.5)') gndot_s
      endif

      call object_modified ('nloc') ! JT
      call object_modified ('numr') ! JT
      call object_modified ('nforce') ! JT

      if(nloc.ge.2 .and. nquad.gt.MPS_QUAD) stop 'nquad > MPS_QUAD'

      read(5,*) nelec,nup
      ndn=nelec-nup
      write(6,'(/,''no. of electrons (all,up,dn) ='',t31,3i5)') nelec,nup,ndn
      if(nup.le.0) stop 'nup must be >=1'
      if(nup.lt.ndn) stop 'nup must be >=ndn'

      nupdn = max(nup, ndn)
      nup_square = nup**2
      ndn_square = ndn**2
      nupdn_square = nupdn**2
      nelec_pair = nelec*(nelec-1)/2
      call object_modified ('nup_square')
      call object_modified ('ndn_square')
      call object_modified ('nupdn_square')
      call object_modified ('nelec_pair')

      if(nloc.eq.-4 .and. iperiodic.eq.0) then ! values used for quantum wire
        wire_length2 = wire_length*wire_length
        wire_radius2 = 1 / (2 * wire_w)
        wire_prefactor = nelec / (wire_radius2 * wire_length)
        wire_root1 = dsqrt((wire_potential_cutoff **2)  * wire_length2 + wire_radius2)
      endif

      call object_modified ('nelec') !JT
      call object_modified ('nup') !JT
      call object_modified ('ndn') !JT


! The foll. does not work because mode is set in maindmc.f and the same maindmc.f is used for all-electron and 1-electron move versions.
!     if(nelec.gt.4 .and. index(mode,'one').eq.0) write(6,'(''Warning: for more'',
!    &'' than 4 electrons, you should use mov1 version of program.'',/,
!    &''Otherwise acceptance will be low and for nelec a few hundred there will be under or overflows'')')

! Geometrical section
      read(5,*) section
      write(6,'(/,a30,/)') section

      read(5,*) ndim
      write(6,'(i1,'' dimensional system'')') ndim
      if(ndim.ne.2.and.ndim.ne.3) stop 'ndim must be 2 or 3'
      if(ndim.lt.iperiodic) stop 'iperiodic must be less than or equal to ndim'
      if(iperiodic.eq.2) stop 'systems periodic in 2d not yet implemented'
      if(ndim.eq.2.and.imetro.ne.1.and.index(mode,'vmc').ne.0) &
     &stop 'imetro!=1 not yet implemented for ndim=2'

      call object_modified ('ndim') ! JT

      if(iperiodic.eq.1) then
        read(5,*) alattice
        if(nloc.eq.-4) wire_length = alattice
        write(6,'(''Length of 1d periodic system, alattice='',t31,f10.5)') alattice
      elseif(iperiodic.eq.3) then

! npoly is the polynomial order for short-range part
        read(5,*) npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big
        write(6,'(/,''Npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big'',2i4,9f8.2)') npoly,np,cutg,cutg_sim,cutg_big,cutg_sim_big
        if(npoly.ne.7) then
          write(6,'(''present version works best with npoly=7'')')
          stop 'present version works best with npoly=7'
        endif

        ncoef=npoly+1

        read(5,*) alattice
        write(6,'(''Length of 3d periodic system, alattice='',t31,f10.5)') alattice
        call alloc ('rlatt', rlatt, 3, 3)
        do 10 i=1,ndim
          read(5,*) (rlatt(k,i),k=1,ndim)
          do 10 k=1,ndim
   10       rlatt(k,i)=rlatt(k,i)*alattice

        write(6,'(/,''Lattice basis vectors'',3(/,3f10.6))') ((rlatt(k,j),k=1,ndim),j=1,ndim)

! read the dimensions of the simulation 'cube'
        call alloc ('rlatt_sim', rlatt_sim, 3, 3)
        do 20 i=1,ndim
          read(5,*) (rlatt_sim(k,i),k=1,ndim)
          do 20 k=1,ndim
   20       rlatt_sim(k,i)=rlatt_sim(k,i)*alattice

        write(6,'(/,''Simulation lattice basis vectors'',3(/,3f10.6))') ((rlatt_sim(k,j),k=1,ndim),j=1,ndim)

! read k-shift for generating k-vector lattice
        call alloc ('rkvec_shift_latt', rkvec_shift_latt, 3)
        read(5,*) (rkvec_shift_latt(k),k=1,ndim)
        do 22 k=1,ndim
   22     if(rkvec_shift_latt(k).ne.0.d0 .and. rkvec_shift_latt(k).ne..5d0) &
     &    stop 'rkvec_shift_latt components must be 0 or 1/2 to have real orbs'

      endif

      read(5,*) nctype,ncent
      write(6,'(/,''nctype,ncent ='',t31,i3,i5)') nctype,ncent
      call object_modified ('nctype')
      call object_modified ('ncent')
      
      if (nloc .eq. -6 .or. nloc .eq. -7) then !GO
        if (nelec .gt. 2 * ncent) then
          stop 'nelec cannot be greater than 2 * ncent for nloc = -6 and -7'
        endif
      endif

      call alloc ('iwctype', iwctype, ncent)
      read(5,*) (iwctype(i),i=1,ncent)
      write(6,'(''iwctype ='',t31,20i3,(20i3))') (iwctype(i),i=1,ncent)
      do 25 ic=1,ncent
   25   if(iwctype(ic).gt.nctype) stop 'iwctype(ic) > nctype'

      call object_modified ('iwctype')  !JT

! Store the number of centers of each center type in ncent_ctype()
      call alloc ('ncent_ctype', ncent_ctype, nctype)
      do 26 ict=1,nctype
   26   ncent_ctype(ict)=0
      do 27 ic=1,ncent
   27   ncent_ctype(iwctype(ic))=ncent_ctype(iwctype(ic))+1
      write(6,'(''Number of centers of each centertype: '',20(i4,1x))') (ncent_ctype(ict),ict=1,nctype) !GO

      call alloc ('znuc', znuc, nctype)

      read(5,*) (znuc(i),i=1,nctype)
      write(6,'(''znuc='',t31,10f5.1,(10f5.1))') (znuc(i),i=1,nctype)
      call object_modified ('znuc') !JT

      if(nloc.eq.-3) then ! Jellium RM
!MS Jellium sphere plus charge at center
        dn_background = nelec - znuc(1)
        rs_jel = 1.d0
        radius_b = (dn_background*(rs_jel)**3)**(1.d0/3.d0)
        zconst = 20 !* 27Aug06
      else
        zconst = 0  ! normal case
      endif

! Read in which is the local component of the potential
      if(nloc.gt.0) then
        call alloc ('lpotp1', lpotp1, nctype)
        read(5,*) (lpotp1(i),i=1,nctype)
        write(6,'(''lpotp1='',t31,20i3,(20i3))') (lpotp1(i),i=1,nctype)
        do 35 i=1,nctype
   35     if(lpotp1(i).gt.MPS_L) stop 'lpotp1(i) > MPS_L'
      endif

!     if(iperiodic.eq.0 .or iperiodic.eq.1) then
!       write(6,'(/,''center positions'')')
!      else
!       write(6,'(/,''center positions in primitive lattice vector units'')')
!     endif
      call alloc ('cent', cent, 3, ncent)
      if(iperiodic.eq.0 .or. iperiodic.eq.1) write(6,'(/,''center positions'')')
      do 50 ic=1,ncent
        read(5,*) (cent(k,ic),k=1,ndim)
!       if(iperiodic.eq.3) then
!         do 40 k=1,ndim
!  40       cent(k,ic)=cent(k,ic)*alattice
!       endif   
   50   if(iperiodic.eq.0 .or. iperiodic.eq.1) write(6,'(''center'',i4,1x,''('',3f13.5,'')'')') ic,(cent(k,ic),k=1,ndim) ! quick fix for large systems !GO

! Convert center positions from primitive lattice vector units to cartesian coordinates
      if(iperiodic.eq.3) then
        write(6,'(/,''center positions in primitive lattice vector units and in cartesian coordinates'')')
        do 66 ic=1,ncent
          do 62 k=1,ndim
   62       cent_tmp(k)=cent(k,ic)
          do 65 k=1,ndim
            cent(k,ic)=0
            do 65 i=1,ndim
   65         cent(k,ic)=cent(k,ic)+cent_tmp(i)*rlatt(k,i)
   66     write(6,'(''center'',i4,1x,''('',3f9.5,'')'',1x,''('',3f9.5,'')'')') ic, (cent_tmp(k),k=1,ndim),(cent(k,ic),k=1,ndim)
      endif
      write(6,*)

      call object_modified ('cent')

      if(nloc.gt.0) then
        write(6,'(/,''pseudopotential calculation'')')
        call alloc ('vps', vps, nelec, ncent, MPS_L)
        call alloc ('npotd', npotd, nctype)
        if(nloc.eq.1) then
          call readps
        elseif(nloc.eq.2.or.nloc.eq.3) then
          call readps_tm
        elseif(nloc.eq.4 .or. nloc.eq.5) then
          call readps_champ
        elseif(nloc.eq.6) then
          call readps_gauss
        else
          stop 'nloc >= 7'
        endif
        do 67 ict=1,nctype
          if(npotd(ict).ge.4 .and. nquad.lt.12) then
            nquad=12
            write(6,'(''Number of quadrature points, nquad, reset to 12 because npotd='',i2)') npotd(ict)
          endif
          if(npotd(ict).ge.5 .and. nquad.lt.24) then
            nquad=24
            write(6,'(''Number of quadrature points, nquad, reset to 24 because npotd='',i2)') npotd(ict)
          endif
          if(npotd(ict).ge.6) write(6,'(''Warning: We are not ensuring the right number of quadrature points for npotd >=6'')')
   67   continue
        call alloc ('xq0', xq0, MPS_QUAD)
        call alloc ('yq0', yq0, MPS_QUAD)
        call alloc ('zq0', zq0, MPS_QUAD)
        call alloc ('xq', xq, MPS_QUAD)
        call alloc ('yq', yq, MPS_QUAD)
        call alloc ('zq', zq, MPS_QUAD)
        call alloc ('wq', wq, MPS_QUAD)
        call gesqua(nquad,xq0,yq0,zq0,wq)
        if(ipr.ge.0) then
          write(6,'(''Quadrature points for nonlocal pseudopotential'')')
          do 68 i=1,nquad
   68       write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)
        endif
      endif

      if(iperiodic.eq.3) call set_ewald
      if(iperiodic.eq.1) call set_ewald_1d

! Compute total nuclear charge and compare to number of electrons
! Warn if not equal, stop if they differ by more than 2.
      znuc_tot=0
      do 69 ic=1,ncent
        ict=iwctype(ic)
   69   znuc_tot=znuc_tot+znuc(ict)
      if(iperiodic.eq.3) znuc_tot=znuc_tot*vcell_sim/vcell
      if(znuc_tot.ne.dfloat(nelec)) write(6,'(''znuc_tot='',f6.1,'' != nelec='',i4)') znuc_tot,nelec
      if(abs(znuc_tot-dfloat(nelec)).gt.3) stop 'abs(znuc_tot - nelec) > 3' ! JT

      if(nloc.ne.-3) then ! RM
        if(abs(znuc_tot-dfloat(nelec)).gt.6) stop 'abs(znuc_tot - nelec) > 6'
      endif

! TEMPORARY: Warning: we are not calling readforce and only using one geometry
      if(index(mode,'fit').ne.0) then
        nforce=1
        nwftype=1
        call alloc ('iwftype', iwftype, nforce)
        iwftype(1)=1
      endif
      
! Determinantal section
      read(5,*) section
      write(6,'(/,a30,/)') section

      if(ibasis.ge.3 .and. ibasis.le.7) then
        read(5,*) inum_orb
        iorb_used=0
        iorb_format='unused'
      else
        read(5,*) inum_orb,iorb_used,iorb_format
      endif
      write(6,'(''inum_orb,iorb_used,iorb_format ='',t31,i10,i5,1x,a16)') inum_orb,iorb_used,iorb_format
      if(iperiodic.gt.0 .and. (inum_orb.ne.0.and.abs(inum_orb).ne.4.and.abs(inum_orb).ne.5.and.abs(inum_orb).ne.6 &
        &.and.abs(inum_orb).ne.8)) then
         stop 'abs(inum_orb) must be 0, 4 or 5 or 6 or 8'
      endif
! If ndim=2 then ngrid_orbx,ngrid_orby are read in from the orbital file itself
      if(inum_orb .ne. 0 .and. ndim.eq.3) then
        read(5,*)ngrid_orbx,ngrid_orby,ngrid_orbz,igrad_lap ! This grid is not used for smoothing B-splines, it is read in from the bwfn.data file
        if(abs(inum_orb)==4 .or. abs(inum_orb)==5 .or. abs(inum_orb)==8) &
     &  write(6,'(''Number of grid points for interpolating orbitals='',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
         if(abs(inum_orb).eq.8) then
           if(igrad_lap .eq. 0) then
             write(6,'(''Interpolate orbitals, calculate gradient and Laplacian from interpolated orbitals'')')
           elseif(igrad_lap .eq. 1) then
             write(6,'(''Interpolate orbitals and Laplacian, calculate gradient from interpolated orbitals'')')
           elseif(igrad_lap .eq. 2) then
             write(6,'(''Interpolate orbitals, gradient and Laplacian'')')
           else
             stop 'igrad_lap must equal 0, 1 or 2'
           endif
         endif
      endif

      read(5,*) ndet,nbasis,norb
      write(6,'(''no. of determinants ='',t31,i10)') ndet
      write(6,'(''no. of orbitals ='',t31,i10)') norb
      write(6,'(''no. of basis states ='',t31,i10)') nbasis

      call object_modified ('ndet') !JT
      call object_modified ('nbasis') !JT
      call object_modified ('norb') !JT

!     total number of orbitals
      orb_tot_nb = norb
      call object_modified ('orb_tot_nb')

! For Lagrange interpolation allocate orb, dorb and ddorb; for interpolating splines allocate just orb
! For the moment, interpolating splines are being done in the bsplines_mode module, so allocate nothing for them.
!     if(abs(inum_orb).eq.4 .or. abs(inum_orb).eq.8) then
      if(abs(inum_orb).eq.4) then
!       call alloc('orb_num',orb_num,norb,0:ngrid_orbx-1,ngrid_orbx-1,ngrid_orbx-1)
!       call alloc('dorb_num',dorb_num,3,norb,0:ngrid_orbx-1,0:ngrid_orby-1,0:ngrid_orbz-1)
!       call alloc('ddorb_num',ddorb_num,norb,0:ngrid_orbx-1,0:ngrid_orby-1,0:ngrid_orbz-1)
        allocate(orb_num(norb,0:ngrid_orbx-1,0:ngrid_orby-1,0:ngrid_orbz-1),stat=istat)
        if(istat.ne.0) then
          stop "Error in allocating orb_num"
        endif
        if(abs(inum_orb).eq.4) then
          allocate(dorb_num(3,norb,0:ngrid_orbx-1,0:ngrid_orby-1,0:ngrid_orbz-1),stat=istat)
          if(istat.ne.0) then
             stop "Error in allocating dorb_num"
          endif
          allocate(ddorb_num(norb,0:ngrid_orbx-1,0:ngrid_orby-1,0:ngrid_orbz-1),stat=istat)
          if(istat.ne.0) then
             stop "Error in allocating ddorb_num"
          endif
        endif
      endif

      if(norb.lt.nup .or. norb.lt.ndn) stop 'norb must be >= nup and ndn'

      call alloc ('coef', coef, nbasis, orb_tot_nb, nwf)
      call alloc ('zex', zex, nbasis, nwf)

!---swapped
! Call setup routines for Ylm's if we are using recursion to generate them
! irecursion_ylm=0 use Cyrus' and John's spherical harmonics (upto g functions (L=4))
! irecursion_ylm=1 use Ryo' spherical harmonics (any L)
! Note that at present it always calculates upto lmax (set in basis_fns.f) and so it takes long if lmax is large.
! Change it to calculate upto largest l actually used.
      if(ibasis.eq.1) then
        if(nloc.eq.-3) then ! Je spheres
          irecursion_ylm=1
        else
          irecursion_ylm=0
!         irecursion_ylm=1
!         write(6,'(''Warning temporarily set irecursion_ylm=1'')')
        endif
        if(irecursion_ylm.eq.0)then
          write(6,*) 'Not using recursion for Ylm'
        elseif(irecursion_ylm.eq.1) then ! Ryo Maezono's recursion for high order Ylm
          write(6,*) 'Using recursion for Ylm'
          call setup_spherical_harmonics
          call setup_coefficients_ylm
        else
          stop 'irecursion_ylm must be 0 or 1'
        endif
      endif

! Read in analytical or numerical orbitals
      if(ibasis.eq.1) then
        call read_orb_loc
      elseif(ibasis.eq.3.and.numr.eq.1) then
        write(6,'(''Warning: ibasis.eq.3.and.numr.eq.1 never tested'')')
        call read_orb_loc
      elseif((ibasis.ge.3.and.ibasis.le.7).and.numr.eq.0) then
        call read_orb_dot
      elseif(ibasis.ge.4) then
        write(6,'(''read_input: This combination of ibasis='',i2,'' and numr='',i2,'' not allowed'')') ibasis,numr
        write(6,'(''read_input: numr must be 0 for ibasis=4,5,6,7'')')
        stop 'read_input: This combination of ibasis and numr not allowed'
      endif

      if(inum_orb.eq.0) then
        call object_modified ('n_bas') !JT
        call object_modified ('l_bas') !JT
        call object_modified ('m_bas') !JT
        call object_modified ('ictype_basis') !JT
        call object_modified ('zex')   !JT
      endif

! Read in numerical radial basis
! This has to be done after reading in the LCAO coefs. in read_orb_loc because we will use information read in read_orb_loc
! to compactify the radial basis functions in read_bas_num.
      if(inum_orb.eq.0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
        if(minval(zex(:,1)).ne.0.d0) then
          l_purely_analytical_basis = .true.
        else
          l_purely_analytical_basis = .false.
        endif
!       if((ibasis.eq.1.or.ibasis.eq.3).and.numr.gt.0.and.inum_orb.eq.0) call read_bas_num(1)
!       if((ibasis.eq.1.or.ibasis.eq.3) .and. minval(zex(:,1)).eq.0.d0 .and. inum_orb.eq.0) call read_bas_num(1)
        if((ibasis.eq.1.or.ibasis.eq.3) .and. inum_orb.eq.0) call read_bas_num(1)
        if(minval(zex(:,1)).ne.0.d0) write(6,'(/,''Purely analytical radial basis functions used'')')
        if(maxval(zex(:,1)).eq.0.d0) write(6,'(/,''Purely numerical radial basis functions used'')')
        if(minval(zex(:,1)).eq.0.d0 .and. maxval(zex(:,1)).ne.0.d0) &
     &  write(6,'(/,''Mixed analytical and numerical radial basis functions used'')')
        write(6,'(''ict,nrbas_analytical,nrbas_numerical,nrbas='',4i5)') &
     &  (ict,nrbas_analytical(ict),nrbas_numerical(ict),nrbas(ict),ict=1,nctype)
      endif

! get normalization for basis functions
! Devrims note: moved down; basis_norm_dot need idot to be defined
!      if(ibasis.eq.3.and.numr.eq.0) then
!        call basis_norm_dot(1,1)
!       else
!        call basis_norm(1,1)
!      endif

!     read(5,*) (cdet(i,1),i=1,ndet)
!     write(6,'(/,''determinant coefficients'')')
!     write(6,'(20f10.6)') (cdet(k,1),k=1,ndet)

      write(6,'(/,''orbitals in determinants'')')
!     norb_used=0
      call alloc ('iworbd', iworbd, nelec, ndet)
      do 72 i=1,ndet
        read(5,*) (iworbd(j,i),j=1,nelec)
        do 70 j=1,nelec
!          norb_used=max(norb_used,iworbd(j,i))
   70     if(iworbd(j,i).gt.norb) stop 'iworbd(j,i) > norb'
        if(nup+ndn.lt.60) then
          if(ndn.gt.0) then
            write(fmt,'(''('',i2,''i3,3x,'',i2,''i3)'')') nup,ndn
            write(6,fmt) (iworbd(j,i),j=1,nup),(iworbd(j+nup,i),j=1,ndn)
          else
            write(fmt,'(''('',i2,''i3)'')') nup
            write(6,fmt) (iworbd(j,i),j=1,nup)
          endif
         else
          write(6,'(30i4)') (iworbd(j,i),j=1,nup)
          write(6,'(30i4)') (iworbd(j+nup,i),j=1,ndn)
        endif
        do 71 j=2,nup
          do 71 jj=1,j-1
   71       if(iworbd(jj,i).eq.iworbd(j,i)) stop 'An up-spin determinant has 2 identical orbitals'
        do 72 j=2,ndn
          do 72 jj=1,j-1
   72       if(iworbd(jj+nup,i).eq.iworbd(j+nup,i)) stop 'A down-spin determinant has 2 identical orbitals'
!      if(norb_used.lt.norb) then
!        write(6,'(''norb reset from'',i4,'' to'',i4)') norb,norb_used
!        norb=norb_used
!      endif

      read(5,*) ncsf
      write(6,'(/,''ncsf='',i5)') ncsf
      call systemflush(6)

      call alloc ('csf_coef', csf_coef, ncsf, nwf)
      read(5,*) (csf_coef(icsf,1),icsf=1,ncsf)
      write(6,'(''CSF coefs='',20f10.6)') (csf_coef(icsf,1),icsf=1,ncsf)
      call alloc ('ndet_in_csf', ndet_in_csf, ncsf)
      read(5,*) (ndet_in_csf(icsf),icsf=1,ncsf)
      write(6,'(''ndet_in_csf='',20i4)') (ndet_in_csf(icsf),icsf=1,ncsf)
      call alloc ('iflag', iflag, ndet)
      do 75 idet=1,ndet
   75   iflag(idet)=0

! If all the cdet_in_csf are inputted in integer format (no dots in those lines) then csf_coef are
! assumed to correspond to normalized CSF's and the cdet_in_csf are renormalized so that the
! CSF's are normalized.
! normalize_csf is reset to 0 if any of the cdet_in_csf's are in floating format.
      call alloc ('iwdet_in_csf', iwdet_in_csf, maxval(ndet_in_csf), ncsf)
      call alloc ('cdet_in_csf', cdet_in_csf, maxval(ndet_in_csf), ncsf)
      normalize_csf=1
      do 85 icsf=1,ncsf
        read(5,*) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        write(6,'(''CSF'',i4,'' iwdet_in_csf='',100i4)') icsf,(iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        read(5,'(a)') input_line
        if(index(input_line,'.').ne.0) normalize_csf=0
        read(input_line,*) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        write(6,'(''CSF'',i4,'' cdet_in_csf='',900f8.5)') icsf,(cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        do 85 idet_in_csf=1,ndet_in_csf(icsf)
          iflag(iwdet_in_csf(idet_in_csf,icsf))=1
   85     if(iwdet_in_csf(idet_in_csf,icsf).gt.ndet) stop 'iwdet_in_csf(idet_in_csf,icsf) > ndet'

      if(normalize_csf.eq.1) then
! First calculate normalization and adjust csf_coef to correspond to that.
        write(6,'(''Normalizing cdet_in_csf'')')
        do 88 icsf=1,ncsf
          csf_norm=0
          do 86 idet_in_csf=1,ndet_in_csf(icsf)
  86       csf_norm=csf_norm+cdet_in_csf(idet_in_csf,icsf)**2
          csf_norm=sqrt(csf_norm)
          do 88 idet_in_csf=1,ndet_in_csf(icsf)
  88       cdet_in_csf(idet_in_csf,icsf)=cdet_in_csf(idet_in_csf,icsf)/csf_norm
      endif

! Check if all the determinants are used in CSFs
      do 90 idet=1,ndet
   90   if(iflag(idet).eq.0) write(6,'(''Warning: determinant'',i3,'' is unused'')') idet

      call sort_iworbd

      call object_modified ('iworbd')  !JT
      write(6,'(/,''Determine unique up and dn determinants'')')
      call determinant_up_dn

      call object_modified ('ncsf')         !JT
      call object_modified ('csf_coef')     !JT
      call object_modified ('ndet_in_csf')  !JT
      call object_modified ('iwdet_in_csf') !JT
      call object_modified ('cdet_in_csf')  !JT

      if(ndim.eq.2 .and. ibasis.eq.1) then
        read(5,*) ltot
        write(6,'(''L_tot='',i3)') ltot
      endif
      if(ndim.eq.2.and.(ibasis.eq.1.or.ibasis.eq.3).and.inum_orb.eq.0) call emagnetic(ltot)
!     if(ndim.eq.2 .and. (iperiodic.eq.0 .and. nloc.ne.-4)) call emagnetic(ltot)
!     if(ibasis.eq.2) call read_orb_pw_real
      if(ibasis.eq.2) call read_orb_pw
!     if(iperiodic.eq.0 .and. inum_orb.gt.0) call read_orb_num

! Jastrow section
      read(5,*) section
      write(6,'(/,a30,/)') section
      call systemflush(6)

      read(5,*) ianalyt_lap
      write(6,'(''ianalyt_lap='',i3)') ianalyt_lap

      read(5,*) ijas,isc,nspin1,nspin2,nord,ifock
      write(6,'(''ijas,isc,nspin1,nspin2,nord,ifock='',9i4)') &
     &ijas,isc,nspin1,nspin2,nord,ifock
      call systemflush(6)

      if(ianalyt_lap.eq.0 .and. nloc.gt.0) stop 'Cannot have numerical Lap. with pseudopot'
      if(ianalyt_lap.eq.0 .and. iperiodic.gt.1) stop 'Cannot have numerical Lap. with periodic system: distances in jastrow_num not correct'
      if(ianalyt_lap.eq.0 .and. iperiodic.eq.1) write(6,'(''Warning: numerical Lap. might not be correct for iperiodic = 1'')')
      if(ijas.ne.4 .and. iperiodic.gt.0) stop 'Only ijas=4 implemented for periodic systems'
      if(ijas.gt.6) stop 'only ijas=1,2,3,4,5,6 implemented'
      if(ifock.lt.0.or.ifock.gt.4) stop 'ifock must be between 0 and 4'
      if(ndn.eq.1.and.nspin2.eq.3) stop '1 spin down and nspin2=3'
      if((ijas.eq.4.or.ijas.eq.5).and. &
     &(isc.ne.2.and.isc.ne.4.and.isc.ne.6.and.isc.ne.7.and. &
     &isc.ne.8.and.isc.ne.10.and. &
     &isc.ne.12.and.isc.ne.14.and.isc.ne.16.and.isc.ne.17)) &
     & stop 'if ijas=4 or 5, isc must be one of 2,4,6,7,8,10,12,14,16,17'
      if((ijas.eq.6).and.(isc.ne.6.and.isc.ne.7)) stop 'if ijas=6, isc must be 6 or 7'

      if(ijas.eq.3.and.nspin2.gt.1) stop 'ijas=3 and nspin2>1'
      nspin2b=iabs(nspin2)
      nocuspb=0
      if(nspin2.lt.0) then
        if(nspin2.eq.-1) nocuspb=1
        nspin2=1
      endif

      if(ijas.eq.1) then
        call alloc ('cjas1', cjas1, nwf)
        call alloc ('cjas2', cjas2, nwf)
        read(5,*) cjas1(1),cjas2(1)
        write(6,'(''jastrow numerator,denominator ='',2f10.5)') cjas1(1),cjas2(1)
      elseif(ijas.eq.2) then
        nparm_read=69
        if(isc.ge.2) then
          call alloc ('scalek', scalek, nwf)
          read(5,*) scalek(1),a21
          write(6,'(''scalek,a21='',t31,9f10.5)') scalek(1),a21
        endif
        call alloc ('a1', a1, nparm_read, nspin2-nspin1+1, nwf)
        do 270 isp=nspin1,nspin2
          read(5,*) (a1(iparm,isp,1),iparm=1,nparm_read)
          if(ncent.gt.1.and.a1(2,isp,1).ne.zero) &
     &    write(6,'(''WARNING e-n cusp condition cannot be imposed'', &
     &    '' for molecules'',/,''with present weighted form of Jastrow'')')
          write(6,'(''a='',x,7f10.6,(8f10.6))') (a1(iparm,isp,1),iparm=1,nparm_read)
  270   continue
        call alloc ('a2', a2, nparm_read, nspin2-nspin1+1, nwf)
        do 275 isp=nspin1,nspin2
          read(5,*) (a2(iparm,isp,1),iparm=1,nparm_read)
  275     write(6,'(''b='',x,7f10.6,(8f10.6))') (a2(iparm,isp,1),iparm=1,nparm_read)
      elseif(ijas.eq.3) then
        nparm_read=2
        nparmc_read=(nord**3+5*nord)/6+nord**2+nord
        write(6,'(''nparm_read,nparmc_read='',3i5)') nparm_read,nparmc_read
        if(isc.ge.2) then
          call alloc ('scalek', scalek, nwf)
          read(5,*) scalek(1),a21
          write(6,'(''scalek(1),a21='',2f10.5)') scalek(1),a21
        endif
        call alloc ('a', a, nparm_read, nwf)
        read(5,*) (a(iparm,1),iparm=1,nparm_read)
        write(6,'(''a='',x,7f10.6,(8f10.6))')(a(iparm,1),iparm=1,nparm_read)
        call alloc ('b', b, nparm_read, nspin2b-nspin1+1,nwf)
        do 280 isp=nspin1,nspin2b
          read(5,*) (b(iparm,isp,1),iparm=1,nparm_read)
  280     write(6,'(''b='',x,7f10.6,(8f10.6))') (b(iparm,isp,1),iparm=1,nparm_read)
        call alloc ('c', c, nparmc_read, nctype, nwf)
        do 290 it=1,nctype
          read(5,*) (c(iparm,it,1),iparm=1,nparmc_read)
  290     write(6,'(''c='',x,7f10.6,(8f10.6))') (c(iparm,it,1),iparm=1,nparmc_read)
        if(ifock.gt.0) then
          nfock=9
          if(ifock.eq.2) nfock=15
          call alloc ('fck', fck, nfock, nctype, nwf)
          do 300 it=1,nctype
            read(5,*) (fck(iparm,it,1),iparm=1,nfock)
            if(ifock.gt.2) then
              call scale3(1,it)
            endif
            write(6,'(''f='',x,7f10.6,(8f10.6))') (fck(iparm,it,1),iparm=1,nfock)
  300     continue
        endif
      elseif(ijas.ge.4.and.ijas.le.6) then
        if(ifock.gt.0) stop 'fock not yet implemented for ijas=4,5,6'
        read(5,*) norda,nordb,nordc
        write(6,'(''norda,nordb,nordc='',3i5)') norda,nordb,nordc
        nparma_read=2+max(0,norda-1)
        nparmb_read=2+max(0,nordb-1)
        nparmc_read=nterms4(nordc)
        write(6,'(''nparma_read,nparmb_read,nparmc_read='',3i5)') nparma_read,nparmb_read, &
     &  nparmc_read
! WAS
        if(iperiodic.gt.0 .and. nordc.gt.0 .and. ijas .le. 3) stop 'J_een only implemented with ijas= 4,5,6'
!cWAS
        if(isc.ge.2) then
          call alloc ('scalek', scalek, nwf)
          read(5,*) scalek(1),a21
          write(6,'(''scalek(1),a21='',2f10.5)') scalek(1),a21
        endif
        if(isc.ne.8 .and. isc.ne.10) then
          parm2min=-scalek(1)
        else
          parm2min=-1.d0
        endif
        call alloc ('a4', a4, nparma_read, nctype, nwf)
        do 301 it=1,nctype
           read(5,*) (a4(iparm,it,1),iparm=1,nparma_read)
           write(6,'(''a='',x,7f10.6,(8f10.6))') (a4(iparm,it,1),iparm=1,nparma_read)
           if(nparma_read.ge.2 .and. a4(2,it,1).lt.parm2min) then
               write(6,'(''Warning: a4(2,it,1) too low, Jastrow denom could become negative'')')
               stop 'a4(2,it,1) too low, Jastrow denom could become negative'
           endif
  301   continue
        call alloc ('b', b, nparmb_read, nspin2b-nspin1+1,nwf)
        do 302 isp=nspin1,nspin2b
          read(5,*) (b(iparm,isp,1),iparm=1,nparmb_read)
          write(6,'(''b='',x,7f10.6,(8f10.6))') (b(iparm,isp,1),iparm=1,nparmb_read)
           if(nparmb_read.ge.2 .and. b(2,isp,1).lt.parm2min) then
             write(6,'(''Warning: b(2,isp,1) too low, Jastrow denom could become negative'')')
             stop 'b(2,isp,1) too low, Jastrow denom could become negative'
           endif
  302   continue
        call alloc ('c', c, nparmc_read, nctype, nwf)
        do 303 it=1,nctype
          read(5,*) (c(iparm,it,1),iparm=1,nparmc_read)
  303     write(6,'(''c='',x,7f10.6,(8f10.6))') (c(iparm,it,1),iparm=1,nparmc_read)
! Note: Fock terms yet to be put in ijas=4,5,6.
      endif
      call systemflush(6)

      call object_modified ('scalek')
      call object_modified('isc')          !fp
      call object_modified ('norda') !JT
      call object_modified ('nordb') !JT
      call object_modified ('nordc') !JT
      call object_modified ('nparma_read') !JT
      call object_modified ('nparmb_read') !JT
      call object_modified ('nparmc_read') !JT
      call object_modified ('a4') !JT
      call object_modified ('b') !JT
      call object_modified ('c') !JT
      call object_modified ('nspin2b') !JT

! Read cutoff for Jastrow4,5,6 and call set_scale_dist to evaluate constants
! that need to be reset if scalek is being varied.
! If cutjas=0, then reset cutjas_en, cutjas_ee to infinity
! Warning: At present we are assuming that the same scalek is used
! for primary and secondary wavefns.  Otherwise c1_jas6i,c1_jas6,c2_jas6
! should be dimensioned to MWF
      if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) then
        if(iperiodic.eq.1)then  ! set maximum allowed value for cutjas
          cutjas_en = alattice/2.
          cutjas_ee = alattice/2.
        endif
        read(5,*) cutjas_en_tmp,cutjas_ee_tmp
        if(iperiodic.ne.0 .and. cutjas_en_tmp.gt.cutjas_en+eps) then
          write(6,'(''Warning: input cutjas > half shortest primitive cell lattice vector; &
     &    cutjas_en reset from'',f9.5,'' to'',f9.5)') cutjas_en_tmp,cutjas_en
!         write(6,'(''Warning: input cutjas_en > half shortest primitive cell lattice vector; cutjas_en='',d12.5)') cutjas_en_tmp
!         write(6,'(''Warning: cutjas_en will NOT be reset.  I hope you know what you are doing!'')')
!         cutjas_en = cutjas_en_tmp
        else
          if(cutjas_en_tmp.lt.cutjas_en-eps) then
             write(6,'(''Warning: Could use larger cutjas_en='',f9.5, &
     &'' instead of the input value='',f9.5)') cutjas_en,cutjas_en_tmp
          endif
          write(6,'(''input cutjas_en='',d12.5)') cutjas_en_tmp
          cutjas_en=cutjas_en_tmp
        endif
        if(iperiodic.ne.0 .and. cutjas_ee_tmp.gt.cutjas_ee+eps) then
          write(6,'(''Warning: input cutjas > half shortest simulation cell lattice vector; &
     &cutjas_ee reset from'',f9.5,'' to'',f9.5)') cutjas_ee_tmp,cutjas_ee
        else
          if(cutjas_ee_tmp.lt.cutjas_ee-eps) then
            write(6,'(''Warning: Could use larger cutjas_ee='',f9.5, '' instead of the input value='',f9.5)') cutjas_ee,cutjas_ee_tmp
          endif
         write(6,'(''input cutjas_ee='',d12.5)') cutjas_ee_tmp
         cutjas_ee=cutjas_ee_tmp
       endif
       if(cutjas_en_tmp.le.0.d0) then
         write(6,'(''cutjas_en reset to infinity'')')
         cutjas_en=1.d99
       endif
       if(cutjas_ee_tmp.le.0.d0) then
         write(6,'(''cutjas_ee reset to infinity'')')
         cutjas_ee=1.d99
       endif
      endif
      call set_scale_dist(1,1)
      call systemflush(6)

      if(ifock.gt.0) then
! Setup for Chris' Fock
!       fflag=7

! Read pars for Chris's wf
!       call wfpars
        if(ifock.eq.4) then
          open(11, file = '/afs/theory.cornell.edu/user/tc/cyrus/qmc/vmc/lob.dat')
          rewind 11
          nsplin = 1001
          call alloc ('rlobx', rlobx, nsplin)
          call alloc ('rloby', rloby, nsplin)
          read(11,*) (rlobx(i),rloby(i),i=1,nsplin)
          call alloc ('rloby2', rloby2, nsplin)
          call spline(rlobx,rloby,nsplin,0.d0,0.d0,rloby2)
        endif
      endif

! Optional section:
!   default values:
      ifixe=0
      xmax=5.d0
      xfix(1)=0.d0
      xfix(2)=0.d0
      xfix(3)=360.d0 !GO
      rring=0.d0
      iring_coulomb=0
      iperturb=0
      ang_perturb=0.d0
      amp_perturb=0.d0
      shrp_perturb=0.d0
      omg_perturb=0.1d0
      ifourier=0
      izigzag=0
      zzdelyr=0.25
      if(nloc.eq.-4 .or. nloc.eq.-1 .or. nloc.eq.-5) izigzag=1
      fmax1=10.d0
      fmax2=1.d0
      nv=0
      idot=0
      rmin=0.d0
      rmax=10.d0
      nmeshr=NAX
      nmesht=NAX
      nmeshk1=NAK1
      icoosys=1
      iper_gaussian_type = 2
      if (nup.eq.ndn .and. ibasis.ge.3 .and. ibasis.le.7) then
        iantiferromagnetic = 1
      else
        iantiferromagnetic = 0
      endif
      gndot_k=0.d0 !GO
      gauss_width_max = -1 ! GO
!     default values of dot_bump_height and dot_bump_radius are set above
!        where w0, etc... are read in

! Read optional variables if any:
! ADG: used only for 2D systems (dots, rings, composite fermions etc..)
! (also for 2D systems with numerical orbitals)
!     if((ibasis.ge.3 .and. ibasis.le.7) .or. (inum_orb .ne. 0 .and. ndim .eq. 2)) then
!  ACM: I think we should read this in for ALL 2d systems.
      if((ibasis.ge.3 .and. ibasis.le.7) .or. (nloc.eq.-1 .and. ndim .eq. 2)) then
        read(5,*) section
        write(6,'(/,a30,/)') section
        read(5,opt_list)
        write(6,opt_list)
      endif
      call systemflush(6)

      if(nloc.eq.-5) then ! ring with barrier in center
         dot_bump_radius_inv2 = 1.0d0 / (dot_bump_radius * dot_bump_radius)
      endif

! make sure that iantiferromagnetic makes sense
      if(iantiferromagnetic.eq.1) then
        if(nup.ne.ndn) then
          stop 'iantiferromagnetic can be 1 only if nup=ndn'
        else
          do i =1,nup
            if (iworbddn(i,1).ne.(iworbdup(i,1)+1))  then
              stop 'Incorrect input format for iantiferromagnetic=1. Orbitals must alternate &
                 &between spin-up and spin-down (ie, up = 13 5 7..., dn = 24 6 8...)'
            endif
          enddo
        endif
      endif

! pair density calculation parameters:
      if(ifixe.lt.-4 .or. ifixe.gt.nelec) stop 'ifixe must be between -4 and nelec'
      if(abs(ifixe).gt.0) then
        if(ifixe.gt.0 .and. index(mode,'vmc').eq.0) stop 'fixed electron not possible in fit or dmc!'
!        if(ifixe.gt.0 .and. nopt_iter.ne.0) stop 'fixed electron not possible with optimization'
!        if(ncent.ne.1) stop 'Pair-density calculation not implemented for ncent.ne.1'
        if(index(mode,'vmc').ne.0 .and. imetro.ne.1) stop 'Pair-density calculation only possible for imetro=1 in vmc'
        if(index(mode,'dmc').ne.0 .and. abs(idmc).ne.2) stop 'Pair-density calculation only possible for idmc=2 in dmc'
        if(ndim.ne.2) stop 'Pair-density calculation not implemented for 3D systems'
      endif
      if(iperiodic.eq.1) then ! we print out densities for all x
        delxi(1) = (2.*NAX + 1.)/alattice
        delxi(2) = NAX/xmax
      else
        delxi(1)=NAX/xmax
        delxi(2)=delxi(1)
      endif
      if (abs(xfix(1)) .gt. xmax) stop 'abs(xfix(1)) should be less equal than xmax' !GO
      if (abs(xfix(2)) .gt. xmax) stop 'abs(xfix(2)) should be less equal than xmax'
      if ((xfix(1) .eq. 0d0) .and.  (xfix(2) .eq. 0d0) .and. xfix(3) .ne. 360) then
        stop 'xfix(3) must be equal to 360 if xfix(1) and xfix(2) equal to 0.' 
      end if  
      if (xfix(3) .le. 0 .or. xfix(3) .gt. 360) stop 'xfix(3) must be 0 < xfix(3) <= 360' 
      allocate(imeshfix1(0))
      allocate(imeshfix2(0))
      allocate(thetafix(0))
      ! Append first fixed mesh point with angle 0
      call append(imeshfix1, nint(delxi(1) * xfix(1)))
      call append(imeshfix2, nint(delxi(2) * xfix(2)))
      call append(thetafix, 0d0)
      ! Append other mesh points according to desired symmetry
      ithetafix = nint(360 / xfix(3))
      do itheta = 1, ithetafix - 1
        thetafixtemp = pi * (itheta * xfix(3)) / 180 
        call rotate(thetafixtemp, xfix(1), xfix(2), xfix1rot, xfix2rot)
        ixfix1rot = nint(delxi(1) * xfix1rot)
        ixfix2rot = nint(delxi(2) * xfix2rot)
        ! Check if this point inside the fixed mesh point list. If not add to the list.
        insidelist = .false.
        do imesh = 1, size(imeshfix1)
          if (ixfix1rot .eq. imeshfix1(imesh) .and. ixfix2rot .eq. imeshfix2(imesh)) then
            insidelist = .true.
            exit
          end if    
        end do
        if (.not. insidelist) then
          call append(imeshfix1, ixfix1rot)
          call append(imeshfix2, ixfix2rot)
          call append(thetafix, thetafixtemp)
        end if
      end do
      ithetafix = size(thetafix) !GO

! fourier transform :
      if(ifourier.lt.0 .or. ifourier.gt.3 ) stop 'ifourier must be 0,1,2, or 3'
      if(ifourier.gt.1) then
        if(index(mode,'vmc').ne.0 .and. imetro.ne.1) stop 'Fourier transform only possible for imetro=1 in vmc'
        if(index(mode,'dmc').ne.0 .and. (abs(idmc).ne.2 .or. (nloc.ne.-1 .and. nloc.ne.-5))) &
     &    stop 'Fourier transform calculation only possible for idmc=2,nloc=-1 or -5 in dmc'
        if(ndim.ne.2) stop 'Fourier transform not implemented for 3D systems'
      endif
      if(nmeshk1.gt.NAK1.or.nmeshk1.lt.1) stop 'we must have 1<nmeshk1<=NAK1'
      delk1=fmax1/nmeshk1
      delk2=fmax2/NAK2

! ZigZag quantities:
      if(izigzag.lt.0) stop 'izigzag must be greater than or equal to 0'
      if(izigzag.gt.0) then
        if(ndim.ne.2) stop 'ZigZag measurements not implemented in 3D'
        call alloc ('xold_sav', xold_sav, 3, nelec)
        call alloc ('xnew_sav', xnew_sav, 3, nelec)
        call alloc ('zzposold', zzposold, 2, nelec)
        call alloc ('zzposnew', zzposnew, 2, nelec)
        call alloc ('iold_indices', iold_indices, nelec)
        call alloc ('inew_indices', inew_indices, nelec)
        call alloc ('zzcum', zzcum, nzzvars)
        call alloc ('zzsum', zzsum, nzzvars)
        call alloc ('zzcm2', zzcm2, nzzvars)
      endif

! composite fermions:
      lv=0
      emagv=0.d0
      if(idot.lt.0 .or. idot.gt.3) stop 'idot must be 0,1,2 or 3'
      if(idot.gt.0) then
        if(numr.ne.0) write(6,*) 'numerical orbitals not tested with comp fermions'
        if(nv.lt.0) stop 'nv must be zero or positive'
        if(idot.eq.2) then
          write(6,'(''Ignoring determinantal part for Laughlin wave functions'')')
          emaglz=0.d0                           ! no determinantal part for laughlin wfs
          lv=((2*nv+1)*nelec*(nelec-1))/2
        else
          lv=nv*nelec*(nelec-1)
        endif
        emagv=-0.5d0*bext*lv
        if(ibasis.ne.3) stop 'ibasis must be 3 for composite fermions. set nv to zero if not dealing with composite fermions.'
        if(ndim.ne.2) stop 'ndim must be 2 for composite fermions'
        write(6,*) 'mode=',mode
        if(index(mode,'mov1').ne.0) stop '1 electron move not yet implemented for composite fermions'
        write(6,'(''vortices angular momentum, lv ='',t31,i10)') lv
        write(6,'(''vortices angular mom. magnetic energy ='',t31,f10.5)') emagv
      endif
      emag=emaglz+emagsz+emagv
      write(6,'(''emaglz,emagsz,emagv,emag='',9f10.6)') emaglz,emagsz,emagv,emag
      call systemflush(6)

! Quantum rings/wires:
      if(bext.ne.0.d0 .and. rring.ne.0.d0) stop 'Quantum rings in magnetic field not yet implemented'
      if(iring_coulomb.ne.0 .and. rring.eq.0.d0) then
        stop 'Modified coulomb potential for rings (iring_coulomb=1) only possible for quantum rings'
      endif
      if(iperturb.ne.0) then
        if(iperturb.eq.3) then
          shrp_perturb = 15.0d0
          if (nloc.eq.-4) then
            ang_perturb = 2.0d0*dsqrt(amp_perturb)/omg_perturb
          else ! Quantum Ring
            ang_perturb = 2.0d0*dsqrt(amp_perturb)/(omg_perturb*rring)
          endif
          write(6,*) 'iperturb = 3; parabolic perturbation:'
          write(6,*) 'Resetting shrp_perturb = ', shrp_perturb
          write(6,*) 'Resetting ang_perturb = ', ang_perturb
        endif
        if (shrp_perturb.le.0 .or. shrp_perturb.gt.100) stop 'shrp_perturb must be between 0 and 100'
        if(nloc.eq.-4) then  ! Quantum Wires
          if (iperiodic.eq.0) then
            if (ang_perturb.lt.0 .or. ang_perturb.gt.wire_length) stop 'ang_perturb must be between 0 and wire_length'
          else
            if (ang_perturb.lt.0 .or. ang_perturb.gt.alattice) stop 'ang_perturb must be between 0 and alattice'
          endif
        else ! Quantum Ring
          if (rring.eq.0.d0) stop 'Perturbation only possible for quantum rings'
          if (ang_perturb.lt.0 .or. ang_perturb.gt.4*pi) stop 'ang_perturb must be between 0 and 4pi'
        endif
      endif
      call systemflush(6)

! circular coordinates
      if(ibasis.ge.3 .and. ibasis.le.7) then ! not quite the right conditions for circular coordinates, but good enough for now
        if(icoosys.lt.1 .or. icoosys.gt.2) stop 'icoosys must be 1 or 2'
        if(rmax.lt.0.d0 .or. rmin.lt.0.d0 .or. rmax.lt.rmin) stop 'we must have 0<rmin<rmax'
        if(nmeshr.gt.NAX .or. nmesht.gt.NAX .or. nmeshr.lt.1 .or. nmesht.lt.1) stop 'we must have 1<nmeshr<=NAX  and 1<nmesht<=NAX'
        delti=(2*nmesht+1)/(2*pi)
        delradi=(2*nmeshr+1)/(rmax-rmin)
        rmean=(rmin+rmax)*0.5d0
        write(6,'(''Value of r used in densities is relative to rmean ='',f10.6)') rmean
      endif

! get normalization for basis functions
! (used to be up)
      if(inum_orb.eq.0) then
        if(ibasis.eq.3.and.numr.eq.0) then
          call basis_norm_dot(1,1)
        else
          call basis_norm(1,1)
        endif
      endif

! get nuclear potential energy
      call pot_nn(cent,znuc,iwctype,ncent,pecent)
      write(6,'(''pecent='',f14.7)') pecent
      call systemflush(6)

! get interparticle distances
! Don't need to call distances if hpsi always calls it.
!c    call distances(xold,rvec_en,r_en,rvec_ee,r_ee,pe)
!     call distances(xold,pe)

! Find the minimum distance of each electron to any nucleus
!     do 386 i=1,nelec
!       rmino(i)=99.d9
!       do 385 ic=1,ncent
!         if(r_en(i,ic).lt.rmino(i)) then
!           rmino(i)=r_en(i,ic)
!           nearesto(i)=ic
!         endif
! 385     continue
!       do 386  k=1,ndim
! 386     rvmino(k,i)=rvec_en(k,i,nearesto(i))


! Optimization section
      read(5,*) section
      write(6,'(/,a,/)') section
      call systemflush(6)

      read(5,*) nopt_iter,nblk_max,add_diag(1),p_var,tol_energy
      write(6,'(/,''nopt_iter,nblk_max,add_diag(1),p_var,tol_energy='',i4,i8,1p,9d12.4)') &
     &nopt_iter,nblk_max,add_diag(1),p_var,tol_energy
      increase_blocks_limit = nblk_max           !JT
      energy_threshold = tol_energy              !JT
      diag_stab = add_diag(1)                    !JT
      call object_modified ('diag_stab')         !JT
      call object_modified ('energy_threshold')  !JT
      call object_modified ('p_var')             !JT
      if(nopt_iter.ne.0) igradhess=1
      if(nopt_iter.ne.0 .and. nforce.gt.1) stop 'nopt_iter != 0 .and. nforce > 1 not allowed. At present can optim 1 wf only'

!     if(add_diag(1).le.0.d0) stop 'add_diag(1) must be >0'
      if(p_var.lt.0.d0 .or. p_var.gt.1.d0) stop 'p_var must be in [0,1]'
      if(tol_energy.le.0.d0) stop 'tol_energy must be >0'

      if(index(mode,'fit').eq.0 .and. nopt_iter.eq.0) return

      if(nopt_iter.ne.0 .and. (nwf.lt.3)) stop 'if nopt_iter!=0 then nwf should be >=3'

      read(5,*) ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt
      call object_modified ('nparm')
      write(6,'(/,''ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt='' &
     &,6i4,i6,20i4)') ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt

! initialize saved configuration indice iconfg (necessary for some compilers)
      isaved=0
      iconfg=1
! if doing fit, allocate memory for saved configurations
      if(index(mode,'fit').ne.0) then
        call alloc('cvd_sav',cvd_sav,ndim,nelec,ndata)
        call alloc('vd_sav',vd_sav,ndim,nelec,ndata)
        call alloc('psid_sav',psid_sav,ndata)
        call alloc('d2d_sav',d2d_sav,ndata)
        call alloc('div_vd_sav',div_vd_sav,nelec,ndata)
        call alloc('cvk_sav',cvk_sav,ndim,nelec,ndata)
        call alloc('psik_sav',psik_sav,ndata)
        call alloc('div_vk_sav',div_vk_sav,nelec,ndata)
        call alloc('d2k_sav',d2k_sav,ndata)
      endif

!JT      if(mod(iopt,10).ne.2 .and. p_var.ne.0.d0) stop 'For Newton method one can optimize linear combination of energy and variance,

      if(index(mode,'mc').ne.0 .and. nopt_iter.gt.0) then
!        if(mod(iopt,10).eq.1) write(6,'(/,''Optimizing wave function using linear method'',/)')
!        if(mod(iopt,10).eq.2) write(6,'(/,''Optimizing wave function using modified Newton method'',/)')
!        if(mod(iopt,10).eq.3) write(6,'(/,''Optimizing wave function using perturbation theory'',/)')
      elseif(index(mode,'fit').ne.0 .and. (iopt.le.1.or.iopt.ge.3)) then
        iopt=2
        write(6,'(''Warning: iopt set to 2 because now fit uses quench (iopt=2) or possibly Transtrums levmar (iopt=3) only; zxssq is obsolete'')')
      endif
      call systemflush(6)

! JT beg: checking e-N cusp condition on orbitals -> moved to orbitals_mod
!      if (index(mode,'vmc').ne.0 .and. icusp.ge.0) then
!
!      call alloc ('imnbas', imnbas, ncent)
!      imnbas(1)=1
!      do i=1,ncent-1
!        it=iwctype(i)
!        imnbas(i+1)=imnbas(i)+nbasis_ctype(it)
!      enddo
!
!       call ie
!       call cusp_en_orb
!      endif
! JT end

      read(5,*) i3body,irewgt,iaver,istrch
      if(mod(irewgt,100).eq.1) then
        write(6,*) '**Warning irewgt=1 reset to irewgt=10'
        irewgt=irewgt+9
      endif
      write(6,'(''i3body,irewgt,iaver,istrch'',9i5)') i3body,irewgt,iaver,istrch
      call systemflush(6)

      read(5,*) ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt.gt.0 .and. (ijas.ne.2)) &
     &stop 'ipos+...>0 checkjas2 exists only be used with Jastrow2'
!     do 404 i=1,ndata
! 404   wght(i)=one
      write(6,'(''ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt'',10i8)') &
     & ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdt,idbdu

      if(isc.eq.2) write(6,'(''dist scaled r=(1-exp(-scalek*r))/scalek'')')
      if(isc.eq.3) write(6,'(''dist scaled r=(1-exp(-scalek*r-(scalek*r)**2/2))/scalek'')')
      if(isc.eq.4) write(6,'(''dist scaled r=r/(1+scalek*r)'')')
      if(isc.eq.5) write(6,'(''dist scaled r=r/(1+(scalek*r)**2)**.5'')')
      if(isc.eq.8) write(6,'(''dist scaled r=(1-exp(-scalek*r))'')')
      if(isc.eq.10) write(6,'(''dist scaled r=scalek*r/(1+scalek*r)'')')
      call systemflush(6)

      if(ijas.eq.1) write(6,'(''Conventional Jastrow'')')
      if(ijas.eq.2) write(6,'(''Exp. Pade + non-anal terms'')')
      if(ijas.eq.3) write(6,'(''Standard form'')')
      if(ijas.eq.4) write(6,'(''New transferable standard form 4'')')
      if(ijas.eq.5) write(6,'(''New transferable standard form 5'')')
      if(ijas.eq.6) write(6,'(''New transferable standard form 6'')')

      if(icusp.ge.0) write(6,'(''Nuclear cusp constraint is imposed'')')

      call alloc ('lo', lo, norb)
      read(5,*) (lo(iorb),iorb=1,norb)
! JT constuct lo internally instead
!      call object_provide ('lo') !JT
      write(6,'(''lo='',20i3)') (lo(iorb),iorb=1,norb)
!     read(5,*) (n(ib),l(ib),ib=1,nbasis)
!     write(6,'(''n,l='',20(2i3,1x))') (n(ib),l(ib),ib=1,nbasis)
      call systemflush(6)

      if(ijas.le.3) then
        na1=nspin1
        na2=nspin2
       else
        na1=1
        na2=nctype
      endif

      nparmot=0
      iconstrain_gauss_orbs=0
      call alloc ('nparma', nparma, na2-na1+1)
      call alloc ('nparmb', nparmb, nspin2b-nspin1+1)
      call alloc ('nparmc', nparmc, nctype)
      call alloc ('nparmf', nparmf, nctype)
      call alloc ('nparmo', nparmo, notype)
      if(ibasis.le.3) then
        read(5,*) nparml,(nparma(ia),ia=na1,na2), &
     &  (nparmb(isp),isp=nspin1,nspin2b),(nparmc(it),it=1,nctype), &
     &  (nparmf(it),it=1,nctype),nparmcsf,nparms,nparmg
      else
        read(5,*) nparml,(nparma(ia),ia=na1,na2), &
     &  (nparmb(isp),isp=nspin1,nspin2b),(nparmc(it),it=1,nctype), &
     &  (nparmf(it),it=1,nctype),nparmcsf,nparms,nparmg, &
     &  (nparmo(it),it=1,notype)
        do it=1,notype
          nparmot=nparmot+iabs(nparmo(it))
          if(nparmo(it).lt.(1-nbasis) .or. nparmo(it).gt.nbasis) then !GO
            write(6, '(''nparmo must be between (-nbasis+1) and nbasis'')')
            stop 'nparmo must be between (-nbasis+1) and nbasis'
          elseif(nparmo(it).lt.0) then  !some orbitals of type 'it' constrained
            iconstrain_gauss_orbs = 1
            write(6,'(''Constraint imposed in optimization over orbital parameters.'')')
          endif
        enddo
      endif

!     JT: nparmd to replace MPARMD
      nparmd = nparmot+nparmcsf
      call object_modified ('nparmd')

      if(nparmcsf.gt.ncsf) then
        write(6,'(a,i5,a,i5)') 'nparmcsf=',nparmcsf,' must be <= ncsf=',ncsf
        stop 'nparmcsf must be <= ncsf'
      endif
      if(nparmcsf.eq.ncsf) then
        write(6,'(a,i5,a,i5)') 'Warning: since normalization of wavefn. is arb. nparmcsf=',nparmcsf,' should be <= ncsf-1=',ncsf-1
      endif
      call systemflush(6)

      if(ijas.ge.4.and.ijas.le.6) then
        do 405 it=1,nctype
          if(nloc.eq.0) then
!     All-electron with analytic slater basis
            if((norda.eq.0.and.nparma(it).gt.0) &
     &      .or.(norda.gt.0 .and. nparma(it).gt.norda+1)) then
              write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
              stop 'nparma too large for norda in all-electron calculation'
            endif
           else
! Pseudopotential with numerical basis (cannot vary a(1) or a(2)
            if(norda.eq.1) stop 'makes no sense to have norda=1 for numr>0'
            if((norda.eq.0.and.nparma(it).gt.0) &
     &      .or.(norda.gt.0 .and. nparma(it).gt.norda-1)) then
              write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
              stop 'nparma too large for norda in pseudopot calculation'
            endif
          endif
          if(isc.le.10 .and. &
     &       ((nordc.le.2.and.nparmc(it).gt.0) &
     &    .or.(nordc.eq.3.and.nparmc(it).gt.2) &
     &    .or.(nordc.eq.4.and.nparmc(it).gt.7) &
     &    .or.(nordc.eq.5.and.nparmc(it).gt.15) &
     &    .or.(nordc.eq.6.and.nparmc(it).gt.27) &
     &    .or.(nordc.eq.7.and.nparmc(it).gt.43))) then
            write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
            stop 'nparmc too large for nordc in J_een with cusp conds'
          endif
          if(isc.gt.10 .and. &
     &       ((nordc.le.1.and.nparmc(it).gt.0) &
     &    .or.(nordc.eq.2.and.nparmc(it).gt.2) &
     &    .or.(nordc.eq.3.and.nparmc(it).gt.6) &
     &    .or.(nordc.eq.4.and.nparmc(it).gt.13) &
     &    .or.(nordc.eq.5.and.nparmc(it).gt.23) &
     &    .or.(nordc.eq.6.and.nparmc(it).gt.37) &
     &    .or.(nordc.eq.7.and.nparmc(it).gt.55))) then
            write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
            stop 'nparmc too large for nordc without cusp conds'
          endif
  405   continue
! For the b coefs. we assume that b(1) is fixed by the cusp-cond.
        do 406 isp=1,nspin1,nspin2b
            if((nordb.eq.0.and.nparmb(isp).gt.0).or.(nordb.gt.0 .and. nparmb(isp).gt.nordb)) then
              write(6,'(''isp,nordb,nparmb(isp)'',3i5)') isp,nordb,nparmb(isp)
              stop 'nparmb too large for nordb'
            endif
  406   continue
      endif

! compute nparmj and nparme
      nparmj=0
      call alloc ('npoint', npoint, nctype)
      call alloc ('npointa', npointa, na2)
      npointa(1)=0
      do 407 ia=na1,na2
        if(ia.gt.1) npointa(ia)=npointa(ia-1)+nparma(ia-1)
  407   nparmj=nparmj+nparma(ia)
      do 408 isp=nspin1,nspin2b
  408   nparmj=nparmj+nparmb(isp)
      npoint(1)=nparmj
      do 409 it=1,nctype
        if(it.gt.1) npoint(it)=npoint(it-1)+nparmc(it-1)
  409   nparmj=nparmj+nparmc(it)+nparmf(it)
      nparme=nparm-nparml-nparmj-nparmcsf-nparms-nparmg-nparmot
      write(6,'(''No of linear coefs, exponents, Jastrow, det, scale parms varied='',9i5)') &
     &nparml, nparme, nparmj, nparmcsf, nparms
      if(nparme.lt.0) stop 'nparme < 0'
      if(nparme.gt.nbasis) stop 'nparme > nbasis'
      if(nparme.gt.0 .and. numr.gt.0) stop 'nparme > 0 and numr > 0'
      if(nparme.gt.0 .and. ibasis.eq.3 .and. idot.ne.0) stop 'for quantum dots, nparme.gt.0 only possible for Fock-Darwin states'
      if(nparme.gt.0 .and. ibasis.eq.4) stop 'nparme > 0' !GO
      if(nparml.lt.0 .or. nparmj.lt.0 .or. nparmcsf.lt.0 .or. nparms.lt.0 .or.nparmg.lt.0) stop 'nparm? must be >= 0'
      if(nparms.gt.1) stop 'nparms must be 0 or 1'
      nparmjs=nparmj+nparms !JT
      call object_modified ('nparmjs')
!JT      if(nparmcsf.ge.ncsf) then
!JT        write(6,'(''Since normalization of wavefunction is arbitrary, nparmcsf must be <= ncsf-1'')')
!JT        stop 'Since normalization of wavefunction is arbitrary, nparmcsf must be <= ncsf-1'
!JT      endif
      call systemflush(6)

      call alloc ('iwo', iwo, nbasis, notype) !GO
      do it=1,notype
        read(5,*) (iwo(iparm,it),iparm=1,iabs(nparmo(it)))
        if(nparmo(it).lt.0) then  !constrained orbital optimization
          write(6,'(''Constraints applied to orbital parameters - constraints printed below'')')
        endif
        write(6,'(''orbital parameters varied='',10(2i3,2x))')(iwo(iparm,it),iparm=1,iabs(nparmo(it)))
        do iparm=1,iabs(nparmo(it))
          if(iwo(iparm,it).lt.0 .or. iwo(iparm,it).gt.nbasis) then
            stop 'Incorrect value for iwo.'
          endif
        enddo
      enddo
      call systemflush(6)

      call alloc ('iworb', iworb, nparml)
      call alloc ('iwbasi', iwbasi, nparml)
      read(5,*) (iworb(iparm),iwbasi(iparm),iparm=1,nparml)
      write(6,'(''lin. coefs. of orbs varied='',/,10(2i3,2x))') (iworb(iparm),iwbasi(iparm),iparm=1,nparml)
      call systemflush(6)

      call alloc ('iwbase', iwbase, nparme)
      read(5,*) (iwbase(iparm),iparm=1,nparme)
      write(6,'(''exponents varied='',20i3)') (iwbase(iparm),iparm=1,nparme)

!     read(5,*) (iwdet(iparm),iparm=1,nparmd)
!     write(6,'(''determinantal coefs varied='',20i3)')
!    &(iwdet(iparm),iparm=1,nparmd)

      call alloc ('iwcsf', iwcsf, nparmcsf)
      read(5,*) (iwcsf(iparm),iparm=1,nparmcsf)
      write(6,'(''CSF coefs varied='',100i3)') (iwcsf(iparm),iparm=1,nparmcsf)
      do 412 iparm=1,nparmcsf
  412   if(iwcsf(iparm).gt.ncsf) stop 'iwcsf(iparm).gt.ncsf'
      call systemflush(6)

      call object_modified ('iwcsf')

      if(ijas.eq.2.or.ijas.eq.3) then
        write(6,'(''Correl. params. that are varied are:'')')
        call alloc ('iwjasa', iwjasa, nparmj, nspin2-nspin1+1)
        do 414 isp=nspin1,nspin2
          read(5,*) (iwjasa(iparm,isp),iparm=1,nparma(isp))
  414     write(6,'(''a: '',30i3)') (iwjasa(iparm,isp),iparm=1,nparma(isp))
        call alloc ('iwjasb', iwjasb, nparmj, nspin2b-nspin1+1)
        do 416 isp=nspin1,nspin2b
          read(5,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
  416     write(6,'(''b: '',30i3)') (iwjasb(iparm,isp),iparm=1,nparmb(isp))
       elseif(ijas.ge.4.and.ijas.le.6) then
        call alloc ('iwjasa', iwjasa, nparmj, nctype)
        do 418 it=1,nctype
          read(5,*) (iwjasa(iparm,it),iparm=1,nparma(it))
  418     write(6,'(''a: '',30i3)') (iwjasa(iparm,it),iparm=1,nparma(it))
        call alloc ('iwjasb', iwjasb, nparmj, nspin2b-nspin1+1)
        do 420 isp=nspin1,nspin2b
          read(5,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
  420     write(6,'(''b: '',30i3)') (iwjasb(iparm,isp),iparm=1,nparmb(isp))
      endif
      if(ijas.ge.3.and.ijas.le.6) then
        call alloc ('iwjasc', iwjasc, nparmj, nctype)
        do 425 it=1,nctype
          read(5,*) (iwjasc(iparm,it),iparm=1,nparmc(it))
  425     write(6,'(''c: '',60i3)') (iwjasc(iparm,it),iparm=1,nparmc(it))
        if(ifock.gt.0) then
          do 430 it=1,nctype
            call alloc ('iwjasf', iwjasf, 15, nctype)
            read(5,*) (iwjasf(iparm,it),iparm=1,nparmf(it))
  430       write(6,'(''f: '',30i3)') (iwjasf(iparm,it),iparm=1,nparmf(it))
        endif
      endif
      call systemflush(6)

      if(icusp2.ge.1 .and. ijas.eq.3 .and. isc.le.7) call cuspinit3(1)
      if(icusp2.ge.1 .and. ijas.eq.4 .and. isc.le.10) call cuspinit4(0)

! Moved here from fit

      write(6,'(''nparm,nparml,nparmj,nparmcsf,nparms,nparmg,nparme='',9i5)') nparm,nparml,nparmj,nparmcsf,nparms,nparmg,nparme
      write(6,'(''total number of parameters, nparm='',i5)') nparm

      if(mbasis_ctype.gt.0) then   ! Needed since nbasis_ctype(it) is not initialized for quantum rings
        call alloc ('imnbas', imnbas, ncent)
        imnbas(1)=1
        do i=1,ncent-1
          it=iwctype(i)
          imnbas(i+1)=imnbas(i)+nbasis_ctype(it)
        enddo
      endif

! If necn < 0, fit will call orb_params to figure out which orbital (LCAO) coefs are varied (iworb, iwbasi), and which are constrained to be equal (ieorb, iebasi)
      read(5,*) necn,nebase
      write(6,'(/,''No of linear coefs, exponents set equal='',3i5)') necn,nebase
      write(6,'(''If necn<0, orb_params will be called, if nebase<0, basis_exp_params will be called'',/)')

      call alloc ('ieorb', ieorb, 2, max(0,necn))
      call alloc ('iebasi', iebasi, 2, max(0,necn))
      read(5,*) ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,max(0,necn))

      do j=1,max(0,necn)
        do i=1,2
          if(ieorb(i,j).gt.norb) stop 'ieorb(i,j).gt.norb'
          if(iebasi(i,j).gt.nbasis) stop 'iebasi(i,j).gt.nbasis'
        enddo
      enddo

      write(6,'(''lin. coefs of orbs (ieorb,iebasi) set equal='',/,5(2(2i3,2x),2x))') ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)

      call alloc ('iebase', iebase, 2, nbasis)
      read(5,*) ((iebase(i,j),i=1,2),j=1,max(0,nebase))
      write(6,'(''expon. (iebase) set equal='',/,10(2i3,2x))') ((iebase(i,j),i=1,2),j=1,max(0,nebase))

      do j=1,max(0,nebase)
        do i=1,2
          if(iebase(i,j).gt.nbasis) stop 'iebase(i,j).gt.nbasis'
        enddo
      enddo

      call systemflush(6)

      call alloc('ipivot', ipivot, norb)
      read(5,*) (ipivot(j),j=1,norb)
      write(6,'(''ipivot='',10i4)') (ipivot(j),j=1,norb)

      read(5,*) eguess
      write(6,'(''eguess='',f12.6)') eguess

      read(5,*) pmarquardt,tau_fit,noutput,nstep_fit,ibold
      write(6,'(''pmarquardt,tau_fit,ibold '',2g9.2,i3)') pmarquardt,tau_fit,ibold

      read(5,'(2l2)') analytic_jacobian,cholesky
      write(6,'(''analytic_jacobian,cholesky'',2l2)') analytic_jacobian,cholesky

      if(analytic_jacobian .and. irewgt.ne.0) write(6,'(''*** Warning *** analytic derivs not correctly implemented when points are reweighted'')')
      if(analytic_jacobian .and. nparmf(1).gt.0) stop 'analytic derivatives not implemented yet for Fock parameters'
! End moved here from fit

      if(iconstrain_gauss_orbs.eq.1) then  !read in constraints for orbital optimization
        call alloc ('norb_constraints', norb_constraints, notype)
        read(5,*) (norb_constraints(it),it=1,notype)
        write(6,'(''Number of constraints applied to each type of orbital: '',4(i4,1x))') (norb_constraints(it),it=1,notype) !GO
        call alloc ('orb_constraints', orb_constraints, notype, nbasis-1, 2) !GO
        do it=1,notype  ! read in constraints
          if(norb_constraints(it).lt.0 .or. norb_constraints(it).gt.(nbasis-1)) then !GO
            write(6, '(''There must be between 0 and (nbasis-1) constraints'')')
            stop 'Invalid number of constraints.'
          endif
          if(norb_constraints(it).eq.0 .and. nparmo(it).lt.0) then
            write(6, '(''nparmo<0, but no constraints found'')')
            stop 'Invalid number of constraints.'
          elseif(norb_constraints(it).ne.0 .and. nparmo(it).ge.0) then
            write(6, '(''Cannot specify constraints if nparmo>=0'')')
            stop 'Invalid number of constraints.'
          endif
          read(5,*) ((orb_constraints(it,i,j),j=1,2),i=1,norb_constraints(it))
          if (norb_constraints(it).eq.0) cycle
          if(ibasis.eq.4) then
            if(it.eq.1) write(6,'(''Applying constraints to floating gaussian x-positions:'')')
            if(it.eq.2) write(6,'(''Applying constraints to floating gaussian y-positions:'')')
            if(it.eq.3) write(6,'(''Applying constraints to floating gaussian widths:'')')
          elseif(ibasis.eq.5) then
            if(it.eq.1) write(6,'(''Applying constraints to floating gaussian radial positions:'')')
            if(it.eq.2) write(6,'(''Applying constraints to floating gaussian angular positions:'')')
            if(it.eq.3) write(6,'(''Applying constraints to floating gaussian radial widths:'')')
            if(it.eq.4) write(6,'(''Applying constraints to floating gaussian angular widths:'')')
          elseif(ibasis.eq.6 .or. ibasis.eq.7) then
            if(it.eq.1) write(6,'(''Applying constraints to floating gaussian x-positions:'')')
            if(it.eq.2) write(6,'(''Applying constraints to floating gaussian y-positions:'')')
            if(it.eq.3) write(6,'(''Applying constraints to floating gaussian x-widths:'')')
            if(it.eq.4) write(6,'(''Applying constraints to floating gaussian y-widths:'')')
          endif
          do icon=1,norb_constraints(it)  ! check that constraints are ok
            write(6,'(''Constraining orbitals: '',2i5)') (orb_constraints(it,icon,j),j=1,2)
            if (orb_constraints(it,icon,1).le.0 .or. orb_constraints(it,icon,1).gt.nbasis) then !GO
              write(6,'(''Constrained orbital must be between 1 and nbasis'')')
              stop 'Constrained orbital out of range'
            endif
            is_first_ok = 1 ! check to see if first orbital of pair is one that we aren't optimizing
            is_second_ok = 0 ! check that second orbital is one that we are optimizing
            do iparm=1,iabs(nparmo(it))
              if (iabs(orb_constraints(it,icon,2)).eq.iwo(iparm,it)) then
                is_second_ok = 1
              endif
              if (orb_constraints(it,icon,1).eq.iwo(iparm,it)) then
                is_first_ok = 0
                exit
              endif
            enddo
            if (is_first_ok.eq.0 .or. is_second_ok.eq.0) then
              write(6,'(''Constraints must be between an orbital that is listed as being optimized and one that is not'')')
              stop 'Invalid constraint'
            endif
            do icon2=icon+1,norb_constraints(it) ! check for duplicate constraints
              if(orb_constraints(it,icon,1).eq.orb_constraints(it,icon2,1)) then
                write(6,'(''Duplicate or Conflicting Constraint'')')
                stop 'Duplicate or conflicting constraint'
              endif
            enddo
            consgn = real(sign(1, orb_constraints(it,icon,2)))
            oparm(it,orb_constraints(it,icon,1),1) = consgn*oparm(it,iabs(orb_constraints(it,icon,2)),1)
          enddo  ! finished check do icon=1,norbconstrain(it)

!         now we make sure all parameters of this type are the same
          write(6,'(''Read-in values reset to enforce constraint:'')')

          if(ibasis.eq.4) then
            if(it.eq.1) write(6,'(''New (constrained) floating gaussian x-positions:'')')
            if(it.eq.2) write(6,'(''New (constrained) floating gaussian y-positions:'')')
            if(it.eq.3) write(6,'(''New (constrained) floating gaussian widths:'')')
          elseif(ibasis.eq.5) then
            if(it.eq.1) write(6,'(''New (constrained) floating gaussian radial positions:'')')
            if(it.eq.2) write(6,'(''New (constrained) floating gaussian angular positions:'')')
            if(it.eq.3) write(6,'(''New (constrained) floating gaussian radial widths:'')')
            if(it.eq.4) write(6,'(''New (constrained) floating gaussian angular widths:'')')
          elseif(ibasis.eq.6 .or. ibasis.eq.7) then
            if(it.eq.1) write(6,'(''New (constrained) floating gaussian x-positions:'')')
            if(it.eq.2) write(6,'(''New (constrained) floating gaussian y-positions:'')')
            if(it.eq.3) write(6,'(''New (constrained) floating gaussian x-widths:'')')
            if(it.eq.4) write(6,'(''New (constrained) floating gaussian y-widths:'')')
          endif
          write(6,'(1000f12.6)') (oparm(it,ib,1),ib=1,nbasis)
        enddo ! do it=1,notype
! if antiferromagnetic constraint is desired, enforce it here.
        if (iantiferromagnetic.eq.1) then
          call sort_af_gauss_orbs(1)
          write(6,'(''Constrained values adjusted to enforce antiferromagnetic arrangement'')')
          if(ibasis.eq.4 .or. ibasis.eq.6 .or. ibasis.eq.7) then
            write(6,'(''Adjusted (antiferromagnetic) floating gaussian x-positions:'')')
            it_af = 1
          elseif(ibasis.eq.5) then
            write(6,'(''Adjusted (antiferromagnetic) floating gaussian angular positions:'')')
            it_af = 2
          endif
          write(6,'('' spin-up: '', 1000f12.6)') (oparm(it_af,iworbdup(ib,1),1),ib=1,nup)
          write(6,'('' spin-down: '', 1000f12.6)') (oparm(it_af,iworbddn(ib,1),1),ib=1,ndn)

        endif
      endif
      
      ! the default value of oparm3_max is the maximum value of a real*8 variable
      ! so that other basis sets using oparm(3,:,:) would not be affected during the optimization.
      oparm3_max = huge(oparm3_max) ! GO
      if (ibasis .eq. 4 .and. gauss_width_max .ge. 0d0) then
        oparm3_max = gauss_width_max
      endif
      write(6,'(/,''oparm3_max = '',G20.8E3,/)') oparm3_max

      write(6,'(''ipr in read_input'',i5)') ipr
      call systemflush(6)

      call object_modified ('nparma') !JT
      call object_modified ('nparmb') !JT
      call object_modified ('nparmc') !JT
      call object_modified ('iwjasa') !JT
      call object_modified ('iwjasb') !JT
      call object_modified ('iwjasc') !JT
      call object_modified ('nparmj') !JT
      call object_modified ('nparmcsf') !JT
      call object_modified ('norb_constraints')
      call object_modified ('orb_constraints')

      return
      end subroutine read_input
!-----------------------------------------------------------------------

      subroutine sort_iworbd
! Written by Cyrus Umrigar
! Order iworbd for each determinant to be monotonically increasing for up and dn electrons separately
! and change signs of cdet_in_csf accordingly.  This is needed for orbital optimization.

      use dorb_mod
      use dets_mod
      implicit real*8(a-h,o-z)


      dimension iodd_permut(ndet)

      do 20 i=1,ndet
        iodd_permut(i)=1
        do 10 j=1,nup
          do 10 k=j+1,nup
            if(iworbd(k,i).lt.iworbd(j,i)) then
              itmp=iworbd(j,i)
              iworbd(j,i)=iworbd(k,i)
              iworbd(k,i)=itmp
              iodd_permut(i)=-iodd_permut(i)
            endif
   10   continue
        do 20 j=nup+1,nup+ndn
          do 20 k=j+1,nup+ndn
            if(iworbd(k,i).lt.iworbd(j,i)) then
              itmp=iworbd(j,i)
              iworbd(j,i)=iworbd(k,i)
              iworbd(k,i)=itmp
              iodd_permut(i)=-iodd_permut(i)
            endif
   20 continue

      do 30 icsf=1,ncsf
        do 30 idet_in_csf=1,ndet_in_csf(icsf)
   30     cdet_in_csf(idet_in_csf,icsf)=iodd_permut(iwdet_in_csf(idet_in_csf,icsf))*cdet_in_csf(idet_in_csf,icsf)

      return
      end subroutine sort_iworbd
!-----------------------------------------------------------------------

      subroutine sort_af_gauss_orbs(iadd_diag)
! Written by Abhijit Mehta
! Order orbitals so that they alternate between spin-up and spin-down
!   (i.e., make floating gaussians antiferromagnetic)
!  This should only work if nup=ndn and iantiferromagnetic=1
!   We use the positions in oparm(it,i_orbital,iadd_diag)
!   So, set iadd_diag = 1 by default
!  Note that for rings, this changes angles so that they are in the range 0.. 2pi rather than -pi..pi
!  We're not worrying about sign of permutation when we reorder orbitals since we
!      don't use CSF's

      use dorb_mod
      use dets_mod
      use orbpar_mod
      use const_mod
      use contrl_per_mod
      use optimo_mod
      implicit real*8(a-h,o-z)

      dimension oparmtemp(notype)   ! used for swap

      if(nup.ne.ndn) then
        write(6,'(''sort_af_gauss_orbs only defined for nup=ndn'')')
        stop 'nup \= ndn in sort_af_gauss_orbs'
      endif

      if(ibasis.eq.5) then  ! rings, so 2nd coordinate is angular position
        it = 2
      else ! wires, so 1st coordinate is x-position (ie, position along length of wire)
        it = 1
      endif

!      do idet=1,ndet
!     Use Shell sort to put orbitals in order of position along wire/ring
!       Adapted from routine written by Cyrus in December 1983
        M=nup+ndn
        LOGNB2=INT(DLOG(DFLOAT(M))/DLOG(2.D0)+1.D-14)
        if(ibasis.eq.5) then
          do ib=1,M
             oparm(it,ib,iadd_diag) = modulo(oparm(it,ib,iadd_diag), 2*pi)
          enddo
        endif
        DO 20 NN=1,LOGNB2
         M=M/2
         K=nup+ndn-M
         DO 20 J=1,K
           DO 10 I=J,1,-M
             L=I+M
             oparmL = oparm(it,L,iadd_diag)
             oparmI = oparm(it,I,iadd_diag)
             IF (oparmL.GT.oparmI)   GOTO 20
             oparmtemp(:) = oparm(:,I,iadd_diag)
             oparm(:,I,iadd_diag) = oparm(:,L,iadd_diag)
             oparm(:,L,iadd_diag) = oparmtemp(:)
!      looping method
!             do it_temp = 1,notype
!               oparmtemp(it_temp) = oparm(it_temp,I,iadd_diag)
!               oparm(it_temp,I,iadd_diag) = oparm(it_temp,L,iadd_diag)
!               oparm(it_temp,L,iadd_diag) = oparmtemp(it_temp)
!             enddo

!      Old code to swap indices instead of values -  ACM
!              itemp=iworbd(I,idet)
!              iworbd(I,idet)=iworbd(L,idet)
!              iworbd(L,idet)=itemp
   10      CONTINUE
   20   CONTINUE

      write(6,'(''sort_af_gauss_orbs: orbs have order:'', 100g13.6)') (oparm(it,ib,iadd_diag),ib=1,(nup+ndn))  ! ACM debug

!     Now make sure that orbitals alternate between up and down - only needed if swapping indices
!        do iorb=1,nup
!          iworbdup(iorb,idet) = iworbd((2*iorb-1), idet)
!          iworbddn(iorb,idet) = iworbd(2*iorb, idet)
!        enddo
!        do iorb=1,nup
!          iworbd(iorb,idet) = iworbdup(iorb,idet)
!          iworbd(iorb+nup, idet) = iworbddn(iorb,idet)
!        enddo
!        write(6,'(a)') 'After sort_af_gauss_orbs, spin-up determinants have orbitals:'
!        write(6,'(a,i5,a,100i4)') ' det # ',idetup, ': ',(iworbdup(iup,idet),iup=1,nup)
!        write(6,'(a)') 'After sort_af_gauss_orbs, spin-down determinants have orbitals:'
!        write(6,'(a,i5,a,100i4)') ' det # ',idetdn, ': ',(iworbddn(idn,idet),idn=1,ndn)
!      enddo

!     Make sure iworbd is properly sorted, then make sure iworbdup and iworbddn are too
!      call sort_iworbd
!
!      do idet=1,ndet
!        do iorb=1,nup
!          iworbdup(iorb,idet) = iworbd(iorb,idet)
!          iworbddn(iorb,idet) = iworbd(iorb+nup, idet)
!        enddo
!      enddo
!
      return
      end subroutine sort_af_gauss_orbs
