MODULE module_simradar_tools
#ifdef SIMRADAR
! Module wiuth the tools for simrad.AR
!   L. Fita, CIMA-IFAECI, November 2025

  USE module_wrf_error

!!!!!!! Subroutines / functions
! beam_height: Function to compute the height of the beam
! calc_tvelocities_1D: Subroutine to compute terminal velocities for water species 1D version
! calc_vt: Function to compute terminal velocity after a gamma function
! cvals_varn: Function to provide the name of the variable accordingly its index from cvals array
! com_gamma: Gamma function
! compute_radar_lonlat_radial: Subroutine to copmute lon, lat, hgt coordinates in the radar matrix 
!   following polar coordinates
! create_logarithmic_scaled_vectorsD: Function to create logarithmic-based equi-spaced scaled double 
!   Real vectors
! distance_between_points: Subroutine tro compute the coordinates at a given distance and angle from a 
!   given point
! DIRCT1: Solution of the geodetic direct problem after T. Vincenty
! inverse: Subroutine to compute distance between 2 points
! is_inside: Subroutine to check if a tile gets a radar
! lin_interpolate_1RR: Function to provide linear interpolation a Real location on a given Real axis
! pertinence_polar: Subroutine to compute the 2-dimensional pertinence of a lonlat polar radial space 
!   (r) into another regular lonlat (WRF)
! Polarimetrics_fromLUT_1R: Subroutine to interpolate the polarimetic characteristics of a given water 
!   specie following T-matrix
! position_from_point: Subroutine tro compute the coordinates at a given distance and angle from a 
!   given point
! read_lut: Subroutine to read the T-MATRIX look-up-table (LUT) from an specific microphysics from a 
!   netCDF file
! radar_inside: Subroutine to determine if a radar fells inside the tile
! real_distances: Subroutine to compute the diferential of coordinattes of a distance from a given 
!   point
! real_distances_ellipsoide: Subroutine to compute the diferential of coordinattes of a distance from 
!   a given point on an ellipsoide
! weights_9pts_passed1coords_ellipsoide_1D: Subroutine to compute the weights of a inverse distance 
!   weighted interpolation using the 3x3 grid points when already passing the coordinates of the 
!   target space at each reference point on an ellipsoide
! weights_hcoordsR27around_polar_ellipsoide_1D:  Subroutine to compute weights of the 27 around grid 
!   points to 3-dimensional interpolate to a new space using the list of cells on the horizontal of 
!   the targeted space that lay within each grid cell of the reference space with the values to 
!   interpolate version on an ellipsoide

! basic: (from PyNCplot)
! ErrMsg: Subroutine to stop execution and provide an error message
! ErrMsg[I/R/S]Avail: Subroutine to stop execution by non exsitent Integer/Real/String value from 
!   available ones and provide an error message
! FileExist: Checks if a file exists
! freeunit: provides the number of a free unit in which open a file
! [I/R/D]toS: Function to transform an Integer/Real/double Real to String
! removeNONnum: Subroutine to remove non numeric characters from a string
! set_dbg: Function to setup debug value
! Sto[I/R/L]: Function to transform a String to an Integer/Real/Boolean

! netCDF: (from PyNCplot)
! close_NCfile: Subroutine to close a netCDF file
! get_dimlength_ncid: Subroutine to get the length of a dimension from an existing to an open netCDF
! get_dims_ncid: Subroutine to get the dimensions from an open netCDF file
! get_Ndims_ncid: Subroutine to get the amount of dimensions from an open netCDF file
! get_Nvars_ncid: Subroutine to get the amount of variables from an open netCDF file
! get_vars_ncid: Subroutine to provide all the ids and names of the variables existing to an open 
!   netCDF
! get_varR[0/1]D_ncid: Subroutine to get a [0/1]D Real variable from a netCDF file unit passing 
!   starting and ending
! get_varD3D_ncid: Subroutine to get a 3D double Real variable from a netCDF file unit passing 
!   starting and ending
! get_varS0D_ncid: Subroutine to get a 0D String variable from a netCDF file unit passing starting 
!   and ending
! handle_errf: Subroutine to provide the error message when something with netCDF went wrong 
!   (including fname)
! isin_ncunit: Function to tell if a given variable is inside a netcdf file unit
! isin_dim_ncid: Function to tell if a given dimension is inside a netcdf file unit
! open_NCfile: Subroutine to open a basic netCDF file

  IMPLICIT NONE
  
  ! Integer 4 bytes 
  INTEGER, PARAMETER                                     :: ii4 = SELECTED_INT_KIND(4)
  ! REal 16 bytes 
  INTEGER, PARAMETER                                     :: rr8 = SELECTED_REAL_KIND(8)
  
  REAL, PARAMETER                                        :: pict = ACOS(-1.)  
  REAL, PARAMETER                                        :: DegRad = 180./pict
  ! double  Real
  REAL(kind=rr8), PARAMETER                              :: dpict = 4.d0*DATAN(1.d0)
  REAL(kind=rr8), PARAMETER                              :: drad = 180.d0/dpict

  ! 16 charcater String
  INTEGER, PARAMETER                                     :: S16 = 16
  ! 32 character String
  INTEGER, PARAMETER                                     :: S32 = 32
  ! 50 character String
  INTEGER, PARAMETER                                     :: Sl = 50
  ! 64 character String
  INTEGER, PARAMETER                                     :: S64 = 64
  ! 128 character String
  INTEGER, PARAMETER                                     :: S128 = 128
  ! 256 character String
  INTEGER, PARAMETER                                     :: S256 = 256
  ! 512 character String
  INTEGER, PARAMETER                                     :: S512 = 512
  ! 1024 character String
  INTEGER, PARAMETER                                     :: S1024 = 1024
  ! 2048 character String
  INTEGER, PARAMETER                                     :: S2048 = 2048

  ! Fill values
  !  !! 'B', NF_BYTE
  !  BYTE, PARAMETER                                      :: fillValueB = 255
  !! 'C', NF_CHAR
  CHARACTER(len=1), PARAMETER                            :: fillValueC = '-' 
  !! 'I', NF_SHORT
  INTEGER, PARAMETER                                     :: fillValueI = -9999
  !!! 'I16', NF_INT
  !INTEGER(kind=ii4), PARAMETER                           :: fillValueID = fI4
  !! 'R', NF_FLOAT
  REAL, PARAMETER                                        :: fillValueR = 1.e20
  !! 'D', NF_DOUBLE 
  REAL(kind=rr8), PARAMETER                              :: fillValueD = 1.d20
  
  CHARACTER(len=33)                                      :: errormsg = 'ERROR -- error -- ERROR -- error'

  ! error tolerance
  REAL(kind=rr8), PARAMETER                              :: tol=1.d-14
  ! Epsilon
  REAL, PARAMETER                                        :: reps=1.e-15
  REAL(kind=rr8), PARAMETER                              :: eps=1.d-15

  ! Maximum rank of a Frotran variable
  INTEGER, PARAMETER                                     :: maxfrank = 7
  
  ! simrad.AR version
  CHARACTER(len=S32)                                     :: simradARversion = '1.0'
  
  ! amount of times to debug values
  INTEGER, PARAMETER                                     :: Ndbg = 1

  ! frequency to debug values
  INTEGER, PARAMETER                                     :: Ndbgfrq = 100

  ! 4 / 3 radar beam height 
  REAL(kind=rr8), PARAMETER                              :: ke = 4.d0 / 3.d0
  
  ! Derived type with all the ranges for all water species for LUT
  TYPE lut_vars
    INTEGER                                              :: Nbands
    INTEGER                                              :: Nlutqxbins
    INTEGER                                              :: Nlutelevs
    CHARACTER(len=S16)                                   :: mp_scheme
  
    REAL, DIMENSION(:), ALLOCATABLE                      :: elevations

    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zh_RAIN
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zdr_RAIN
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: LDR_RAIN
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aih_RAIN
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aiv_RAIN
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: KDP_RAIN
    
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zh_SNOW
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zdr_SNOW
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: LDR_SNOW
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aih_SNOW
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aiv_SNOW
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: KDP_SNOW

    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zh_GRAU
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zdr_GRAU
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: LDR_GRAU
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aih_GRAU
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aiv_GRAU
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: KDP_GRAU

    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zh_ICE
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zdr_ICE
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: LDR_ICE
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aih_ICE
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aiv_ICE
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: KDP_ICE

    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zh_HAIL
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Zdr_HAIL
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: LDR_HAIL
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aih_HAIL
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: Aiv_HAIL
    REAL(kind=rr8), DIMENSION(:,:,:), ALLOCATABLE        :: KDP_HAIL    
    
    ! Which variables need to be interpolated / used depending on the microphysics scheme
    INTEGER                                              :: N_tocompute
    INTEGER, DIMENSION(12)                               :: tocompute
    !! Indices of the water species to compute / use
    INTEGER                                              :: Nmoist
    INTEGER, DIMENSION(12)                               :: imoist
    
    REAL                                                 :: q_rain_min
    REAL                                                 :: q_rain_max
    REAL                                                 :: qn_rain_min
    REAL                                                 :: qn_rain_max
    REAL                                                 :: q_snow_min
    REAL                                                 :: q_snow_max
    REAL                                                 :: qn_snow_min
    REAL                                                 :: qn_snow_max
    REAL                                                 :: q_grau_min
    REAL                                                 :: q_grau_max
    REAL                                                 :: qn_grau_min
    REAL                                                 :: qn_grau_max
    REAL                                                 :: q_ice_min
    REAL                                                 :: q_ice_max
    REAL                                                 :: qn_ice_min
    REAL                                                 :: qn_ice_max
    REAL                                                 :: q_hail_min
    REAL                                                 :: q_hail_max
    REAL                                                 :: qn_hail_min
    REAL                                                 :: qn_hail_max
    REAL                                                 :: q_clou_min
    REAL                                                 :: q_clou_max
    INTEGER                                              :: dim_ice
    INTEGER                                              :: dimv

    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_rain_vec
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_snow_vec
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_grau_vec
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_ice_vec
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_hail_vec
    
    ! Denisities
    REAL                                                 :: rhoo
    REAL                                                 :: rho_r
    REAL                                                 :: n0r
    REAL                                                 :: rho_s
    REAL                                                 :: n0s
    REAL                                                 :: rho_g
    REAL                                                 :: n0g
    REAL                                                 :: rho_h
    REAL                                                 :: n0h

  END TYPE lut_vars

  CONTAINS

  SUBROUTINE weights_9pts_passed1coords_ellipsoide_1D(i1, j1, id1x, ed1x, id1y, ed1y, a, f, elips,    &
    lon1, lat1, MaxNcoords, ij21Ncoords, Nc2, ij21coords, id2x, ed2x, id2y, ed2y, lon2, lat2, Na2,    &
    wghts, idbg)
  ! Subroutine to compute the weights of a inverse distance weighted interpolation using the
  !   3x3 grid points when already passing the coordinates of the target space at each reference point 
  !   on an ellipsoide
  !     lon1,lat1: space with data to interpolate (reference)
  !     lon2,lat2: space towards which interpolate (target)
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: id1x, ed1x, id1y, ed1y, MaxNcoords, Nc2,  &
      Na2
    REAL(kind=rr8), INTENT(in)                           :: a, f
    CHARACTER(len=S32), INTENT(in)                       :: elips
    INTEGER, INTENT(in)                                  :: id2x, ed2x, id2y, ed2y
    REAL, DIMENSION(id1x:ed1x,id1y:ed1y), INTENT(in)     :: lon1, lat1
    INTEGER, INTENT(in)                                  :: ij21Ncoords
    INTEGER, DIMENSION(MaxNcoords,Nc2), INTENT(in)       :: ij21coords
    REAL, DIMENSION(id2x:ed2x,id2y:ed2y), INTENT(in)     :: lon2, lat2
    LOGICAL, OPTIONAL                                    :: idbg
    REAL, DIMENSION(MaxNcoords,Na2), INTENT(out)         :: wghts
    
    ! Local
    INTEGER                                              :: i1, j1, ic, i2, j2, ij, icc
    INTEGER                                              :: iwi, iwe, jwi, jwe
    REAL                                                 :: sumdists, totwght
    INTEGER, DIMENSION(Nc2)                              :: imaxwght, izero
    REAL, DIMENSION(3,3)                                 :: ijlon1, ijlat1, dists, wgs
    REAL, DIMENSION(Na2)                                 :: dist9
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! id1x, ed1x, id1y, ed1y: Extent of the reference space
    ! lon1, lat1: coordinates of the reference space
    ! MaxNcoords: maximum amount of target space grid points in a single reference grid point
    ! ij21Ncoords: amount of target grid points in a single reference grid point
    ! Nc2: amount of coordinates in a 2-dimensional space
    ! ij21coords: coordinates of the target grid points in a single reference grid point
    ! id2x, ed2x, id2y, ed2y: Extent of the target space
    ! lon2, lat2: coordinates of the target space
    ! Na2: amount of coordinates around in a 2-dimensional space
    ! wghts: weights to use to interpolate values of a target grid point to its closest 3x3 grid 
    !   points in the reference space
    ! ij21dists: distances of a target grid point to its closest 3x3 grid points in the reference space
    fname = 'weights_9pts_passed1coords_ellipsoide_1D'
    
    wghts = 0.
    
    dbg = set_dbg(idbg)
    
    icc = 1
    IF (ij21Ncoords == 0) RETURN
    DO ic=1, ij21Ncoords
      iwi = -1
      iwe = -1
      jwi = -1
      jwe = -1
    
      IF (i1 == id1x) iwi = 0
      IF (i1 == ed1x) iwe = 0
      IF (j1 == id1y) jwi = 0
      IF (j1 == ed1y) jwe = 0
      ijlon1 = 0.
      ijlat1 = 0.
      ! 3x3 array as 
      !   -1,-1,(1,1): SW; 0,-1,(2,1): SC, 1,-1,(3,1): SE; 
      !    -1,0,(1,2): CW;  0,0,(2,2): CC;  1,0,(3,2): CE;
      !    -1,1,(1,3): NW;  0,1,(2,3): NC;  1,1,(3,3): NE;
      DO i2=iwi, iwe
        DO j2=jwi, jwe
          ijlon1(i2+2,j2+2) = lon1(i1+i2,j1+j2)
          ijlat1(i2+2,j2+2) = lat1(i1+i2,j1+j2)
        END DO
      END DO
      CALL real_distances_ellipsoide(1, 3, 1, 3, a, f, elips, ijlon1, ijlat1,                         &
        lon2(ij21coords(ic,1),ij21coords(ic,2)), lat2(ij21coords(ic,1),ij21coords(ic,2)), dists)
       ! Filtering for dist == 0!!
      IF (ANY(dists == 0.)) THEN
        izero = MINLOC(dists)
        wgs = 0.
        wgs(izero(1), izero(2)) = 1.
        sumdists = 1.
      ELSE
        sumdists = SUM(1./dists)
        wgs = 0.
        wgs = 1./dists/sumdists
        totwght = SUM(wgs)

        ! Fixing
        IF (totwght < 1.) THEN
          imaxwght = MAXLOC(wgs)
          wgs(imaxwght(1),imaxwght(2)) = wgs(imaxwght(1),imaxwght(2)) + 1. - totwght
        END IF
      END IF
          
      ij = 1
      dist9 = 0.          
      wghts(ic,:) = 0.
      DO i2=1,3
        DO j2=1,3
          dist9(ij) = dists(i2,j2)
          wghts(ic,ij) = wgs(i2,j2)
          ij = ij + 1
        END DO
      END DO
          
      IF (dbg .AND. icc < 10) PRINT *, '  ' // TRIM(fname) // ": i,j,ic ", i1,j1, ic, ':',            &
        lon1(i1,j1), lat1(i1,j1), ' 2: ', lon2(ij21coords(ic,1),ij21coords(ic,2)),                    &
        lat2(ij21coords(ic,1),ij21coords(ic,2)), ' dists: ', dist9, ' wghts: ',                       &
        wghts(ic,:), ' totw ', SUM(wghts(ic,:))
      icc = icc + 1
                      
    END DO
    
    RETURN
    
  END SUBROUTINE weights_9pts_passed1coords_ellipsoide_1D

  CHARACTER(len=S32) FUNCTION cvals_varn(icval)
  ! Function to provide the name of the variable accordingly its index from cvals array

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: icval
    
    ! Local
    INTEGER                                              :: icvalmax
    CHARACTER (LEN=S2048)                                :: diag_message
    CHARACTER(len=S32)                                   :: fname
    !!!!!!! Variables
    ! icval: index of the variable in 'cvals' (see phys/module_simradar.F) 

    fname = 'cvals_varn'

    icvalmax = 12

    cvals_varn = ''

    IF (icval == 1) THEN
      cvals_varn = 'QVAPOR'
    ELSE IF (icval == 2) THEN
      cvals_varn = 'QCLOUD'
    ELSE IF (icval == 3) THEN
      cvals_varn = 'QRAIN'
    ELSE IF (icval == 4) THEN
      cvals_varn = 'QICE'
    ELSE IF (icval == 5) THEN
      cvals_varn = 'QSNOW'
    ELSE IF (icval == 6) THEN
      cvals_varn = 'QGRAUPEL'
    ELSE IF (icval == 7) THEN
      cvals_varn = 'QHAIL'
    ELSE IF (icval == 8) THEN
      cvals_varn = 'U'
    ELSE IF (icval == 9) THEN
      cvals_varn = 'V'
    ELSE IF (icval == 10) THEN
      cvals_varn = 'W'
    ELSE IF (icval == 11) THEN
      cvals_varn = 'T_PHY'
    ELSE IF (icval == 12) THEN
      cvals_varn = 'PRESS (P+PB)'
    ELSE
      WRITE(diag_message, *) &
      TRIM(errormsg) // ":   " //TRIM(fname)// ": icval: " // TRIM(ItoS(icval)) // " not ready !!"//  &
        " It has to be between 1 and " // TRIM(ItoS(icvalmax))
      CALL wrf_error_fatal ( diag_message )
    END IF

    RETURN    

  END FUNCTION cvals_varn
  
  SUBROUTINE create_logarithmic_scaled_vectorsD(vec_min, vec_max, Nvals, logvec, idbg)
  ! Function to create logarithmic-based equi-spaced scaled double Real vectors
  
    IMPLICIT NONE 
    
    INTEGER, INTENT(in)                                  :: Nvals
    REAL, INTENT(in)                                     :: vec_min, vec_max
    LOGICAL, OPTIONAL                                    :: idbg
    REAL(kind=rr8), DIMENSION(NVals), INTENT(out)        :: logvec

    ! Local
    INTEGER                                              :: i
    REAL(kind=rr8)                                       :: log10n, log10x, ddvec
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! vec_min: minimum value
    ! vec_max: maximum value
    ! Nvals: amount of values
    fname = 'create_logarithmic_scaled_vectorsD'

    dbg = set_dbg(idbg)
    
    log10n = DLOG10(vec_min*1.D0)
    log10x = DLOG10(vec_max*1.D0)
    
    ddvec = (log10x - log10n)/(Nvals-1)*1.D0
    
    DO i=1, Nvals
      logvec(i) = (10.D0)**(log10n + (i-1)*ddvec)
    END DO

    IF (dbg) THEN
      PRINT *, '    ' // TRIM(fname) // ": log10 min: ", log10n, ' max: ', log10x,                    &
        ' estimated distance: ', ddvec       
      PRINT *, '      LOG10 distance between consecutive log10(1), log10(2), dist:',                  &
        DLOG10(logvec(1)), ', ', DLOG10(logvec(2)), ', ', DLOG10(logvec(2)) - DLOG10(logvec(1))
      PRINT *, '      LOG10 distance between consecutive log10(Ndim-1), log10(Ndim), dist:',          &
        DLOG10(logvec(Nvals-1)), ', ',  DLOG10(logvec(Nvals)), ', ', DLOG10(logvec(Nvals)) -          &
        DLOG10(logvec(Nvals-1))
    END IF

    RETURN

  END SUBROUTINE create_logarithmic_scaled_vectorsD
  
  REAL(kind=rr8) FUNCTION lin_interpolate_1RR(pt, Ndim, coords, vals)
! If willing to manage error as function
!  REAL(kind=rr8) FUNCTION lin_interpolate_1RD(pt, Ndim, coords, vals, statv)
  ! Function to provide linear interpolation a Real location on a given real axis

    IMPLICIT NONE
   
    INTEGER, INTENT(in)                                  :: Ndim
    REAL, INTENT(in)                                     :: pt
    REAL, DIMENSION(Ndim), INTENT(in)                    :: coords, vals
    ! 0: no error; 1: error
!    INTEGER, INTENT(out)                                 :: statv
   
    ! Local
    INTEGER                                              :: icor
    INTEGER                                              :: refi
    REAL                                                 :: minc, maxc
    REAL                                                 :: signv, a, b
    REAL, DIMENSION(Ndim)                                :: origsing
    CHARACTER (LEN=S2048)                                :: diag_message
    CHARACTER(len=Sl)                                    :: fname

    !!!!!!! Variables
    ! pt: pt to interpolate
    ! coords: coordinates along the axis (incremental)
    ! vals: values at the coordinates
    fname = 'lin_interpolate_1RR'
  
    lin_interpolate_1RR = fillValueD  
    ! Sign of the axis
    origsing = ABS(coords(Ndim)-coords(1))/(coords(Ndim)-coords(1))
  
    minc = MINVAL(coords)
    maxc = MAXVAL(coords)

    a = fillValueR
    refi = -1
    IF (pt < minc) THEN
      a = vals(1)
      refi = 1
      signv = -1.
    ELSE IF (pt > maxc) THEN
      a = vals(Ndim)
      refi = Ndim-1
      signv = 1.
    ELSE IF (pt == minc) THEN
      lin_interpolate_1RR = vals(1)
      refi = 1
      signv = 1.
      
      RETURN 
    ELSE IF (pt == maxc) THEN
      lin_interpolate_1RR = vals(Ndim)
      refi = Ndim-1
      signv = 1.

      RETURN     
    ELSE
      DO icor=1, Ndim-1
        IF (pt == coords(icor)) THEN
          lin_interpolate_1RR = vals(icor)
          RETURN
        ELSE IF (pt >= coords(icor) .AND. pt < coords(icor+1)) THEN
          a = vals(icor)
          refi = icor
          signv = 1.
          EXIT
        END IF
      END DO
    END IF

!    statv = 0
    IF (refi == -1) THEN
      WRITE (diag_message , * ) &
      TRIM(errormsg) // ":   " //TRIM(fname)// ": No interpolation found for: " // TRIM(RtoS(pt)) //  &
        ' within range: [' // TRIM(RtoS(minc)) // ', ' // TRIM(RtoS(maxc))
      CALL wrf_error_fatal ( diag_message )
    END IF

    b = ABS(pt-coords(refi))/(coords(refi+1)-coords(refi))
    
    lin_interpolate_1RR = a + signv*(vals(refi+1)-a)*b

    RETURN

  END FUNCTION lin_interpolate_1RR

  SUBROUTINE com_gamma(x, ga)
  ! GAMMA FUNCTION
  !-----------------------------------------------------------------------
  
  !==================================================
  !       Purpose: Compute the gamma function â(x)
  !       Input :  x  --- Argument of â(x)
  !                       ( x is not equal to 0,-1,-2,úúú )
  !       Output:  GA --- â(x)
  !==================================================
  ! Proff. Jianming Jin
  ! Department of Electrical and Computer Engineering
  ! University of Illinois
  
    IMPLICIT NONE

    REAL(kind=rr8), INTENT(in)                           :: X
    REAL(kind=rr8), INTENT(out)                          :: GA

    ! Local
    REAL(kind=rr8), DIMENSION(26)                        :: G
    REAL(kind=rr8)                                       :: Z , R , GR
    INTEGER                                              :: M1, K , M
    CHARACTER(len=Sl)                                    :: fname
    
    DATA G /1.0D0, 0.5772156649015329D0, -0.6558780715202538D0, -0.420026350340952D-1,                &
      0.1665386113822915D0,-.421977345555443D-1, -.96219715278770D-2, .72189432466630D-2,             &
      -.11651675918591D-2, -.2152416741149D-3, .1280502823882D-3, -.201348547807D-4,                  &
      -.12504934821D-5, .11330272320D-5, -.2056338417D-6, .61160950D-8, .50020075D-8, -.11812746D-8,  &
      .1043427D-9, .77823D-11, -.36968D-11, .51D-12, -.206D-13, -.54D-14, .14D-14, .1D-15 /

    !!!!!!! Variables
    ! x: argument 
    fname = 'com_gamma'
    
    IF (X.EQ.INT(X)) THEN
      IF (X.GT.0.0D0) THEN
        GA=1.0D0
        M1=X-1
        DO K=2,M1
          GA=GA*K
        END DO
      ELSE
        GA=1.0D+300
      ENDIF

    ELSE
      IF (DABS(X).GT.1.0D0) THEN
        Z = DABS(X)
        M = INT(Z)
        R = 1.0D0
        DO K=1, M
          R = R*(Z-K)
        END DO
        Z=Z-M
      ELSE
        Z=X
      ENDIF
      GR=G(26)
      DO K = 25, 1, -1
        GR = GR*Z + G(K)
      END DO

      GA=1.0D0/(GR*Z)
      IF (DABS(X).GT.1.0D0) THEN
        GA=GA*R
        IF (X.LT.0.0D0) GA=-DPICT/(X*GA*DSIN(DPICT*X))
      ENDIF
    ENDIF

    RETURN

  END SUBROUTINE com_gamma

  SUBROUTINE Polarimetrics_fromLUT_1R(qx, iz, ibnd, dimv, q_log_vec, Nbins, Zh_x, Zdr_x, KDP_x,       &
    Zh_out, Zdr_out, KDP_out, idbg)
  ! Subroutine to interpolate the polarimetic characteristics of a given water specie following T-matrix
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: iz, ibnd, dimv, Nbins
    REAL, INTENT(in)                                     :: qx
    REAL, DIMENSION(dimv), INTENT(in)                    :: q_log_vec
    REAL, DIMENSION(Nbins), INTENT(in)                   :: Zh_x, Zdr_x, KDP_x
    LOGICAL, OPTIONAL                                    :: idbg
    REAL, INTENT(out)                                    :: Zh_out, Zdr_out, KDP_out
    
    ! Local
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! qx: value of the water specie
    ! iz: lut azimuth index
    ! ibnd: band
    !   1: 'C'
    !   2: 'S'
    !   3: 'X'
    ! dimv: length of the logarithmic scale
    ! q_log_vec: logarithmic scale
    ! Nbins: amount of bins 
    ! Zh_x: values of 'Zh' for the specific water specie
    ! Zdr_x: values of 'Zdr' for the specific water specie
    ! KDP_x: values of 'KDP' for the specific water specie
    fname = 'Polarimetrics_fromLUT_1R'
 
    dbg = set_dbg(idbg)
    
    IF (dbg) THEN
      PRINT *, '  ' // TRIM(fname) // ": q_log_vec 1 a 5:", q_log_vec(1:5)
      PRINT *, '  ' // TRIM(fname) // ": Zh_x at ibnd= ", ibnd ," iz= ", iz, " 1 a 5:", Zh_x(1:5)
    END IF
    
    dbg = set_dbg(idbg)

    Zh_out = fillValueR
    Zdr_out = fillValueR
    KDP_out = fillValueR
    
    Zh_out = lin_interpolate_1RR(qx, dimv, q_log_vec, Zh_x)
    Zdr_out = lin_interpolate_1RR(qx, dimv, q_log_vec, Zdr_x)
    KDP_out = lin_interpolate_1RR(qx, dimv, q_log_vec, KDP_x)
                  
    RETURN
  
  END SUBROUTINE Polarimetrics_fromLUT_1R
 
  REAL(kind=rr8) FUNCTION calc_vt(a, b, lambdax)
  ! Function to compute terminal velocity after a gamma function
  
    IMPLICIT NONE
    
    REAL, INTENT(in)                                     :: a, b
    REAL(kind=rr8), INTENT(in)                           :: lambdax
    
    ! Local
    REAL(kind=rr8)                                       :: db, gammav, lambda6
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! a, b: gamma parameters
    ! lambdax: lambda of the water specie
    
    fname = 'calc_vt' 
    
    calc_vt = fillValueD
    
    db = b*1.D0
    CALL com_gamma(4.D0+db, gammav)
    
    calc_vt = a*gammav/6
    lambda6 = lambdax**db
    
    calc_vt = calc_vt / lambda6

    RETURN
    
  END FUNCTION calc_vt
  
  SUBROUTINE calc_tvelocities_1D(qr, qs, qg, qh, p, T, lutgrid, v_rain, v_snow, v_grau, v_hail)
  ! Subroutine to compute terminal velocities for water species 1D version
  !
  ! Calculate the MASS-WEIGHTED MEAN terminal velocities following Srivastava (1967).
  !    int vt l(D) dD/l, where vt is the terminal velocity of diamater D, l(D) is the mixing
  !    ratio of the precipitation particle of diameter D, and l is the mixing ratio of a precip. field
  !    see eq. 11-13 in Lin et al. 1983. (same as wrf2radar)

    IMPLICIT NONE
    
    REAL, INTENT(in)                                     :: qr, qs, qg, qh, p, T
    TYPE(lut_vars), INTENT(in)                           :: lutgrid
    REAL(kind=rr8), INTENT(out)                          :: v_rain, v_snow, v_grau, v_hail
    
    ! Local
    REAL                                                 :: rhoo, rho_air, rhofac
    REAL                                                 :: rho_r, n0r, rho_s, n0s, rho_g, n0g,       &
      rho_h, n0h
    REAL(kind=rr8)                                       :: lambdar, lambdas, lambdag, lambdah
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! qr: mixing ratio for rain [kgkg-1]
    ! qs: mixing ratio for snow [kgkg-1]
    ! qg: mixing ratio for graupel [kgkg-1]
    ! qh: mixing ratio for hail [kgkg-1]
    ! p: pressure [Pa]
    ! T: air temperature [K]
    ! lutgrid: derived type with all the T-MATRIX look-up-table values
    fname = 'calc_tvelocities_1D'
    
    !! Attention bring herer rl rhofactor! 
    !!!! THIS must be generalized from lutgrid !!!!
    !!! reference air density used in some WSM6 docs, kg m^-3
    rhoo = lutgrid%rhoo
    !!! kg/m3
    rho_air  = p / (287.058 * T)
    rhofac  = (rhoo/rho_air)**0.5
    
    rho_r = lutgrid%rho_r
    n0r   = lutgrid%n0r
    rho_s = lutgrid%rho_s
    n0s   = lutgrid%n0s
    rho_g = lutgrid%rho_g
    n0g   = lutgrid%n0g
    rho_h = lutgrid%rho_h
    n0h   = lutgrid%n0g
    
    v_rain = fillValueD
    v_snow = fillValueD
    v_grau = fillValueD

    lambdar = fillValueR
    lambdas = fillValueR
    lambdag = fillValueR

    ! Rain
    IF (qr > 0.) THEN
      lambdar = (dpict*rho_r*n0r / (qr))**(1./4.)
      v_rain = calc_vt(841.9, 0.8, lambdar)*rhofac
    END IF
    
    ! Snow
    IF (qs > 0.) THEN
      lambdas = (dpict*rho_s*n0s / (qs))**(1./4.)
      v_snow = calc_vt(11.72, 0.41, lambdas)*rhofac
    END IF

    ! Grau
    IF (qg > 0.) THEN
      lambdag = (dpict*rho_g*n0g / (qg))**(1./4.)
      v_grau = calc_vt(330., 0.8, lambdag)*rhofac
    END IF

    ! Hail (L. Fita not known range !!)
    IF (qh > 0.) THEN
      lambdah = (dpict*rho_h*n0h / (qh))**(1./4.)
      v_hail = calc_vt(-9999., -9999., lambdah)*rhofac
    END IF

    RETURN

  END SUBROUTINE calc_tvelocities_1D
  
  SUBROUTINE weights_hcoordsR27around_polar_ellipsoide_1D(i1, j1, id1x, ed1x, id1y, ed1y, id1z,       &
    ed1z, a, f, elips, lon1, lat1, vert1, miss_v1, maxNcoords, ij21Ncoords, Nc2, ij21coords,          &
    Na2, id2x, ed2x, id2y, ed2y, id2z, ed2z, lon2, lat2, vert2, miss_v2, maxN3coords, ijk21Ncoords,   &
    Nc3, ijk21coords, totngridradar, Na3, wghts, idbg)
  ! Subroutine to compute weights of the 27 around grid points to 3-dimensional interpolate to a new 
  !   space using the list of cells on the horizontal of the targeted space that lay within each grid 
  !   cell of the reference space with the values to interpolate version on an ellipsoide
  !    lon1,lat1,vert1: space with data to interpolate (reference)
  !    lon2,lat2,vert2: polar space towards which interpolate (target)
  !    id1x= i1-1, edx1= i1+1; id1y= j1-1, edy1= j1+1; id1z=1, ed1z=dimz
  !   NOTE in polar space:
  !     dx: distance    dy: angle    dz: azimuth
  !     3Dvariables (dy,dz,dx),        lon,lats(dy,dx),        hgt(dx,dz)
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: i1, j1, id1x, ed1x, id1y, ed1y, id1z,     &
      ed1z, MaxNcoords, MaxN3coords, Nc2, Na2, Nc3, Na3
    INTEGER, INTENT(in)                                  :: id2x, ed2x, id2y, ed2y, id2z, ed2z
    REAL, INTENT(in)                                     :: miss_v1, miss_v2
    REAL, INTENT(in)                                     :: a, f
    CHARACTER(len=30), INTENT(in)                        :: elips
    REAL, DIMENSION(id1x:ed1x,id1y:ed1y), INTENT(in)     :: lon1, lat1
    REAL, DIMENSION(id1x:ed1x,id1y:ed1y,id1z:ed1z),                                                   &
      INTENT(in)                                         :: vert1
    INTEGER, INTENT(in)                                  :: ij21Ncoords
    INTEGER, DIMENSION(MaxNcoords,Nc2), INTENT(in)       :: ij21coords
    REAL, DIMENSION(id2y:ed2y,id2x:ed2x), INTENT(in)     :: lon2, lat2
    REAL, DIMENSION(id2z:ed2z,id2x:ed2x), INTENT(in)     :: vert2
    LOGICAL, OPTIONAL                                    :: idbg
    INTEGER, INTENT(out)                                 :: totngridradar
    INTEGER, DIMENSION(id1z:ed1z), INTENT(out)           :: ijk21Ncoords
    REAL, DIMENSION(id1z:ed1z,maxN3coords,Na3),                                                       &
      INTENT(out)                                        :: wghts
    INTEGER, DIMENSION(id1z:ed1z,maxN3coords,Nc3),                                                    &
      INTENT(out)                                        :: ijk21coords
    
    ! Local
    INTEGER                                              :: k1, ic, i2, j2, k2, kk2, ij, ijk
    INTEGER                                              :: ici, icj, ick
    INTEGER                                              :: ii1i, ii1e, ii2i, ii2e, ii3i, ii3e
    INTEGER                                              :: ijkfound
    REAL                                                 :: ijkvert1, sumdists, totwght
    REAL(kind=rr8)                                       :: da, df
    INTEGER, DIMENSION(1)                                :: kmindist
    INTEGER, DIMENSION(1)                                :: imaxwght, izero
    REAL, DIMENSION(Na2)                                 :: ijlon1, ijlat1, ijdists
    REAL, DIMENSION(Na3)                                 :: zdists, dists, wgs
    REAL, DIMENSION(id1z:ed1z)                           :: vert1mask, zdist1
    LOGICAL, DIMENSION(Na3)                              :: ijkmask
    LOGICAL, DIMENSION(id1z:ed1z)                        :: firstijk21Ncoords
    LOGICAL                                              :: printdbg
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! id1x, ed1x, id1y, ed1y, id1z, ed1z: Extent of the reference space
    ! a: semimajor axis equatorial (in meters)
    ! f: flattening
    ! elips: name of the ellipsis being used
    ! lon1, lat1, vert1: coordinates of the reference space [degrees_east, degrees_north, m]
    ! MaxNcoords: maximum amount of target space grid points in a single reference grid point
    ! ij21Ncoords: amount of target grid points in a single reference grid pointv
    ! Nc2: amount of coordinates of 2-dimensional space
    ! ij21coords: coordinates of the target grid points in a single reference grid point
    ! ij21dists: horizonatal distancess of the target grid points in a single reference grid point
    ! id2x, ed2x, id2y, ed2y, id2z, ed2z: Extent of the target space
    ! lon2, lat2, vert2: coordinates of the target space [degrees_east, degrees_north, m]
    ! ijk21Ncoords: amount of target grid 3-dimensional points in a single reference grid pointv 
    ! Nc3: amount of coordinates of 3-dimensional space
    ! ijk21coords: 3-dimensional coordinates of the target grid points in a single reference grid point
    ! Na3: amount of points around in a 3-dimensional space
    ! wghts: weights to use to interpolate values of a target grid point to its closest 3x3x3 grid 
    !   points in the reference space
    ! Rdists: distances of a target grid point to its closest 3x3x3 grid points in the reference space
    fname = 'weights_hcoordsR27around_polar_ellipsoide_1D'
     
    dbg = set_dbg(idbg)
        
    ijk21Ncoords = 0
    ijk21coords = fillValueI
    wghts = 0.
    firstijk21Ncoords = .TRUE.
    
    da = a*1.D0
    df = f*1.D0

    ! Capping run of 9 weights if they lay in the border
    ii1i = -1
    ii1e = 1
    ii2i = -1
    ii2e = 1
    IF (i1 == id1x) ii1i = 0
    IF (i1 == ed1x) ii1e = 0
    IF (j1 == id1y) ii2i = 0
    IF (j1 == ed1y) ii2e = 0

    IF (dbg) PRINT *,'  ' // TRIM(fname) // ": Compute 3-dimensional weights at ", i1, ', ', j1 ,     &
      ' with ', ij21Ncoords, " radar cells horizontally"
    
    totngridradar = 0
    ijkfound = 0
    
    gridpts: DO ic=1, ij21Ncoords
      ! Indices of target space within reference grid 
      !! angle
      ici = ij21coords(ic,1)
      !! dist
      icj = ij21coords(ic,2)

      ! computing distances from 3x3 horizontal grid (the already computed is only for i1,j1)
      ! horizontal 3x3 array as 
      !   -1,-1,(1,1)[1]: SW; 0,-1,(2,1)[2]: SC, 1,-1,(3,1)[3]: SE; 
      !    -1,0,(1,2)[4]: CW;  0,0,(2,2)[5]: CC;  1,0,(3,2)[6]: CE;
      !    -1,1,(1,3)[7]: NW;  0,1,(2,3)[8]: NC;  1,1,(3,3)[9]: NE;
      ij = 1
      ijlon1 = 0.
      ijlat1 = 0.
      DO i2=ii1i,ii1e
        DO j2=ii2i,ii2e
          ijlon1(ij) = lon1(i1+i2,j1+j2)
          ijlat1(ij) = lat1(i1+i2,j1+j2)
          ij = ij + 1
        END DO
      END DO

      
      CALL real_distances(1, 3, 1, 3, ijlon1, ijlat1, lon2(ici,icj), lat2(ici,icj), ijdists)

      ! Computing vertical distances (from closest i,j)
      zloop: DO k2=id2z, ed2z
        IF (vert2(k2,icj) /= miss_v2) THEN
          
          printdbg = .FALSE.         
          IF (dbg) THEN
            IF (MOD(ic+k2, Ndbgfrq) == 0 .AND. ijkfound < Ndbg) THEN
              printdbg = .TRUE.
              ijkfound = ijkfound + 1
            END IF
          END IF

          ! vertical level at center of 3x3
          DO k1=id1z, ed1z
            IF (vert1(i1,j1,k1) == miss_v1) THEN
              vert1mask(k1) = fillValueR
            ELSE
              vert1mask(k1) = vert1(i1,j1,k1)
            END IF
          END DO

          zdist1 = ABS( vert1mask - vert2(k2,icj) )
          kmindist = MINLOC(zdist1)
          ick = kmindist(1)
              
          IF (printdbg) PRINT *,'    ' // TRIM(fname) // ": i1,j1:", i1, j1, " i2, j2, k2: ",         &
            ici, icj, k2, ' vert1 masked:', vert1mask, ' vert2: ', vert2(k2,icj), ' zdist1: ',        &
            zdist1, ' ick:', ick, ' found: ', totngridradar+1

          !! This is not properly workking, ijk21Ncoords = 0 at the begining is not woring ¿?
          IF (firstijk21Ncoords(ick)) THEN
            ijk21Ncoords(ick) = 1
            firstijk21Ncoords(ick) = .FALSE.
          ELSE 
            ijk21Ncoords(ick) = ijk21Ncoords(ick) + 1
          END IF
          IF (ijk21Ncoords(ick) > maxN3coords) THEN
            CALL ErrMsg("Too much radar cells: " // TRIM(ItoS(ijk21Ncoords(ick))) //                  &
              " in a WRF grid. Increase 'ngrid_radar3' in namelist current value: " //                &
              TRIM(ItoS(maxN3coords)), fname, -1)
          END IF 
          ijk21coords(ick,ijk21Ncoords(ick),1) = ici
          ijk21coords(ick,ijk21Ncoords(ick),2) = icj
          ijk21coords(ick,ijk21Ncoords(ick),3) = k2
          totngridradar = totngridradar + 1
          
          ! vertical 27 array as 
          !   SW
          !   -1,-1,-1,(1,1,1),[1]: bottom; -1,-1,0,(1,1,2),[2]: middle, -1,-1,1,(1,1,3),[3]: top; 
          !   (...)
          !   NE
          !   1,1,-1,(3,3,1),[25]: bottom; 1,1,0,(3,3,2),[26]: middle, 1,1,1,(3,3,3),[27]: top; 
          ij = 1
          ijk = 1
          dists = fillValueR
          ijkmask = .FALSE.
          ii3i = -1
          ii3e = -1
          IF (ick == 1) ii3i = 0
          IF (ick == ed1z) ii3e = 0
          DO i2=ii1i,ii1e
            DO j2=ii2i,ii2e
              DO kk2=ii3i,ii3e
                ijkvert1 = vert1(i1+i2,j1+j2,ick+kk2)
                IF (ijkvert1 /= miss_v1) THEN
                  ! 3-dimensional distance
                  !   dist = SQRT(hdist**2 + zdist**2)
                  dists(ijk) = SQRT(ijdists(ij)**2 + (ijkvert1-vert2(k2,icj))**2)
                  ijkmask(ijk) = .TRUE.
                END IF
                ijk = ijk + 1
              END DO
              ij = ij + 1
            END DO
          END DO              
              
          ! Filtering for dist == 0!!
          wgs = 0.
          IF (ANY(dists == 0.)) THEN
            izero = MINLOC(dists)
            wgs = 0.
            wgs(izero(1)) = 1.
            sumdists = 1.
          ELSE
            sumdists = SUM(1./dists, mask=ijkmask)
            ! Applying mask in weights
            DO k1=1, Na3
              IF (.NOT.ijkmask(k1)) THEN
                wgs(k1) = 0.
              ELSE                 
                wgs(k1) = 1./dists(k1)/sumdists
              END IF
            END DO
            totwght = SUM(wgs)
            ! Fixing
            IF (totwght < 1.) THEN
              imaxwght = MAXLOC(wgs)
              wgs(imaxwght(1)) = wgs(imaxwght(1)) + 1. - totwght
            END IF
          END IF
          wghts(ick,ijk21Ncoords(ick),:) = wgs

          IF (printdbg) THEN
            PRINT *, '    '//TRIM(fname) // ": i,j,k radar cell: ",                                   &
              ijk21coords(ick,ijk21Ncoords(ick),:), ' hgt: ', vert2(k2,icj), ' i,j,k reference: ',    &
              i1,j1,ick, ' dists: ', dists, ' mask: ', ijkmask, ' wghts: ',                           &
              wghts(ick,ijk21Ncoords(ick),:)
          END IF
        END IF
      END DO zloop

    END DO gridpts
    
    IF (totngridradar > ij21Ncoords*(ed2z-id2z+1) ) THEN
      CALL ErrMsg('Too many radar cells were found! Theoretical maximum: ' //                         &
        TRIM(ItoS(ij21Ncoords*(ed2z-id2z+1))) // ', but I found: ' // TRIM(ItoS(totngridradar)),      &
        fname, -1)
    END IF
    
    
    IF (dbg) THEN
      PRINT *,'    ' // TRIM(fname) // ": amount of vertical radar cells in WRF grid: ", i1, j1,      &
        ' = ', totngridradar
      PRINT *, '    found along the reference column [amount (k_lev)]:',                              &
        (TRIM(ItoS(ijk21Ncoords(k1))) // ' (' // TRIM(ItoS(k1)) // '), ', k1=id1z, ed1z)
    END IF
            
    RETURN
  
  END SUBROUTINE weights_hcoordsR27around_polar_ellipsoide_1D

  SUBROUTINE compute_radar_lonlat_radial(Nangle, Nazimuth, Ndist, lon_pt, lat_pt, hgt_pt, ddist,      &
    max_dist, max_hgt, azimuths, a, f, elips, angle1d, dist1d, azimuth1d, lons, lats, hgts, lhgts,    &
    lon_bnds, lat_bnds, hgt_bnds, idbg)
  ! Subroutine to copmute lon, lat, hgt coordinates in the radar matrix following polar coordinates
    
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: Nangle, Nazimuth, Ndist
    REAL, INTENT(in)                                     :: lon_pt, lat_pt, hgt_pt, ddist, max_dist,  &
      max_hgt
    REAL, DIMENSION(Nazimuth), INTENT(in)                :: azimuths
    REAL(kind=rr8), INTENT(in)                           :: a, f
    CHARACTER(len=S32), INTENT(in)                       :: elips
    LOGICAL, OPTIONAL                                    :: idbg
    REAL, DIMENSION(Nangle), INTENT(out)                 :: angle1d
    REAL, DIMENSION(Ndist), INTENT(out)                  :: dist1d
    REAL, DIMENSION(Nazimuth), INTENT(out)               :: azimuth1d
    REAL, DIMENSION(Nangle,Ndist), INTENT(out)           :: lons, lats
    REAL, DIMENSION(Nazimuth,Ndist), INTENT(out)         :: hgts
    REAL, DIMENSION(Nazimuth,Ndist), INTENT(out)         :: lhgts
    REAL, DIMENSION(Nangle,Ndist,4),INTENT(out)          :: lon_bnds, lat_bnds
    REAL, DIMENSION(Nazimuth,Ndist,4),INTENT(out)        :: hgt_bnds
        
    ! Local
    INTEGER                                              :: i, j, k, ijfound
    REAL                                                 :: dangle, ddist2
    REAL                                                 :: angle, dist, azimuth, hgt, lhgt
    REAL                                                 :: lonv, latv
    REAL, DIMENSION(Nazimuth)                            :: radazimuths, dhgt, dhgt2
    LOGICAL                                              :: printdbg
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname
    
    !!! Test inverse distance
    REAL                                                 :: faz, s
    
    !!!!!!! Variables
    ! Nangle, Nazimuth, Ndist: dimension of the radar space
    ! lon_pt, lat_pt, hgt_pt: location of the radar
    ! ddist: horizontal resolution [m]
    ! max_dist: maximum distance of the radar [km]
    ! max_hgt: maximum height of the radar [km]
    ! azimuths: angles of the azimuths [degree]
    ! a, f: axis / semi-axis of the ellipsoide
    ! offset: offset to apply to the radar array
    ! iis, jjs, kks: indices of the radar array
    ! hdists: horizontal distances of the radar 
    ! zdists: vertical distances of the radar 
    ! lons, lats, hgts: coordinates of the radar array
    ! lon_bnds, lat_bnds, hgt_bnds: vertices of the grids of the radar array
    fname = 'compute_radar_lonlat_radial'
    
    dbg = set_dbg(idbg)
    
    ! Initializing
    angle1d = fillValueR
    dist1d = fillValueR
    azimuth1d = azimuths
    lhgts = fillValueR
    lons = fillValueR
    lats = fillValueR
    hgts = fillValueR
    lat_bnds = fillValueR
    lon_bnds = fillValueR
    hgt_bnds = fillValueR
    
    ijfound = 0
    ! Angle at the center of the ray
    dangle = 360. / (Nangle)
    ! Angle at the center of the ray
    radazimuths = azimuths / DegRad
    dhgt = 0.
    dhgt(1:Nazimuth-1) = radazimuths(2:Nazimuth) - radazimuths(1:Nazimuth-1)
    dhgt(Nazimuth) = dhgt(Nazimuth-1)
    
    ddist2 = ddist / 2.
    dhgt2 = dhgt / 2.
    
    angl: DO i=1, Nangle
      angle = dangle*(i-1)
      angle1d(i) = angle
            
      DO j=1, Ndist
        printdbg = .FALSE.
        IF (dbg) THEN
          IF (MOD(i+j,Ndbgfrq) == 0) THEN
            ijfound = ijfound + 1
            IF (ijfound < Ndbg) printdbg = .TRUE.
          END IF
        END IF

        ! Centering at half distance
        dist = ddist*(j-1) + ddist2
        IF (i == 1) dist1d(j) = dist

        ! Vertices
        ! Cell boundaries
        ! SW
        dist = ddist*(j-1.)
        angle = dangle*(i-1.-0.5)
        CALL position_from_point(lon_pt, lat_pt, a, f, elips, dist, angle, lonv, latv)
        lon_bnds(i,j,1) = lonv
        lat_bnds(i,j,1) = latv
          
        ! SE
        dist = ddist*(j+1.)
        angle = dangle*(i-1.-0.5)
        CALL position_from_point(lon_pt, lat_pt, a, f, elips, dist, angle, lonv, latv)
        lon_bnds(i,j,2) = lonv
        lat_bnds(i,j,2) = latv
        
        ! NE
        dist = ddist*(j+1.)
        angle = dangle*(i-1.+0.5)
        CALL position_from_point(lon_pt, lat_pt, a, f, elips, dist, angle, lonv, latv)
        lon_bnds(i,j,3) = lonv
        lat_bnds(i,j,3) = latv
  
        ! NW
        dist = ddist*(j-1.)
        angle = dangle*(i-1.+0.5)
        CALL position_from_point(lon_pt, lat_pt, a, f, elips, dist, angle, lonv, latv)
        lon_bnds(i,j,4) = lonv
        lat_bnds(i,j,4) = latv

        ! Getting the distance at the center of the radar cell for height computations
        dist = ddist*(j-1) + ddist2
        CALL position_from_point(lon_pt, lat_pt, a, f, elips, dist, angle, lonv, latv)
        
        lons(i,j) = lonv
        lats(i,j) = latv

        IF (printdbg) PRINT *, '    ' // TRIM(fname), ' i,j ',  i, j, ' <> distance ', dist,          &
          ' angle ', angle, ' lon: ', lons(i,j), ' lat: ', lats(i,j)
                            
        DO k=1, Nazimuth
          ! Centered at half angle
          azimuth = radazimuths(k)

          CALL beam_height(hgt_pt, a, angle, azimuth, dist, hgt, lhgt)
          IF (hgt > max_hgt) CYCLE
          
          hgts(k,j) = hgt
          lhgts(k,j) = lhgt
          
          IF (printdbg) PRINT *, '      ', i, j, k, ' angle:', angle, ' dist:', dist, ' azimuth: ',   &
            azimuth, ' height: ', hgt

          ! Only to be computed at first angle
          IF (i == 1 ) THEN
     
            ! Vertices
            ! Cell boundaries
            ! bottomW
            CALL beam_height(hgt_pt, a, angle, azimuth-dhgt2(k), dist-ddist2, hgt, lhgt)
            hgt_bnds(k,j,1) = hgt
          
            ! bottomE
            CALL beam_height(hgt_pt, a, angle, azimuth-dhgt2(k), dist+ddist2, hgt, lhgt)
            hgt_bnds(k,j,2) = hgt

            ! topE
            CALL beam_height(hgt_pt, a, angle, azimuth+dhgt2(k), dist+ddist2, hgt, lhgt)
            hgt_bnds(k,j,3) = hgt

            ! topW
            CALL beam_height(hgt_pt, a, angle, azimuth+dhgt2(k), dist-ddist2, hgt, lhgt)
            hgt_bnds(k,j,4) = hgt
            
          END IF

        END DO
      END DO
      
      ! debugging at N, W, S, E axes
      IF (dbg) THEN
        IF (i== 1 .OR. i == Nangle/4 .OR. i == Nangle/2 .OR. i == 3*Nangle/4 .OR. i == Nangle) THEN
          PRINT *, '  ' // TRIM(fname) // ": values at major AXIS N;E;S;W angle = ", angle, ' _______'
          PRINT *, '  lon first:', lons(i,1), ' last: ', lons(i,Ndist)
          PRINT *, '  lat first:', lats(i,1), ' last: ', lats(i,Ndist)
          PRINT *, '  hgt at first:', hgts(1,1), ' top: ', hgts(Nazimuth,1)
          PRINT *, '  hgt at last:', hgts(1,Ndist), ' top: ', hgts(Nazimuth,Ndist)        
          PRINT *, '      values idist dist lon lat invangle invdist _______'
          DO j=1, Ndist, Ndbgfrq
            CALL distance_between_points(lon_pt, lat_pt, lons(i,j), lats(i,j), a, f, elips, s, faz)
            PRINT *,'      ',  j, ddist*(j-1) + ddist2, lons(i,j), lats(i,j), faz, s
        
          END DO
        END IF
      END IF

    END DO angl
        
    RETURN

  END SUBROUTINE compute_radar_lonlat_radial

  SUBROUTINE pertinence_polar(iwdx, ewdx, iwdy, ewdy, WRFlon, WRFlat, WRFlon_bnds, WRFlat_bnds,       &
    lonn, lonx, latn, latx, irdx, erdx, irdy, erdy, rlon, rlat, crx, cry, rmax, miss_v, maxNcoords,   &
    distradar, inradar, radarnum, ngridradar, Nc2, ijradar, ijdistradar, inside, idbg)
  ! Subroutine to compute the 2-dimensional pertinence of a lonlat polar radial space (r) into 
  !   another regular lonlat (WRF)
  !    NOTE:
  !      - In radial space
  !         dx: distance          dy: angle
  !         lon(dy,dx), lat(dy,dx)
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: iwdx, ewdx, iwdy, ewdy, maxNcoords, Nc2
    INTEGER, INTENT(in)                                  :: irdx, erdx, irdy, erdy
    REAL, INTENT(in)                                     :: lonn, lonx, latn, latx
    REAL, INTENT(in)                                     :: crx, cry, rmax, miss_v
    REAL, DIMENSION(iwdx:ewdx,iwdy:ewdy), INTENT(in)     :: WRFlon, WRFlat
    REAL, DIMENSION(irdy:erdy,irdx:erdx), INTENT(in)     :: rlon, rlat
    REAL, DIMENSION(iwdx:ewdx,iwdy:ewdy,4), INTENT(in)   :: WRFlon_bnds, WRFlat_bnds
    LOGICAL, OPTIONAL                                    :: idbg
    INTEGER, DIMENSION(iwdx:ewdx,iwdy:ewdy), INTENT(out) :: radarnum
    INTEGER, DIMENSION(iwdx:ewdx,iwdy:ewdy), INTENT(out) :: ngridradar
    INTEGER, DIMENSION(iwdx:ewdx,iwdy:ewdy),                                                          &
      INTENT(inout)                                      :: inradar
    REAL, DIMENSION(iwdx:ewdx,iwdy:ewdy), INTENT(out)    :: distradar
    INTEGER, DIMENSION(iwdx:ewdx,iwdy:ewdy,maxNcoords,Nc2),                                           &
      INTENT(out)                                        :: ijradar
    REAL, DIMENSION(iwdx:ewdx,iwdy:ewdy,maxNcoords),                                                  &
      INTENT(out)                                        :: ijdistradar      
    LOGICAL, INTENT(out)                                 :: inside
    
    ! Local
    INTEGER                                              :: i, j, ir
    INTEGER                                              :: ic, jc, ijfound
    INTEGER, DIMENSION(1)                                :: kmindist
    INTEGER                                              :: rdx2, rdy2, rdx
    INTEGER, DIMENSION(2)                                :: ijtile_rad
    REAL                                                 :: mindist
    REAL                                                 :: lonrad, latrad, distrad
    INTEGER                                              :: totngridradar
    REAL, DIMENSION(iwdx:ewdx,iwdy:ewdy)                 :: ijdist
    LOGICAL, DIMENSION(iwdx:ewdx,iwdy:ewdy)              :: firstngridradar
    LOGICAL                                              :: centerin
    LOGICAL                                              :: dbg, printdbg
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! iwdx, ewdx, iwdy, ewdy: size of the WRF lonlathgt
    ! rnum: number of radar from the list of radars to process
    ! WRFlon, WRFlat: coordinates of the WRF lonlathgt [degrees_east, degrees_North]
    ! WRFlon_bnds, WRFlat_bnds: coordinates of the vertices of the grids of the WRF lonlat
    ! lonn, lonx, latn, latx: actual vertices of the lon, lat values within the patch from its SW, NW, 
    !   NE, SE corners
    ! irdx, erdx, irdy, erdy: size of the polar r lonlathgt that we want to know how much fits inside
    !   WRF lonlat
    ! rlon, rlat: coordinates of the r lonlathgt [degrees_east, degrees_North]
    ! crx, cry: coordinates of the center of the radial space
    ! rmax: maximum radii from the center of the radial space [m]
    ! miss_v: missing value
    ! rmax: maximum possible amount of r lonlathgt cells in a single WRFlonlat grid
    !
    ! OUTPUT
    ! distradar: distance in WRF lonlat grids to the center of r lonlat
    ! inradar: whether a WRF lonlat grid holds inside any r lonlat cell (horizontal only)
    ! radarnum: whether a WRF lonlat grid holds inside the r lonlat cell (horizontal only)
    ! ngridradar: amount of r lonlat cells wihtin a WRF lonlat grid
    ! ijkradar: indices of the polar r lonlat cells wihtin a WRF lonlat grid
    ! ijdistradar: distance to the WRF lonlat grid of the r lonlat cells wihtin a WRF lonlat grid
    fname = 'pertinence_polar'
  
    dbg = set_dbg(idbg)
  
    IF (dbg) THEN
      PRINT *,'  ' // TRIM(fname) // ": looking for radar at ", crx, cry, ' max dist: ', rmax
      PRINT *,'    size of WRF lonlat:', iwdx, ewdx, iwdy, ewdy
      PRINT *,'    size of r lonlat:', irdx, erdx, irdy, erdy
    END IF
      
    rdx2 = (erdx - irdx + 1) / 2
    rdy2 = (erdy - irdy + 1) / 2
  
    ! Getting actual areas encompassed by a radar
    distradar = fillValueR
    inradar = fillValueI
    radarnum = fillValueI
    ngridradar = 0
    ijradar = fillValueI
    ijdistradar = fillValueR
    firstngridradar = .TRUE.
    
    !irng = 0    
    centerin = .FALSE.
    CALL radar_inside (iwdx, ewdx, iwdy, ewdy, WRFlon, WRFlat, WRFlon_bnds, WRFlat_bnds, crx, cry,    & 
      rmax, miss_v, distradar, inradar, radarnum, inside, dbg)
    
    ! Getting lon,lat for the radar matrix
    ! Can we know from the radar perspective which of its grid points are inside this tile?
    ! Even if I could know that, radar matrix is not (by now) parallelized !
    ! Therefore, I need to compute the whole radar matrix at each tile
    
    totngridradar = 0
    ijfound = 1
    !IF (rlon > lonn)) .AND. rlon < lonx) .AND.                     &
    !  rlat > latn) .AND. rlat < latx) THEN
    insideradar: IF (inside) THEN
    
      IF (dbg) PRINT *,'    ' // TRIM(fname) // ": radar inside !!"

      centerin = .TRUE.

      ! Looping all over the radar matrix
      lonlat_radar: DO i=irdy, erdy
        DO j=irdx, erdx
!      lonlat_radar: DO i=1,erdy
!        DO j=irdx, 10
          printdbg = .FALSE.
          IF (dbg) THEN
            IF (MOD(i+j,Ndbgfrq) == 0 .AND. ijfound < Ndbg) THEN
              printdbg = .TRUE.
              ijfound = ijfound + 1
            END IF
          END IF
                         
          ! Avoid values outside radar range
          IF (rlon(i,j) >= fillValueR*0.95 .OR. rlat(i,j) >= fillValueR*0.95) CYCLE

          ins: IF (rlon(i,j) >= lonn .AND. rlon(i,j) < lonx .AND.                                     &
            rlat(i,j) >= latn .AND. rlat(i,j) < latx) THEN
              
            ! Closest WRF grid point to the radar cell
            ijdist=SQRT((WRFlon-rlon(i,j))**2+(WRFlat-rlat(i,j))**2)
            ijtile_rad = MINLOC(ijdist)
            ic = ijtile_rad(1) + iwdx - 1
            jc = ijtile_rad(2) + iwdy - 1
            totngridradar = totngridradar + 1
            
            IF ( (ic < iwdx .OR. ic > ewdx) .OR. (jc < iwdy .OR. jc > ewdy) ) THEN
              CALL ErrMsg("For radar's i,j: " // TRIM(ItoS(i)) //', ' // TRIM(ItoS(j)) //             &
                " wrong WRF grid found!! ic, jc " //  TRIM(ItoS(ic)) //', ' // TRIM(ItoS(jc)) //      &
                " outside patch range its,ite ; jts,jte: " // TRIM(ItoS(iwdx)) // ", " //             &
                TRIM(ItoS(ewdx)) // " ; " // TRIM(ItoS(iwdy)) // ", " // TRIM(ItoS(ewdy)), fname, -1)
            END IF
            
            inradar(ic,jc) = 1
            radarnum(ic,jc) = 1
            
            ! Amount of r cells in this WRF grid point
            !! This is not properly workking, ngridradar = 0 at the begining is not woring ¿?
            IF (firstngridradar(ic,jc)) THEN
              ngridradar(ic,jc) = 1
              firstngridradar(ic,jc) = .FALSE.
            ELSE 
              ngridradar(ic,jc) = ngridradar(ic,jc) + 1
            END IF
            IF (ngridradar(ic,jc) > maxNcoords) THEN
              CALL ErrMsg("For WRF's i,j: " // TRIM(ItoS(ic)) //', ' // TRIM(ItoS(jc)) //             &
                " too much radar cells: " // TRIM(ItoS(ngridradar(ic,jc))) // " in a WRF grid. " //   &
                " distance to the radar: " // TRIM(RtoS(MINVAL(ijdist))) //                           &
                " 'ngrid_radar' in namelist current value: " // TRIM(ItoS(maxNcoords)), fname, -1)
            END IF 
                              
            ! Coordinates of the radar cells within this WRF grid
            ijradar(ic,jc,ngridradar(ic,jc), 1) = i
            ijradar(ic,jc,ngridradar(ic,jc), 2) = j
            ijdistradar(ic,jc,ngridradar(ic,jc)) = MINVAL(ijdist)
               
            IF (printdbg) THEN
              PRINT *,'  ' // TRIM(fname) // ": check total amount of radar cells wihtin " //         &
                "WRF grid: ", totngridradar
              PRINT *, '        > WRF grid:', ijtile_rad
              PRINT *, '        > amount radar cells in WRF grid:', ngridradar(ic,jc)
              PRINT *, '        > radar cell coordinates of last found:', i,j
            END IF

          END IF ins
        END DO
      END DO lonlat_radar

      rdx = erdx-irdx+1
      IF (dbg) THEN
        PRINT *,'    WRF lonlat box ...'
        PRINT *,'      SW:', lonn, latn
        PRINT *,'      NE:', lonx, latx
        PRINT *,'  ' // TRIM(fname) // ": with radar at ", crx,', ',cry,  "; amount of " //           &
          "radar cells in WRF grid: ", totngridradar, ' theoretical value:',                          &
          INT(pict*((rdx)/2)**2), INT(SUM(ngridradar))
        ijtile_rad = MAXLOC(ngridradar)
        ic = ijtile_rad(1) + iwdx - 1
        jc = ijtile_rad(2) + iwdy - 1
        ijdist=SQRT((WRFlon(ic,jc)-crx)**2+(WRFlat(ic,jc)-cry)**2)
        PRINT *,'      largest amount radar cells in a WRF grid:', MAXVAL(ngridradar), ' found at: ', &
          ijtile_rad, ' distance to the radar:', ijdist(iwdx,iwdy), ' degree'
      END IF
        
!! This can not longer be verified ...
!      IF (totngridradar < INT((pict*(rdx/2)**2)/10)) THEN
!        CALL ErrMSg(' found radar cells in WRF grid: ' // TRIM(ItoS(totngridradar)) //                &
!          ' theoretical: ' // TRIM(ItoS(INT(pict*(rdx/2)**2))) // ' Not enough ' //                   &
!          'radar cells were found',fname,-1)
!      END IF
      
    ELSE
      IF (dbg) THEN
        PRINT *,'  ' // TRIM(fname) // ": radar at ", crx,', ', cry, "; outside WRF grid !!"
        PRINT *,'    WRF lonlat box ...'
        PRINT *,'      SW:', lonn, latn
        PRINT *,'      NE:', lonx, latx
      END IF
      
    END IF insideradar

    RETURN
    
  END SUBROUTINE pertinence_polar
     
  SUBROUTINE read_lut(lutfilen, lutgrid, idbg)
  ! Subroutine to read the T-MATRIX look-up-table (LUT) from an specific microphysics from a netCDF 
  !   file
  
    IMPLICIT NONE
  
    CHARACTER(len=S256), INTENT(in)                      :: lutfilen
    LOGICAL, OPTIONAL                                    :: idbg
    TYPE(lut_vars), INTENT(inout)                        :: lutgrid
    
    ! Local
    INTEGER                                              :: lutid, ierr
    INTEGER                                              :: Nlutdims, Nlutvars
    INTEGER                                              :: Nbnds, Nlevs, Nqbins
    CHARACTER(len=S16)                                   :: S16a
    CHARACTER(len=S32)                                   :: S32a
    CHARACTER(len=S256)                                  :: S256a
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! lutfilen: name of the netCDF file with the values
    ! lutgrid: derived type where to kept the variables
    
    fname = 'read_lut'
    
    dbg = set_dbg(idbg)
    
    CALL open_NCfile(lutfilen, lutid, dbg)

    CALL get_Ndims_ncid(lutid, Nlutdims)
    CALL get_Nvars_ncid(lutid, Nlutvars)
    IF (dbg) PRINT *,'  ' // TRIM(fname) // ": getting from file '" // TRIM(lutfilen) //              &
      "' Ndims: ", Nlutdims, ' Nvars: ', Nlutvars

    S32a = 'WRF_mp_physics'
    CALL get_varS0D_ncid(lutid, 256, S32a, S256a, dbg)
    S16a = TRIM(S256a)
    lutgrid%mp_scheme = S16a
    
    ! LUT space
    S32a = 'elevation'
    CALL get_dimlength_ncid(lutid, S32a, Nlevs)
    S32a = 'radar_band'
    CALL get_dimlength_ncid(lutid, S32a, Nbnds)
    S32a = 'qxbin'
    CALL get_dimlength_ncid(lutid, S32a, Nqbins)

    ! I do think it makes the code more readable
    lutgrid%Nlutelevs = Nlevs
    lutgrid%Nbands = Nbnds
    lutgrid%Nlutqxbins = Nqbins
    
    !! Variables
    !!! RAIN
    IF (ALLOCATED(lutgrid%Zh_RAIN)) DEALLOCATE(lutgrid%Zh_RAIN)
    ALLOCATE(lutgrid%Zh_RAIN(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_RAIN', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zh_RAIN'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zh_RAIN, dbg)
    IF (ALLOCATED(lutgrid%Zdr_RAIN)) DEALLOCATE(lutgrid%Zdr_RAIN)
    ALLOCATE(lutgrid%Zdr_RAIN(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_RAIN', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zdr_RAIN'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zdr_RAIN, dbg)
    IF (ALLOCATED(lutgrid%LDR_RAIN)) DEALLOCATE(lutgrid%LDR_RAIN)
    ALLOCATE(lutgrid%LDR_RAIN(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_RAIN', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'LDR_RAIN'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%LDR_RAIN, dbg)
    IF (ALLOCATED(lutgrid%Aih_RAIN)) DEALLOCATE(lutgrid%Aih_RAIN)
    ALLOCATE(lutgrid%Aih_RAIN(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_RAIN', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aih_RAIN'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aih_RAIN, dbg)
    IF (ALLOCATED(lutgrid%Aiv_RAIN)) DEALLOCATE(lutgrid%Aiv_RAIN)
    ALLOCATE(lutgrid%Aiv_RAIN(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_RAIN', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aiv_RAIN'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aiv_RAIN, dbg)
    IF (ALLOCATED(lutgrid%KDP_RAIN)) DEALLOCATE(lutgrid%KDP_RAIN)
    ALLOCATE(lutgrid%KDP_RAIN(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_RAIN', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'KDP_RAIN'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%KDP_RAIN, dbg)

    !!! SNOW
    IF (ALLOCATED(lutgrid%Zh_SNOW)) DEALLOCATE(lutgrid%Zh_SNOW)
    ALLOCATE(lutgrid%Zh_SNOW(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_SNOW', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zh_SNOW'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zh_SNOW, dbg)
    IF (ALLOCATED(lutgrid%Zdr_SNOW)) DEALLOCATE(lutgrid%Zdr_SNOW)
    ALLOCATE(lutgrid%Zdr_SNOW(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_SNOW', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zdr_SNOW'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zdr_SNOW, dbg)
    IF (ALLOCATED(lutgrid%LDR_SNOW)) DEALLOCATE(lutgrid%LDR_SNOW)
    ALLOCATE(lutgrid%LDR_SNOW(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_SNOW', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'LDR_SNOW'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%LDR_SNOW, dbg)
    IF (ALLOCATED(lutgrid%Aih_SNOW)) DEALLOCATE(lutgrid%Aih_SNOW)
    ALLOCATE(lutgrid%Aih_SNOW(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_SNOW', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aih_SNOW'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aih_SNOW, dbg)
    IF (ALLOCATED(lutgrid%Aiv_SNOW)) DEALLOCATE(lutgrid%Aiv_SNOW)
    ALLOCATE(lutgrid%Aiv_SNOW(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_SNOW', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aiv_SNOW'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aiv_SNOW, dbg)
    IF (ALLOCATED(lutgrid%KDP_SNOW)) DEALLOCATE(lutgrid%KDP_SNOW)
    ALLOCATE(lutgrid%KDP_SNOW(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_SNOW', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'KDP_SNOW'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%KDP_SNOW, dbg)

    !!! GRAU
    IF (ALLOCATED(lutgrid%Zh_GRAU)) DEALLOCATE(lutgrid%Zh_GRAU)
    ALLOCATE(lutgrid%Zh_GRAU(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_GRAU', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zh_GRAU'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zh_GRAU, dbg)
    IF (ALLOCATED(lutgrid%Zdr_GRAU)) DEALLOCATE(lutgrid%Zdr_GRAU)
    ALLOCATE(lutgrid%Zdr_GRAU(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_GRAU', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zdr_GRAU'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zdr_GRAU, dbg)
    IF (ALLOCATED(lutgrid%LDR_GRAU)) DEALLOCATE(lutgrid%LDR_GRAU)
    ALLOCATE(lutgrid%LDR_GRAU(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_GRAU', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'LDR_GRAU'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%LDR_GRAU, dbg)
    IF (ALLOCATED(lutgrid%Aih_GRAU)) DEALLOCATE(lutgrid%Aih_GRAU)
    ALLOCATE(lutgrid%Aih_GRAU(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_GRAU', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aih_GRAU'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aih_GRAU, dbg)
    IF (ALLOCATED(lutgrid%Aiv_GRAU)) DEALLOCATE(lutgrid%Aiv_GRAU)
    ALLOCATE(lutgrid%Aiv_GRAU(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_GRAU', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aiv_GRAU'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aiv_GRAU, dbg)
    IF (ALLOCATED(lutgrid%KDP_GRAU)) DEALLOCATE(lutgrid%KDP_GRAU)
    ALLOCATE(lutgrid%KDP_GRAU(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_GRAU', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'KDP_GRAU'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%KDP_GRAU, dbg)

    !!! ICE
    IF (ALLOCATED(lutgrid%Zh_ICE)) DEALLOCATE(lutgrid%Zh_ICE)
    ALLOCATE(lutgrid%Zh_ICE(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_ICE', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zh_ICE'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zh_ICE, dbg)
    IF (ALLOCATED(lutgrid%Zdr_ICE)) DEALLOCATE(lutgrid%Zdr_ICE)
    ALLOCATE(lutgrid%Zdr_ICE(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_ICE', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Zdr_ICE'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Zdr_ICE, dbg)
    IF (ALLOCATED(lutgrid%LDR_ICE)) DEALLOCATE(lutgrid%LDR_ICE)
    ALLOCATE(lutgrid%LDR_ICE(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_ICE', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'LDR_ICE'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%LDR_ICE, dbg)
    IF (ALLOCATED(lutgrid%Aih_ICE)) DEALLOCATE(lutgrid%Aih_ICE)
    ALLOCATE(lutgrid%Aih_ICE(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_ICE', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aih_ICE'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aih_ICE, dbg)
    IF (ALLOCATED(lutgrid%Aiv_ICE)) DEALLOCATE(lutgrid%Aiv_ICE)
    ALLOCATE(lutgrid%Aiv_ICE(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_ICE', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'Aiv_ICE'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%Aiv_ICE, dbg)
    IF (ALLOCATED(lutgrid%KDP_ICE)) DEALLOCATE(lutgrid%KDP_ICE)
    ALLOCATE(lutgrid%KDP_ICE(Nbnds,Nlevs,Nqbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_ICE', 3, (/Nbnds,Nlevs,Nqbins/), fname, ierr)
    S32a = 'KDP_ICE'
    CALL get_varD3D_ncid(lutid, 1, 1, 1, Nbnds, Nlevs, Nqbins, S32a, lutgrid%KDP_ICE, dbg)

! From namelist or from file?
!    S32a = 'q_rain_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_rain_min, dbg)
!    S32a = 'q_rain_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_rain_max, dbg)
!    S32a = 'qn_rain_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_rain_min, dbg)
!    S32a = 'qn_rain_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_rain_max, dbg)
!    S32a = 'q_snow_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_snow_min, dbg)
!    S32a = 'q_snow_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_snow_max, dbg)
!    S32a = 'qn_snow_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_snow_min, dbg)
!    S32a = 'qn_snow_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_snow_max, dbg)
!    S32a = 'q_grau_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_grau_min, dbg)
!    S32a = 'q_grau_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_grau_max, dbg)
!    S32a = 'qn_grau_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_grau_min, dbg)
!    S32a = 'qn_grau_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_grau_max, dbg)
!    S32a = 'q_ice_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_ice_min, dbg)
!    S32a = 'q_ice_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_ice_max, dbg)
!    S32a = 'qn_ice_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_ice_min, dbg)
!    S32a = 'qn_ice_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%qn_ice_max, dbg)
!    S32a = 'q_clou_min'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_clou_min, dbg)
!    S32a = 'q_clou_max'
!    CALL get_varR0D_ncid(lutid, S32a, lutgrid%q_clou_max, dbg)
!    S32a = 'dim_ice'
!    CALL get_varN0D_ncid(lutid, S32a, lutgrid%dim_ice, dbg)
!    S32a = 'dimv'
!    CALL get_varN0D_ncid(lutid, S32a, lutgrid%dimv, dbg)

    !! Getting densities
    S32a = 'rhoo'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%rhoo, dbg)
    S32a = 'rho_r'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%rho_r, dbg)
    S32a = 'n0r'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%n0r, dbg)
    S32a = 'rho_s'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%rho_s, dbg)
    S32a = 'n0s'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%n0s, dbg)
    S32a = 'rho_g'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%rho_g, dbg)
    S32a = 'n0g'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%n0g, dbg)
    S32a = 'rho_h'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%rho_h, dbg)
    S32a = 'n0h'
    CALL get_varR0D_ncid(lutid, S32a, lutgrid%n0h, dbg)

    IF (ALLOCATED(lutgrid%elevations)) DEALLOCATE(lutgrid%elevations)
    ALLOCATE(lutgrid%elevations(Nlevs), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'elevations', 1, (/Nlevs/), fname, ierr)

    S32a = 'elevation'
    CALL get_varR1D_ncid(lutid, 1, Nlevs, S32a, lutgrid%elevations, dbg)

    CALL close_NCfile(lutfilen, lutid, dbg)

    RETURN
 
  END SUBROUTINE read_lut
  
  SUBROUTINE distance_between_points(lon1_pt, lat1_pt, lon2_pt, lat2_pt, a, f, elips, dist, angle)
  ! Subroutine tro compute the coordinates at a given distance and angle from a given point
 
    IMPLICIT NONE
  
    REAL, INTENT(in)                                       :: lon1_pt, lat1_pt, lon2_pt, lat2_pt
    REAL(kind=rr8), INTENT(in)                             :: a, f
    CHARACTER(len=30), INTENT(in)                          :: elips
    REAL, INTENT(out)                                      :: dist, angle

    ! Local
    REAL(kind=rr8)                                         :: ddist
    INTEGER(kind=ii4)                                      :: lon_sign, lat_sign, angle_sign
    REAL(kind=rr8)                                         :: lon_deg, lat_deg, angle_deg, lon_min,     &
      lat_min, angle_min, lon_sec, lat_sec, angle_sec
    REAL(kind=rr8)                                         :: glat1, glon1, faz, distd
    REAL(kind=rr8)                                         :: glat2, glon2, baz
    INTEGER                                                :: Niter, kindv
    REAL(kind=rr8)                                         :: signv, lamv
    CHARACTER(len=S32)                                     :: fname

    !!!!!!! Variables
    ! lon1_pt, lat1_pt: coordinates of the first point [deg]
    ! lon2_pt, lat2_pt: coordinates of the second point [deg]
    ! a: semimajor axis equatorial (in meters)
    ! f: flattening
    ! elips: name of the ellipsis being used
    ! dist: distance [m]
    ! angle: angle [deg]

    fname = 'distance_between_points'
      
    glat1 = lat1_pt / drad
    glon1 = lon1_pt / drad
    glat2 = lat2_pt / drad
    glon2 = lon2_pt / drad
   
    CALL inverse(a, 1.d0/f, glat1, glon1, glat2, glon2, faz, baz, distd, Niter, signv, lamv, kindv)
   
    angle = faz*DegRad
    dist = distd*1.

    RETURN
  
  END SUBROUTINE distance_between_points
        
  SUBROUTINE inverse(a, rf, b1, l1, b2, l2, faz, baz, s, it, sig, lam, kindv)
  ! Subroutine to compute distance between 2 points
  !
  !   Following software by NOAA
  !      https://www.ngs.noaa.gov/TOOLS/Inv_Fwd/Inv_Fwd.html
  !
  !    Bessel, F. W., (1826) On the computation of geographical longitude and latitude from geodetic 
  !      measurements, Astronomische Nachrichten (Astronomical Notes), Band 4, No. 86, Spalten 241-254, 
  !      Altona 1826
  !    Vincenty, T., (1975) Direct and inverse solutions on the ellipsoid with application of nested 
  !      equations, Survey Review, Vol. 23, No. 176, pp. 88-93.
  !    Vincenty, T., (1976) Correspondence: solutions of geodesics, survey Review, Vol. 23, No. 180, 
  !      p. 294.

  !
  !*** inverse for long-line and antipodal cases.
  !*** latitudes may be 90 degrees exactly.
  !*** latitude positive north, longitude positive east, radians.
  !*** azimuth clockwise from north, radians.
  !*** original programmed by thaddeus vincenty, 1975, 1976
  !*** removed back side solution option, debugged, revised -- 2011may01 -- dgm
  !*** this version of code is interim -- antipodal boundary needs work

  !*** output (besides faz,baz, and s):
  !*** it,   iteration count
  !*** sig,  spherical distance on auxiliary sphere
  !*** lam,  longitude dIFference on auxiliary sphere
  !*** kind, solution flag:  kind=1, long-line;  kind=2, antipodal

    IMPLICIT NONE

    REAL(kind=rr8), INTENT(in)                           :: a, rf 
    REAL(kind=rr8), INTENT(in)                           :: b1, b2, l1, l2
    REAL(kind=rr8), INTENT(out)                          :: faz, baz, s
    INTEGER, INTENT(out)                                 :: it, kindv
    REAL(kind=rr8), INTENT(out)                          :: sig, lam
      
    ! Local
    INTEGER                                              :: i, iter, Niter
    REAL(kind=rr8)                                       :: boa, l
    REAL(kind=rr8)                                       :: biga, bigb, bige, bigf, c, d, z
    REAL(kind=rr8)                                       :: cosal2, coslam, cossig, costm, costm2
    REAL(kind=rr8)                                       :: dsig, ep2
    REAL(kind=rr8)                                       :: prev, tem1, tem2, temp, test
    REAL(kind=rr8)                                       :: sinal, sinlam, sinsig
    REAL(kind=rr8)                                       :: beta1, sinu1, cosu1
    REAL(kind=rr8)                                       :: beta2, sinu2, cosu2
    LOGICAL                                              :: converged, keeplooping
    CHARACTER(len=S32)                                   :: fname
       
    fname = 'inverse'

    boa =  1.d0-1.d0/rf
    !***** sinu1 = boa*dsin(b1)/dsqrt((boa*dsin(b1))**2+dcos(b1)**2)
    !***** cosu1 = dsqrt(-sinu1**2+1.d0)                               !*** rounDOff
    !***** sinu2 = boa*dsin(b2)/dsqrt((boa*dsin(b2))**2+dcos(b2)**2)
    !***** cosu2 = dsqrt(-sinu2**2+1.d0)                               !*** rounDOff

    !*** better reduced latitude
    beta1 = datan(boa*dtan(b1))             
    sinu1 = dsin(beta1)
    cosu1 = dcos(beta1)

    !*** better reduced latitude
    beta2 = datan(boa*dtan(b2))             
    sinu2 = dsin(beta2)
    cosu2 = dcos(beta2)

    !*** longitude dIFference [-pi,pi]
    l = l2-l1                               
    IF (l .gt.  dpict) l = l - dpict - dpict
    IF (l .lt. -dpict) l = l + dpict + dpict
    prev = l
    test = l
    kindv = 1
    !*** v13   (rapp part II)
    lam = l

    !*** top of the long-line loop (kind=1)
    it = 0
    Niter = 100

    DO i=1, Niter

      sinlam = dsin(lam)

      !*** no--troublesome
      !***** IF(dabs(pi-dabs(l)).lt.0.2d-10) sinlam=0.d0
      coslam = dcos(lam)
      temp = cosu1*sinu2-sinu1*cosu2*coslam
      !*** v14   (rapp part II)
      sinsig = dsqrt((cosu2*sinlam)**2+temp**2)
      !*** v15   (rapp part II)
      cossig = sinu1*sinu2 + cosu1*cosu2*coslam
      sig = datan2(sinsig,cossig)

      !*** v17   (rapp part II)
      !***** sinal = cosu1*cosu2*sinlam/sinsig
      IF(dabs(sinsig).lt.eps) THEN
        !*** avoid div 0
        IF (sinsig > 0.) THEN
          sinal = cosu1*cosu2*sinlam/eps
        ELSE IF (sinsig < 0.) THEN
          sinal = -cosu1*cosu2*sinlam/eps
        ELSE
          sinal = cosu1*cosu2*sinlam/eps
        END IF
        !sinal = cosu1*cosu2*sinlam/dsign(eps,sinsig)
      ELSE
        sinal = cosu1*cosu2*sinlam/sinsig
      ENDIF

      cosal2 = -sinal**2+1.d0
      !*** v18   (rapp part II)
      !***** costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig
      IF(dabs(cosal2).lt.eps) THEN
        !*** avoid div 0
        IF (cosal2 > 0.) THEN
          costm = -2.d0*(sinu1*sinu2/eps)+cossig
        ELSE IF (cosal2 < 0.) THEN
          costm = +2.d0*(sinu1*sinu2/eps)+cossig
        ELSE
          costm = -2.d0*(sinu1*sinu2/eps)+cossig
        END IF
        !costm = -2.d0*(sinu1*sinu2/dsign(eps,cosal2))+cossig
      ELSE
        costm = -2.d0*(sinu1*sinu2/cosal2)+cossig
      ENDIF
      costm2=costm*costm
      !*** v10   (rapp part II)
      c=((-3.d0*cosal2+4.d0)/rf+4.d0)*cosal2/rf/16.d0

      converged = .FALSE.
      DO it=1, Niter
        
        !*** v11
        d=(((2.d0*costm2-1.d0)*cossig*c+costm)*sinsig*c+sig)*(1.d0-c)/rf
        keeplooping = .FALSE.

        IF (kindv .eq. 1) THEN
          lam = l+d*sinal
          IF (DABS(lam-test) .lt. tol) THEN
            converged = .TRUE.
            EXIT
            
          END IF

          IF (DABS(lam) .gt. dpict) THEN
          
            kindv = 2
            lam = dpict
            IF (l .lt. 0.d0) lam = -lam
            sinal = 0.d0
            cosal2 = 1.d0
            test = 2.d0
            prev = test
            sig = dpict - DABS(DATAN(sinu1/cosu1) + DATAN(sinu2/cosu2))
            sinsig = DSIN(sig)
            cossig = DCOS(sig)
            c = ((-3.d0*cosal2+4.d0)/rf+4.d0)*cosal2/rf/16.d0
            IF (DABS(sinal-prev) .lt. tol) THEN
              converged = .TRUE.
              EXIT
            END IF
            !*****     costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig
            IF (DABS(cosal2) .lt. eps) THEN
              !*** avoid div 0
              IF (cosal2 > 0.) THEN
                costm= -2.d0*(sinu1*sinu2/eps)+cossig
              ELSE IF (cosal2 < 0.) THEN
                costm= +2.d0*(sinu1*sinu2/eps)+cossig
              ELSE
                costm= -2.d0*(sinu1*sinu2/eps)+cossig
              END IF
              !costm= -2.d0*(sinu1*sinu2/dsign(eps,cosal2))+cossig
            ELSE
              costm= -2.d0*(sinu1*sinu2/cosal2)+cossig
            ENDIF
            costm2=costm*costm            
            keeplooping = .TRUE.
          ENDIF
          
          IF ( .NOT. keeplooping) THEN

            IF ( ((lam-test)*(test-prev)) .lt. 0.d0 .AND. it .gt. 5) lam = (2.d0*lam+3.d0*test +      &
              prev) / 6.d0

            prev = test
            test = lam

            EXIT
            
          END IF

        ELSE
          KEEPLOOPING = .TRUE.
    
          sinal= (lam-l)/d
          !*** refined converge.
          IF ( ((sinal-test)*(test-prev)) .lt. 0.d0 .AND. it.gt.5) sinal = (2.d0*sinal+3.d0*test +    &
            prev)/6.d0

          prev = test
          test = sinal

          cosal2 = -sinal**2+1.d0
          sinlam = sinal*sinsig/(cosu1*cosu2)
          coslam = -DSQRT(DABS(-sinlam**2 +1.d0))
          lam = DATAN2(sinlam,coslam)
          temp = cosu1*sinu2 - sinu1*cosu2*coslam
          sinsig = DSQRT((cosu2*sinlam)**2+temp**2)
          cossig = sinu1*sinu2 + cosu1*cosu2*coslam
          sig = DATAN2(sinsig,cossig)
          c = ((-3.d0*cosal2+4.d0)/rf+4.d0)*cosal2/rf/16.d0
          IF (DABS(sinal-prev) .lt. tol) THEN
            converged = .TRUE.
            EXIT
          END IF

          !*****   costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig
          IF (DABS(cosal2) .lt. eps) THEN
            !*** avoid div 0
            IF (cosal2 > 0.) THEN
              costm = -2.d0*(sinu1*sinu2/eps) + cossig
            ELSE IF (cosal2 < 0.) THEN
              costm = +2.d0*(sinu1*sinu2/eps) + cossig
            ELSE
              costm = -2.d0*(sinu1*sinu2/eps) + cossig
            END IF
            !costm = -2.d0*(sinu1*sinu2/dsign(eps,cosal2)) + cossig
          ELSE
            costm = -2.d0*(sinu1*sinu2/cosal2)+cossig
          ENDIF

          costm2=costm*costm

        ENDIF
    
      END DO

      IF (.NOT.keeplooping) EXIT

      IF (converged) EXIT
    
    END DO

    !*** convergence
    IF (kindv .eq. 2) THEN                   
      !*** antipodal
      faz= sinal/cosu1
      baz= DSQRT(-faz**2+1.d0)
      IF (temp .lt. 0.d0) baz = -baz
      faz = DATAN2(faz,baz)
      tem1 = -sinal
      tem2 = sinu1*sinsig - cosu1*cossig*baz
      baz = DATAN2(tem1,tem2)
    ELSE
      !*** long-line
      tem1 = cosu2*sinlam
      tem2 = cosu1*sinu2 - sinu1*cosu2*coslam
      faz = DATAN2(tem1,tem2)
      tem1 = -cosu1*sinlam
      tem2 = sinu1*cosu2 - cosu1*sinu2*coslam
      baz = DATAN2(tem1,tem2)
    ENDIF

    IF (faz .lt. 0.d0) faz = faz + dpict + dpict
    IF (baz .lt. 0.d0) baz = baz + dpict + dpict

    !*** Helmert 1880 from Vincenty "Geodetic inverse solution between antipodal points"
    ep2 = 1.d0/(boa*boa)-1.d0
    !*** 15
    bige = DSQRT(1.d0+ep2*cosal2)
    !*** 16
    bigf = (bige-1.d0)/(bige+1.d0)
    !*** 17
    biga = (1.d0+bigf*bigf/4.d0)/(1.d0-bigf)
    !*** 18
    bigb = bigf*(1.d0-0.375d0*bigf*bigf)

    z =bigb/6.d0*costm*(-3.d0+4.d0*sinsig**2)*(-3.d0+4.d0*costm2)

    !*** 19
    dsig = bigb*sinsig*(costm+bigb/4.d0*(cossig*(-1.d0+2.d0*costm2)-z))
    !*** 20
    s = (boa*a)*biga*(sig-dsig)

    RETURN

  END SUBROUTINE inverse

  SUBROUTINE position_from_point(lon_pt, lat_pt, a, f, elips, dist, angle, lon_dist, lat_dist)
  ! Subroutine to compute the coordinates at a given distance and angle from a given point
  !
  !   Following software by NOAA
  !      https://www.ngs.noaa.gov/TOOLS/Inv_Fwd/Inv_Fwd.html
  !
  !    Bessel, F. W., (1826) On the computation of geographical longitude and latitude from geodetic 
  !      measurements, Astronomische Nachrichten (Astronomical Notes), Band 4, No. 86, Spalten 241-254, 
  !      Altona 1826
  !    Vincenty, T., (1975) Direct and inverse solutions on the ellipsoid with application of nested 
  !      equations, Survey Review, Vol. 23, No. 176, pp. 88-93.
  !    Vincenty, T., (1976) Correspondence: solutions of geodesics, survey Review, Vol. 23, No. 180, 
  !      p. 294.
  
    IMPLICIT NONE
  
    REAL, INTENT(in)                                     :: lon_pt, lat_pt, dist, angle
    REAL(kind=rr8), INTENT(in)                           :: a, f
    CHARACTER(len=30), INTENT(in)                        :: elips
    REAL, INTENT(out)                                    :: lat_dist, lon_dist

    ! Local
    REAL(kind=rr8)                                       :: ddist
    INTEGER(kind=ii4)                                    :: lon_sign, lat_sign, angle_sign
    REAL(kind=rr8)                                       :: lon_deg, lat_deg, angle_deg, lon_min,     &
      lat_min, angle_min, lon_sec, lat_sec, angle_sec
    REAL(kind=rr8)                                       :: glat1, glon1, faz
    REAL(kind=rr8)                                       :: glat2, glon2, baz
    INTEGER(kind=ii4)                                    :: glon2_sign, glat2_sign
    REAL(kind=rr8)                                       :: dms_angle
    INTEGER(kind=ii4)                                    :: glon2_deg, glat2_deg, glon2_min,          &
      glat2_min
    REAL(kind=rr8)                                       :: glon2_sec, glat2_sec
    REAL(kind=rr8)                                       :: lat2, lon2
    CHARACTER(len=S32)                                   :: fname

    !!!!!!! Variables
    ! lon_pt, lat_pt: coordinates of the point
    ! a: semimajor axis equatorial (in meters)
    ! f: flattening
    ! elips: name of the ellipsis being used
    ! dist: distance [m]
    ! angle: angle [deg]
    fname = 'position_from_point'
      
    glat1 = lat_pt / drad
    glon1 = lon_pt / drad
    faz = angle / drad
  
    ddist = dist*1.d0
  
    ! Assuming we are well behaved...
    !  ! Checking distance
    !  IF ( ddist > dd_max ) THEN
    !    PRINT (errormsg)
    !    PRINT *, ' Invalid Distance: ', ddist, ' ... Please re-enter it '
    !    PRINT *, "   maximum distance for projection '" // TRIM(elips) // "' is ", dd_max, ' m'
    !    STOP -1
    !  END IF

    !PRINT *,'  ' // TRIM(fname) // ": ", a, f, glat1, glon1, ddist
    CALL dirct1(a, f, glat1, glon1, glat2, glon2, faz, baz, ddist)
   
    lat_dist = glat2*drad
    lon_dist = glon2*drad
  
    RETURN
  
  END SUBROUTINE position_from_point

  SUBROUTINE DIRCT1(A,F,GLAT1,GLON1,GLAT2,GLON2,FAZ,BAZ,S)
!
! *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY
! *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
! *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
!
! *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
! *** F IS THE FLATTENING OF THE REFERENCE ELLIPSOID
! *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
! *** AZIMUTHS IN RADIANS CLOCKWISE FROM NORTH
! *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A
!
! *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75
! *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608
!
!   Following software by NOAA
!      https://www.ngs.noaa.gov/TOOLS/Inv_Fwd/Inv_Fwd.html
!
!    Bessel, F. W., (1826) On the computation of geographical longitude and latitude from geodetic 
!      measurements, Astronomische Nachrichten (Astronomical Notes), Band 4, No. 86, Spalten 241-254, 
!      Altona 1826
!    Vincenty, T., (1975) Direct and inverse solutions on the ellipsoid with application of nested 
!      equations, Survey Review, Vol. 23, No. 176, pp. 88-93.
!    Vincenty, T., (1976) Correspondence: solutions of geodesics, survey Review, Vol. 23, No. 180, 
!      p. 294.

    IMPLICIT NONE

    REAL(kind=rr8), INTENT(in)                           :: GLAT1, GLON1, FAZ, A, F, S
    REAL(kind=rr8), INTENT(out)                          :: GLAT2, GLON2, BAZ

    ! Local
    INTEGER                                              :: i
    REAL(kind=rr8)                                       :: C, D, EPS
    REAL(kind=rr8)                                       :: R, TU, SF, CF, CU, SU, SA, C2A, X, E 
    REAL(kind=rr8)                                       :: Y, SY, CY, CZ
    LOGICAL                                              :: converged
    
    CHARACTER(len=S32)                                   :: fname

    DATA EPS /0.5D-13/
    
    fname = 'DIRCT1'

    R=1.-F
    TU=R*DSIN(GLAT1)/DCOS(GLAT1)
    SF=DSIN(FAZ)
    CF=DCOS(FAZ)
    BAZ=0.
    IF (CF.NE.0.) BAZ=DATAN2(TU,CF)*2.

    CU=1./DSQRT(TU*TU+1.)
    SU=TU*CU
    SA=CU*SF
    C2A=-SA*SA+1.
    X=DSQRT((1./R/R-1.)*C2A+1.)+1.
    X=(X-2.)/X
    C=1.-X
    C=(X*X/4.+1)/C
    D=(0.375D0*X*X-1.)*X
    TU=S/R/A/C
    Y=TU
    
    converged = .FALSE.
    DO i=1, 100
      SY=DSIN(Y)
      CY=DCOS(Y)
      CZ=DCOS(BAZ+Y)
      E=CZ*CZ*2.-1.
      C=Y
      X=E*CY
      Y=E+E-1.
      Y=(((SY*SY*4.-3.)*Y*CZ*D/6.+X)*D/4.-CZ)*SY*D+TU
      IF (DABS(Y-C).LE.EPS) THEN
        converged = .TRUE.
        EXIT
      END IF
    END DO
    IF (.NOT.converged) THEN
      PRINT *, TRIM(errormsg)
      PRINT *, TRIM(fname) // ': computation of Y,C did not converged !!'
      PRINT *, '    DABS(Y-C)> EPS:', DABS(Y-C), ' > ', EPS
      STOP -1
    END IF

    BAZ=CU*CY*CF-SU*SY
    C=R*DSQRT(SA*SA+BAZ*BAZ)
    D=SU*CY+CU*SY*CF
    GLAT2=DATAN2(D,C)
    C=CU*CY-SU*SY*CF
    X=DATAN2(SY*SF,C)
    C=((-3.*C2A+4.)*F+4.)*C2A*F/16.
    D=((E*CY*C+CZ)*SY*C+Y)*SA
    GLON2=GLON1+X-(1.-C)*D*F
    BAZ=DATAN2(SA,BAZ)+dpict
    
    RETURN

  END SUBROUTINE DIRCT1
  
  SUBROUTINE radar_inside (ims, ime, jms, jme, lon, lat, bnds_lon, bnds_lat, rlon, rlat, rdist,       &
    missing_value, gatehdist, inradar, radnum, inside, idbg)
  ! Subroutine to determine if a radar fells inside the tile
  
    IMPLICIT NONE
  
    INTEGER, INTENT(in)                                  :: ims, ime, jms, jme
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: lon, lat
    REAL, DIMENSION(ims:ime,jms:jme,4), INTENT(in)       :: bnds_lon, bnds_lat
    REAL, INTENT(in)                                     :: rlon, rlat, rdist, missing_value
    LOGICAL, OPTIONAL                                    :: idbg
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(out)        :: gatehdist
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(out)     :: inradar, radnum
    LOGICAL, INTENT(out)                                 :: inside
  
    ! Local
    INTEGER                                              :: i, j
    INTEGER                                              :: ii, ie, ji, je
    INTEGER                                              :: ibnd
    LOGICAL, DIMENSION(ims:ime,jms:jme)                  :: mask
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! lon, lat: 2D array longitude, latitude matrices at the center 
    ! bnds_lon, bnds_lat: 2D array longitude, latitude matrices of the vertices
    ! rlon, rlat: location of the radar
    ! rdist: maximum distance to the radar
    ! gatehdist: distance to the radar in the model grid
    ! inradar: whether the model grid lays inside the measurement area of any radar
    ! radnum: pertinence of the model grid to a given radar
    fname = 'radar_inside'

    dbg = set_dbg(idbg)
  
    ! Is this radar inside the range of the tile?
    inside = .FALSE.
    
    gatehdist = missing_value
    inradar = 0
    radnum = fillValueI
      
    ! Do we need to be efficient? It is not working....
    ii = ims
    ie = ime
    ji = jms
    je = jme
    CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside, dbg)

    ! Going to the borders
    IF (.NOT.inside) THEN
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside, dbg)  
      END DO
    END IF

    !! S border
    IF (.NOT.inside) THEN
      ii = ims
      ie = ime
      ji = jms
      je = jms
      CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside, dbg)

      ! Going to the borders
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside, dbg)  
      END DO
    END IF
    
    !! W border
    IF (.NOT.inside) THEN
      ii = ims
      ie = ims
      ji = jms
      je = jme
      CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside, dbg)

      ! Going to the borders
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside, dbg)  
      END DO
    END IF
    
    !! N border
    IF (.NOT.inside) THEN
      ii = ims
      ie = ime
      ji = jme
      je = jme
      CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside, dbg)

      ! Going to the borders
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside, dbg)  
      END DO
    END IF
        
    !! E border
    IF (.NOT.inside) THEN
      ii = ime
      ie = ime
      ji = jms
      je = jme
      CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside, dbg)

      ! Going to the borders
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside, dbg)  
      END DO
    END IF

    ! Filling other values with a single loop
    IF (inside) THEN
      CALL real_distances(ims, ime, jms, jme, lon, lat, rlon, rlat, gatehdist)
      DO i=ims, ime
        DO j=jms, jme
          IF (gatehdist(i,j) <= rdist) THEN
            inradar(i,j) = 1
            radnum(i,j) = 1
          ELSE
            gatehdist(i,j) = missing_value
          END IF
        END DO
      END DO
      mask = radnum == 1
      IF (dbg) PRINT *,'    ' // TRIM(fname) // ": found ", SUM(radnum, mask=mask), " WRF grids inside"
    END IF

    RETURN

  END SUBROUTINE radar_inside
  
  SUBROUTINE is_inside(ims, ime, jms, jme, lon, lat, rlon, rlat, rdist, inside, idbg)
  ! Subroutine to check if a tile gets a radar
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: ims, ime, jms, jme
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: lon, lat
    REAL, INTENT(in)                                     :: rlon, rlat, rdist
    LOGICAL, OPTIONAL                                    :: idbg
    LOGICAL, INTENT(inout)                               :: inside
  
    ! Local
    REAL                                                 :: mindist
    REAL, DIMENSION(:,:), ALLOCATABLE                    :: distborder  
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! lon, lat: 2D array longitude, latitude matrices at the center 
    ! rlonv, rlat: location of the radar
    ! rdist: maximum distance to the radar
    ! gatehdist: distance to the radar in the model grid
    fname = 'is_inside'
    
    dbg = set_dbg(idbg)
      
    IF (ALLOCATED(distborder)) DEALLOCATE(distborder)
    ALLOCATE(distborder(ims:ime,jms:jme))
    CALL real_distances(ims, ime, jms, jme,                                                           &
      lon(ims:ime,jms:jme), lat(ims:ime,jms:jme),                                                     &
      rlon, rlat,                                                                                     &
      distborder)
    mindist = MINVAL(distborder)
   
    IF (mindist < rdist) inside = .TRUE.

!    IF (dbg) PRINT *, '   ' //TRIM(fname)//' for ims,ime;jms,jme ',ims,' , ',ime,' ; ',jms,' , ',jme, &
!      ' min,max lon ', MINVAL(lon), ',', MAXVAL(lon), ' min,max lat ', MINVAL(lat), ',', MAXVAL(lat), &
!      ' mindist: ', mindist, ' inside ', inside

    IF (ALLOCATED(distborder)) DEALLOCATE(distborder)
    
    RETURN

  END SUBROUTINE is_inside

  SUBROUTINE beam_height(hgt_pt, a, angle, azimuth, dist, hgt, lhgt)
  ! Function to compute the height of the beam
  
  ! After
  !
  ! [PURPOSE:] Common radar procedures and variables.
  !
  ! [HISTORY:]
  !   02/19/2014 Juan Ruiz created
  !

    IMPLICIT NONE
    
    REAL, INTENT(in)                                     :: hgt_pt, angle, azimuth, dist
    REAL(kind=rr8), INTENT(in)                           :: a
    REAL, INTENT(out)                                    :: hgt, lhgt
    
    ! Local
    REAL                                                 :: tmp1, tmp2
    REAL(kind=rr8)                                       :: height, local_elevation
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! hgt_pt: height of the radar station [m]
    ! a: radar ellipsoide semimajor axis equatorial [m]
    ! angle: angle of the beam [rad]
    ! azimuth: azimuth of the beam [rad]
    ! dist: distance respect the radar [m]
    fname = 'beam_height'
    
    ! Initializing
    height = 0.d0
    local_elevation = 0.d0
  
    tmp1 = dist**2 + (ke*a)**2 + 2.*dist*ke*a*sin(azimuth)
    height = hgt_pt + SQRT( tmp1 ) - ke*a

    !Compute the local elevation angle tacking into account the efective earth radius.
    tmp1= dist * cos( azimuth )
    tmp2= dist * sin( azimuth ) + ke*a
    local_elevation = azimuth + atan(tmp1/tmp2)
   
    IF (height < 0.) THEN
      PRINT *,'  ' // TRIM(ErrorMsg)
      PRINT *,'    dist: ', dist, ' azimuth: ', azimuth
      PRINT *,'    radar height: ', height
      CALL ErrMsg('Negative radar height !', fname, -1)
    END IF
   
    hgt = REAL(height)
    lhgt = REAL(local_elevation)
   
   RETURN

  END SUBROUTINE beam_height

  SUBROUTINE real_distances_ellipsoide(ims, ime, jms, jme, a, f, elips                                &
    ,lon, lat                                                                                         &
    ,lonv, latv                                                                                       &
    ,dist, idbg)
!    ,dlonpt, dlatpt)
    ! Subroutine to compute the diferential of coordinattes of a distance from a given point on an 
    !   ellipsoide
    !   Following software by NOAA:
    !      https://www.ngs.noaa.gov/TOOLS/Inv_Fwd/Inv_Fwd.html
    !
    !    Bessel, F. W., (1826) On the computation of geographical longitude and latitude from geodetic 
    !      measurements, Astronomische Nachrichten (Astronomical Notes), Band 4, No. 86, Spalten 241-254, 
    !      Altona 1826
    !    Vincenty, T., (1975) Direct and inverse solutions on the ellipsoid with application of nested 
    !      equations, Survey Review, Vol. 23, No. 176, pp. 88-93.
    !    Vincenty, T., (1976) Correspondence: solutions of geodesics, survey Review, Vol. 23, No. 180, 
    !      p. 294.
    
    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ims, ime, jms, jme
    REAL(kind=rr8), INTENT(in)                           :: a, f
    CHARACTER(len=30), INTENT(in)                        :: elips
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(in)     :: lon, lat
    REAL, INTENT(in)                                     :: lonv, latv
    LOGICAL, OPTIONAL                                    :: idbg
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(out)    :: dist
!    REAL, INTENT(out)                                    :: dlonpt, dlatpt

    ! Local 
    INTEGER                                              :: i,j
    REAL                                                 :: angleab
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! lon, lat: 2D array longitude, latitude matrices
    ! lonv, latv: point to compute the distance
    ! a: semimajor axis equatorial (in meters)
    ! f: flattening
    ! elips: name of the ellipsis being used

    fname = 'real_distance_ellipsoide'
  
    dbg = set_dbg(idbg)
      
    dist = 0.

    DO i=ims, ime
      DO j=jms, jme
        CALL distance_between_points(lonv, latv, lon(i,j), lat(i,j), a, f, elips, dist(i,j), angleab)
      END DO
    END DO
   
    RETURN
 
  END SUBROUTINE real_distances_ellipsoide
  
  SUBROUTINE real_distances(ims, ime, jms, jme                                                        &
    ,lon, lat                                                                                         &
    ,lonv, latv                                                                                       &
    ,dist)
!    ,dlonpt, dlatpt)
    ! Subroutine to compute the diferential of coordinattes of a distance from a given point
    !   Following: https://en.wikipedia.org/wiki/Geographical_distance
    !
    !     D = 2*R * SQRT( [sin(dlat /2)*cos(dlon/2)]**2 + [cos(latm)*sin(dlon/2)]**2 )
    
    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ims, ime, jms, jme
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(in)     :: lon, lat
    REAL, INTENT(in)                                     :: lonv, latv
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(out)    :: dist
!    REAL, INTENT(out)                                    :: dlonpt, dlatpt

    ! Local 
!    REAL, PARAMETER                                      :: R = 6371.009 ! km
    REAL, PARAMETER                                      :: R = 6371009.
    REAL, DIMENSION(ims:ime,jms:jme)                     :: dlat, dlon, latm, sec1, sec2
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! lon, lat: 2D array longitude, latitude matrices
    ! lonv, latv: point to compute the distance
    fname = 'real_distance'
  
    dist = 0.
    
    dlat = (latv - lat)/DegRad
    dlon = (lonv - lon)/DegRad
    latm = (lat + latv)/DegRad / 2.
    
    sec1 = sin(dlat /2.)*cos(dlon/2.)
    sec2 = cos(latm)*sin(dlon/2.)
    dist = 2.*R * SQRT( (sec1)**2 + (sec2)**2 )
   
    RETURN
 
  END SUBROUTINE real_distances

!!!! basic
  
  LOGICAL FUNCTION set_dbg(idbg)
  ! Function to setup debug value
  
    IMPLICIT NONE
    
    LOGICAL, OPTIONAL                                    :: idbg
    
    set_dbg = .FALSE.
    
    IF (.NOT.PRESENT(idbg)) THEN
      set_dbg = .FALSE.
    ELSE
      set_dbg = idbg
    END IF
    
    RETURN
      
  END FUNCTION set_dbg

  SUBROUTINE ErrMsgIAvail(msg, Navail, avail, funcn)
  ! Subroutine to stop execution by non exsitent integer value from available ones and provide an 
  !   error message

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in)                           :: msg, funcn
    INTEGER, INTENT(in)                                    :: Navail
    INTEGER, DIMENSION(Navail), INTENT(in)                 :: avail

    ! Local
    INTEGER                                                :: i
    CHARACTER (LEN=1000)                                   :: diag_message

    !!!!!!! Variables
    ! msg: message to print with the error
    ! funcn: name of the funciton, section to localize the error
    ! Navail: amount of available values
    ! avail: available values

    WRITE (diag_message , * ) &
      TRIM(errormsg) // ":   " // TRIM(funcn) // ": " // TRIM(msg) //                                 &
      ' not ready !!  Available ones: ', (TRIM(ItoS(avail(i))) // ', ', i=1, Navail)
      CALL wrf_error_fatal ( diag_message )

    RETURN

  END SUBROUTINE ErrMsgIAvail

  SUBROUTINE ErrMsgSAvail(msg, Navail, avail, funcn)
  ! Subroutine to stop execution by non exsitent string value from available ones and provide an error 
  !   message

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in)                           :: msg, funcn
    INTEGER, INTENT(in)                                    :: Navail
    CHARACTER(len=*), DIMENSION(Navail), INTENT(in)        :: avail

    ! Local
    INTEGER                                                :: i
    CHARACTER (LEN=1000)                                   :: diag_message

    !!!!!!! Variables
    ! msg: message to print with the error
    ! funcn: name of the funciton, section to localize the error
    ! Navail: amount of available values
    ! avail: available values

    WRITE (diag_message , * ) &
      TRIM(errormsg) // ":   " // TRIM(funcn) // ": " // TRIM(msg) //                                 &
      ' not ready !!  Available ones: ', ("'" // TRIM(avail(i)) // "', ", i=1, Navail)
      CALL wrf_error_fatal ( diag_message )

    RETURN

  END SUBROUTINE ErrMsgSAvail

  SUBROUTINE ErrMsgAlloc(msg, arrayn, rank, shapev, funcn, errN)
  ! Subroutine to stop execution and provide an error message when allocation problems

    IMPLICIT NONE

    INTEGER, INTENT(in)                                    :: errN, rank
    CHARACTER(LEN=*), INTENT(in)                           :: msg, funcn, arrayn
    INTEGER, DIMENSION(rank), INTENT(in)                   :: shapev

    ! Local
    INTEGER                                                :: i
    CHARACTER (LEN=1000)                                   :: diag_message

    !!!!!!! Variables
    ! msg: message to print with the error
    ! rank: rank of the array to allocate
    ! shapev: shape of the array to allocate
    ! funcn: name of the funciton, section to localize the error
    ! errN: number of the error provided for a given FORTRAN function

    IF (errN /= 0) THEN
      WRITE (diag_message , * ) &
      TRIM(errormsg) // ":   " // TRIM(funcn) // ": " // TRIM(msg) //                                 &
      ' Problem allocating : ' // "'" // TRIM(arrayn) // "' (", (TRIM(ItoS(shapev(i))) // ', ',i=1,   &
        rank), " )' !!"
      CALL wrf_error_fatal ( diag_message )
    END IF


    RETURN

  END SUBROUTINE ErrMsgAlloc

  SUBROUTINE ErrMsg(msg, funcn, errN)
  ! Subroutine to stop execution and provide an error message

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in)                         :: msg, funcn
    INTEGER, INTENT(in)                                  :: errN

    ! Local
    CHARACTER (LEN=S2048)                                :: diag_message

    !!!!!!! Variables
    ! msg: message to print with the error
    ! funcn: name of the funciton, section to localize the error
    ! errN: number of the error provided for a given FORTRAN function

    IF (errN /= 0) THEN
      PRINT *,TRIM(errormsg)
      PRINT *,'  ' // TRIM(funcn) // ': ' // TRIM(msg) // ' !!'
      PRINT *,'    error number:', errN
      STOP

      WRITE (diag_message , * ) &
      TRIM(errormsg) // ":   " // TRIM(funcn) // ": " // TRIM(msg) // ' !!'
      CALL wrf_error_fatal ( diag_message )
    END IF

    RETURN

  END SUBROUTINE ErrMsg

  CHARACTER(len=S32) FUNCTION DtoS(Dval)
  ! DtoS: Function to transform a double real to String

    IMPLICIT NONE

    REAL(kind=rr8), INTENT(IN)                           :: Dval

    WRITE(DtoS,'(F32.10)')Dval

  END FUNCTION DtoS

  CHARACTER(len=S32) FUNCTION RtoS(Rval)
  ! RtoS: Function to transform a real to String

    IMPLICIT NONE

    REAL, INTENT(IN)                                     :: Rval

    WRITE(RtoS,'(F32.10)')Rval

  END FUNCTION RtoS

  CHARACTER(len=S32) FUNCTION ItoS(Ival)
  ! ItoS: Function to transform an integer to String

    IMPLICIT NONE

    INTEGER, INTENT(IN)                                  :: Ival

    ! Local
    CHARACTER(len=S32)                                   :: itoS0

    WRITE(ItoS0,'(I32)')Ival
    CALL  removeNONnum(ItoS0, ItoS)

  END FUNCTION ItoS

  INTEGER FUNCTION StoI(String)
  ! Function to transform a String to an integer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                         :: String
    ! Local
    INTEGER                                              :: Lstring

    Lstring = LEN_TRIM(String)
    
    READ(String,'(I' // TRIM(ItoS(Lstring)) // ')')StoI

  END FUNCTION StoI

  REAL FUNCTION StoR(String)
  ! Function to transform a String to a real

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                         :: String

    ! Local
    INTEGER                                              :: Lstring

    Lstring = LEN_TRIM(String)

    READ(String,'(F' // TRIM(ItoS(Lstring)) // '.0)')StoR

  END FUNCTION StoR

  LOGICAL FUNCTION StoL(String)
  ! Function to transform a String to a boolean

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                         :: String

    IF (TRIM(String) == '.T.' .OR. TRIM(String) == '.true.' .OR. TRIM(String) == '.TRUE.'             &
      .OR. TRIM(String) == 'yes' .OR. TRIM(String) == 'YES' .OR. TRIM(String) == 'y' ) THEN
      StoL = .TRUE.
    ELSE
      StoL = .FALSE.
    END IF

  END FUNCTION StoL

  LOGICAL FUNCTION FileExist(filen)
  ! Checks if a file exists

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in)                         :: filen

    FileExist = .FALSE.
    INQUIRE(file=TRIM(filen), exist=FileExist)

    RETURN

  END FUNCTION FileExist

  INTEGER FUNCTION freeunit()
  ! provides the number of a free unit in which open a file

    IMPLICIT NONE

    LOGICAL                                              :: is_used

    is_used = .true.
    DO freeunit=10,100
      INQUIRE(unit=freeunit, opened=is_used)
      IF (.not. is_used) EXIT
    END DO

    RETURN

  END FUNCTION freeunit
  
  SUBROUTINE removeNONnum(String, newString)
  ! Subroutine to remove non numeric characters from a string

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                         :: String
    CHARACTER(LEN=*), INTENT(OUT)                        :: newString

    ! Local
    INTEGER                                              :: ic, inc, Lstring
    CHARACTER(len=S32)                                   :: fname

    !!!!!!! Variables
    ! String: string to remove non-numeric characters
    ! newString: resultant string
    fname = 'removeNONnum'

    Lstring = LEN_TRIM(String)
    newString = ''
    inc = 1
    DO ic=1, Lstring
      IF (ICHAR(String(ic:ic)) >= ICHAR('0') .AND. ICHAR(String(ic:ic)) <= ICHAR('9')) THEN
        newString(inc:inc) = String(ic:ic)
        inc = inc + 1
      END IF
    END DO

  END SUBROUTINE removeNONnum

!!!! netCDF

  SUBROUTINE open_NCfile(filename, ncid, idbg)
  ! Subroutine to open a basic netCDF file

    USE netcdf

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'
    
    CHARACTER(len=S128), INTENT(in)                      :: filename
    LOGICAL, OPTIONAL                                    :: idbg
    INTEGER, INTENT(out)                                 :: ncid
    
    ! Local
    INTEGER                                              :: rcode
    LOGICAL                                              :: fexist
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! filename: name of the file to create
    
    fname = 'open_NCfile'
    
    dbg = set_dbg(idbg)
    
    ncid = -1
    
    fexist = FileExist(filename)
    IF (.NOT.fexist) THEN
      CALL ErrMsg("netCDF file '" // TRIM(filename) // "' does not exist !!", fname, -1)
    END IF
    
    rcode = NF90_OPEN(filename, mode=NF90_NOWRITE, ncid=ncid)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)

    RETURN
    
  END SUBROUTINE open_NCfile

  SUBROUTINE get_Ndims_ncid(ncid, Ndims)
  ! Subroutine to get the amount of dimensions from an open netCDF file
  
    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ncid
    INTEGER, INTENT(out)                                 :: Ndims
    
    ! Local
    INTEGER                                              :: idd
    INTEGER                                              :: rcode
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! ncid: netCDF open object
    fname = 'get_Ndims_ncid'
    
    rcode = nf90_inquire(ncid, nDimensions=Ndims)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    
    RETURN

  END SUBROUTINE get_Ndims_ncid

  SUBROUTINE get_Nvars_ncid(ncid, Nvars)
  ! Subroutine to get the amount of variables from an open netCDF file
  
    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ncid
    INTEGER, INTENT(out)                                 :: Nvars
    
    ! Local
    INTEGER                                              :: idd
    INTEGER                                              :: rcode
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! ncid: netCDF open object
    fname = 'get_Nvars_ncid'
    
    rcode = nf90_inquire(ncid, nVariables=Nvars)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    
    RETURN

  END SUBROUTINE get_Nvars_ncid

  LOGICAL FUNCTION isin_ncunit(nid, varname)
  ! Function to tell if a given variable is inside a netcdf file unit

    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: nid
    CHARACTER(LEN=*), INTENT(in)                         :: varname

    ! Local
    INTEGER                                              :: vid, Ndims, Nvars
    INTEGER                                              :: iv, rcode
    CHARACTER(LEN=1000)                                  :: varinfile
    CHARACTER(len=S32)                                   :: fname

    !!!!!!! Variables
    ! nid: number of the open netCDF
    ! varname: name of the variable 

    fname = 'isin_ncunit'

    ! Initializing
    isin_ncunit = .FALSE.

    rcode = nf90_inquire(nid, Ndims, Nvars)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)

    DO iv=1, Nvars
      rcode = nf90_inquire_variable(nid, iv, name=varinfile)
      IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
      IF (TRIM(varinfile) == TRIM(varname)) THEN
        isin_ncunit = .TRUE.
        EXIT
      ELSE
        isin_ncunit = .FALSE.
      END IF
    END DO
    
    RETURN

  END FUNCTION isin_ncunit
  
  LOGICAL FUNCTION isin_dim_ncid(nid, dimname)
  ! Function to tell if a given dimension is inside a netcdf file unit

    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: nid
    CHARACTER(LEN=*), INTENT(in)                         :: dimname

    ! Local
    INTEGER                                              :: vid, Ndims, Nvars
    INTEGER                                              :: iv, id, rcode
    CHARACTER(len=S16)                                   :: diminfile
    CHARACTER(len=S32)                                   :: fname

    !!!!!!! Variables
    ! nid: number of the open netCDF
    ! dimname: name of the dimension

    fname = 'isin_dim_ncid'

    ! Initializing
    isin_dim_ncid = .FALSE.

    rcode = nf90_inquire(nid, Ndims, Nvars)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)

    DO id=1, Ndims
      rcode = nf90_inquire_dimension(nid, id, name=diminfile)
      IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
      IF (TRIM(diminfile) == TRIM(dimname)) THEN
        isin_dim_ncid = .TRUE.
        EXIT
      END IF
    END DO
    
    RETURN

  END FUNCTION isin_dim_ncid

  SUBROUTINE get_dimlength_ncid(ncid, dimn, dimlength)
  ! Subroutine to get the length of a dimension from an existing to an open netCDF
  
    USE netcdf

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER, INTENT(in)                                  :: ncid
    CHARACTER(len=*), INTENT(in)                         :: dimn
    INTEGER, INTENT(out)                                 :: dimlength
    
    ! Local
    INTEGER                                              :: idd
    INTEGER                                              :: Ndims
    INTEGER, DIMENSION(:), ALLOCATABLE                   :: dimids
    CHARACTER(len=S32), DIMENSION(:), ALLOCATABLE        :: dimns
    INTEGER                                              :: rcode
    LOGICAL                                              :: infile
    CHARACTER(len=S32)                                   :: fname
    
    !!!!!!! Variables
    ! ncid: netCDF open object
    ! dimn: name of dimension
    fname = 'get_dimlength_ncid'
    
    dimlength = -1
  
    infile = isin_dim_ncid(ncid, dimn)
    
    IF (.NOT.infile) THEN
      PRINT *, Errormsg
      CALL get_Ndims_ncid(ncid, Ndims)
      
      IF (ALLOCATED(dimns)) DEALLOCATE(dimns)
      ALLOCATE(dimns(Ndims))
      CALL get_dims_ncid(ncid, Ndims, dimns)
      
      CALL ErrMsgSAvail("ncid does not contain dimension named '" // TRIM(dimn) // "' !!", Ndims,     &
        dimns, fname)
    END IF
  
    rcode = nf90_inq_dimid(ncid, dimn, idd)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)

    rcode = nf90_inquire_dimension(ncid, idd, len=dimlength)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    
    RETURN

  END SUBROUTINE get_dimlength_ncid

  SUBROUTINE get_dims_ncid(ncid, Ndims, dimns, dimlengths)
  ! Subroutine to get the dimensions from an open netCDF file
  
    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ncid, Ndims
    CHARACTER(len=S32), DIMENSION(Ndims), INTENT(out)    :: dimns
    INTEGER, DIMENSION(Ndims), OPTIONAL                  :: dimlengths
    
    ! Local
    INTEGER                                              :: idd
    INTEGER                                              :: rcode    
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! ncid: netCDF open object
    ! Ndims: amount of fimensions in the file
    fname = 'get_dims_ncid'
    
    dimns = ''
    IF (PRESENT(dimlengths)) dimlengths = -1
    
    DO idd=1, Ndims
      IF (PRESENT(dimlengths)) THEN
        rcode = nf90_inquire_dimension(ncid, idd, name=dimns(idd), len=dimlengths(idd))
      ELSE
        rcode = nf90_inquire_dimension(ncid, idd, name=dimns(idd))
      END IF
      IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
      IF (PRESENT(dimlengths)) THEN
        IF (dimlengths(idd) == NF90_UNLIMITED) PRINT *,'  ' // TRIM(fname) // ": dimension '" //      &
          TRIM(dimns(idd)) // "' is unlimited ", NF90_UNLIMITED
      END IF
    END DO
    
    RETURN

  END SUBROUTINE get_dims_ncid
  
  SUBROUTINE get_varS0D_ncid(ncid, StrL, vname, vals, idbg)
  ! Subroutine to get a 0D String variable from a netCDF file unit passing starting and ending
  ! First dimension accounts for length of the string values
  
    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ncid, StrL
    CHARACTER(LEN=S32), INTENT(in)                       :: vname
    LOGICAL, OPTIONAL                                    :: idbg
    CHARACTER(len=1), DIMENSION(StrL), INTENT(out)       :: vals
  
    ! Local
    INTEGER                                              :: i, Ndims, Nvars
    INTEGER                                              :: rcode, varid, vartype, rank
    LOGICAL                                              :: vfound
    INTEGER, DIMENSION(:), ALLOCATABLE                   :: invarids
    CHARACTER(len=S32), DIMENSION(:), ALLOCATABLE        :: invarns
    LOGICAL                                              :: infile
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname

    !!!!!!! Variables
    ! ncid: netCDF file identifier
    ! StrL: length of the String values
    ! vals: values to get
    ! vname: name of the variable to get
    fname = 'get_varS0D_ncid'
        
    dbg = set_dbg(idbg)
    
    ! Checking existence of variable
    infile = isin_ncunit(ncid, vname)
    IF (.NOT.infile) THEN
       rcode = nf90_inquire(ncid, Ndims, Nvars)
       IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
       IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
       ALLOCATE(invarids(Nvars))
       IF (ALLOCATED(invarns)) DEALLOCATE(invarns)
       ALLOCATE(invarns(Nvars))
       
       CALL get_vars_ncid(ncid, Nvars, invarids, invarns)
       PRINT *,'  Available variales in netcdf unit: ', ("'" // TRIM(invarns(i)) // "', ", i=1, Nvars)
       CALL ErrMsg("Passed variable '" // TRIM(vname) // "' not found", fname, -1)
    END IF

    rcode = nf90_inq_varid(ncid, vname, varid)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    
    ! Checking variable
    rcode = nf90_inquire_variable(ncid, varid, xtype=vartype, ndims=rank)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    ! First dimension accounts for length of the string values
    rank = rank - 1 
    IF (vartype /= NF90_CHAR) CALL ErrMSg("Wrong type of variable '" // TRIM(vname) //                &
      "'. I expected NF90_CHAR= " // TRIM(ItoS(NF90_DOUBLE)) // ' instead got ' //                    &
      TRIM(ItoS(vartype)), fname, -1)
    IF (rank /= 0) CALL ErrMSg("Wrong rank of variable '" //TRIM(vname)// "'. I expected rank= 0 " // &
      ' instead got ' // TRIM(ItoS(rank)), fname, -1)

    rcode = nf90_get_var(ncid, varid, vals)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
  
    IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
    IF (ALLOCATED(invarns)) DEALLOCATE(invarns)

    RETURN

  END SUBROUTINE get_varS0D_ncid

  SUBROUTINE get_varR0D_ncid(ncid, vname, vals, idbg)
  ! Subroutine to get a 0D Real variable from a netCDF file unit passing starting and ending

    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ncid
    CHARACTER(LEN=S32), INTENT(in)                       :: vname
    LOGICAL, OPTIONAL                                    :: idbg
    REAL, INTENT(out)                                    :: vals
  
    ! Local
    INTEGER                                              :: i, Ndims, Nvars
    INTEGER                                              :: rcode, varid, vartype, rank
    LOGICAL                                              :: vfound
    INTEGER, DIMENSION(:), ALLOCATABLE                   :: invarids
    CHARACTER(len=S32), DIMENSION(:), ALLOCATABLE        :: invarns
    LOGICAL                                              :: infile
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname

    !!!!!!! Variables
    ! ncid: netCDF file identifier
    ! vals: values to get
    ! vname: name of the variable to get
    fname = 'get_varR0D_ncid'
        
    dbg = set_dbg(idbg)
    
    ! Checking existence of variable
    infile = isin_ncunit(ncid, vname)
    IF (.NOT.infile) THEN
       rcode = nf90_inquire(ncid, Ndims, Nvars)
       IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
       IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
       ALLOCATE(invarids(Nvars))
       IF (ALLOCATED(invarns)) DEALLOCATE(invarns)
       ALLOCATE(invarns(Nvars))
       
       CALL get_vars_ncid(ncid, Nvars, invarids, invarns)
       PRINT *,'  Available variales in netcdf unit: ', ("'" // TRIM(invarns(i)) // "', ", i=1, Nvars)
       CALL ErrMsg("Passed variable '" // TRIM(vname) // "' not found", fname, -1)
    END IF

    rcode = nf90_inq_varid(ncid, vname, varid)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    
    ! Checking variable
    rcode = nf90_inquire_variable(ncid, varid, xtype=vartype, ndims=rank)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    IF (vartype /= NF90_FLOAT) CALL ErrMSg("Wrong type of variable '" // TRIM(vname) //               &
      "'. I expected NF90_FLOAT= " // TRIM(ItoS(NF90_FLOAT)) // ' instead got ' //                    &
      TRIM(ItoS(vartype)), fname, -1)
    IF (rank /= 0) CALL ErrMSg("Wrong rank of variable '" //TRIM(vname)// "'. I expected rank= 0 " // &
      ' instead got ' // TRIM(ItoS(rank)), fname, -1)

    rcode = nf90_get_var(ncid, varid, vals)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
  
    IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
    IF (ALLOCATED(invarns)) DEALLOCATE(invarns)

    RETURN

  END SUBROUTINE get_varR0D_ncid

  SUBROUTINE get_varR1D_ncid(ncid, id1, ed1, vname, vals, idbg)
  ! Subroutine to get a 1D Real variable from a netCDF file unit passing starting and ending

    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ncid, id1, ed1
    CHARACTER(LEN=S32), INTENT(in)                       :: vname
    LOGICAL, OPTIONAL                                    :: idbg
    REAL, DIMENSION(id1:ed1),INTENT(out)                 :: vals
  
    ! Local
    INTEGER                                              :: i, Ndims, Nvars
    INTEGER                                              :: rcode, varid, vartype, rank
    INTEGER                                              :: d1
    LOGICAL                                              :: vfound
    INTEGER, DIMENSION(:), ALLOCATABLE                   :: invarids
    CHARACTER(len=S32), DIMENSION(:), ALLOCATABLE        :: invarns
    LOGICAL                                              :: infile
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname

    !!!!!!! Variables
    ! ncid: netCDF file identifier
    ! id1,ed1: shape of the matrix
    ! vals: values to get
    ! vname: name of the variable to get
    fname = 'get_varR1D_ncid'
        
    dbg = set_dbg(idbg)
    
    ! Checking existence of variable
    infile = isin_ncunit(ncid, vname)
    IF (.NOT.infile) THEN
       rcode = nf90_inquire(ncid, Ndims, Nvars)
       IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
       IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
       ALLOCATE(invarids(Nvars))
       IF (ALLOCATED(invarns)) DEALLOCATE(invarns)
       ALLOCATE(invarns(Nvars))
       
       CALL get_vars_ncid(ncid, Nvars, invarids, invarns)
       PRINT *,'  Available variales in netcdf unit: ', ("'" // TRIM(invarns(i)) // "', ", i=1, Nvars)
       CALL ErrMsg("Passed variable '" // TRIM(vname) // "' not found", fname, -1)
    END IF

    rcode = nf90_inq_varid(ncid, vname, varid)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    
    ! Checking variable
    rcode = nf90_inquire_variable(ncid, varid, xtype=vartype, ndims=rank)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    IF (vartype /= NF90_FLOAT) CALL ErrMSg("Wrong type of variable '" // TRIM(vname) //               &
      "'. I expected NF90_FLOAT= " // TRIM(ItoS(NF90_FLOAT)) // ' instead got ' //                    &
      TRIM(ItoS(vartype)), fname, -1)
    IF (rank /= 1) CALL ErrMSg("Wrong rank of variable '" //TRIM(vname)// "'. I expected rank= 1 " // &
      ' instead got ' // TRIM(ItoS(rank)), fname, -1)

    d1 = ed1-id1+1
    rcode = nf90_get_var(ncid, varid, vals, start=(/id1/), count=(/d1/))
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
  
    IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
    IF (ALLOCATED(invarns)) DEALLOCATE(invarns)

    RETURN

  END SUBROUTINE get_varR1D_ncid

  SUBROUTINE get_varD3D_ncid(ncid, id1, id2, id3, ed1, ed2, ed3, vname, vals, idbg)
  ! Subroutine to get a 3D double Real variable from a netCDF file unit passing starting and ending

    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: ncid, id1, ed1, id2, ed2, id3, ed3
    CHARACTER(LEN=S32), INTENT(in)                       :: vname
    LOGICAL, OPTIONAL                                    :: idbg
    REAL(kind=rr8), DIMENSION(id1:ed1,id2:ed2,id3:ed3),                                               &
      INTENT(out)                                        :: vals
  
    ! Local
    INTEGER                                              :: i, Ndims, Nvars
    INTEGER                                              :: rcode, varid, vartype, rank
    INTEGER                                              :: d1, d2, d3
    LOGICAL                                              :: vfound
    INTEGER, DIMENSION(:), ALLOCATABLE                   :: invarids
    CHARACTER(len=S32), DIMENSION(:), ALLOCATABLE        :: invarns
    LOGICAL                                              :: infile
    LOGICAL                                              :: dbg
    CHARACTER(len=S32)                                   :: fname

    !!!!!!! Variables
    ! ncid: netCDF file identifier
    ! id1,ed1,id2,ed2,id3,ed3: shape of the matrix
    ! vals: values to get
    ! vname: name of the variable to get
    fname = 'get_varD3D_ncid'
        
    dbg = set_dbg(idbg)
    
    ! Checking existence of variable
    infile = isin_ncunit(ncid, vname)
    IF (.NOT.infile) THEN
       rcode = nf90_inquire(ncid, Ndims, Nvars)
       IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
       IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
       ALLOCATE(invarids(Nvars))
       IF (ALLOCATED(invarns)) DEALLOCATE(invarns)
       ALLOCATE(invarns(Nvars))
       
       CALL get_vars_ncid(ncid, Nvars, invarids, invarns)
       PRINT *,'  Available variales in netcdf unit: ', ("'" // TRIM(invarns(i)) // "', ", i=1, Nvars)
       CALL ErrMsg("Passed variable '" // TRIM(vname) // "' not found", fname, -1)
    END IF

    rcode = nf90_inq_varid(ncid, vname, varid)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    
    ! Checking variable
    rcode = nf90_inquire_variable(ncid, varid, xtype=vartype, ndims=rank)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    IF (vartype /= NF90_DOUBLE) CALL ErrMSg("Wrong type of variable '" // TRIM(vname) //              &
      "'. I expected NF90_DOUBLE= " // TRIM(ItoS(NF90_DOUBLE)) // ' instead got ' //                  &
      TRIM(ItoS(vartype)), fname, -1)
    IF (rank /= 3) CALL ErrMSg("Wrong rank of variable '" //TRIM(vname)// "'. I expected rank= 3 " // &
      ' instead got ' // TRIM(ItoS(rank)), fname, -1)

    d1 = ed1-id1+1
    d2 = ed2-id2+1
    d3 = ed3-id3+1
    rcode = nf90_get_var(ncid, varid, vals, start=(/id1,id2,id3/), count=(/d1,d2,d3/))
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
  
    IF (ALLOCATED(invarids)) DEALLOCATE(invarids)
    IF (ALLOCATED(invarns)) DEALLOCATE(invarns)

    RETURN

  END SUBROUTINE get_varD3D_ncid
  
  ! From UNIDATA: https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  SUBROUTINE handle_errf(st, fn)
  ! Subroutine to provide the error message when something with netCDF went wrong (including fname)

    USE netcdf

    INTEGER, INTENT(in)                                  :: st
    CHARACTER(len=*), INTENT(in)                         :: fn
    CHARACTER (LEN=S2048)                                :: diag_message

    !!!!!!! Variables
    ! st: netCDF status number
    ! fn: function name from which it is used

    IF (st /= nf90_noerr) THEN
      WRITE (diag_message , * ) &
      TRIM(errormsg) // ":   " // TRIM(fn) // ": " // TRIM(nf90_strerror(st))
      CALL wrf_error_fatal ( diag_message )

    END IF

  END SUBROUTINE handle_errf

  SUBROUTINE close_NCfile(filename, ncid, idbg)
  ! Subroutine to close a netCDF file

    USE netcdf

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'
    
    CHARACTER(len=S128), INTENT(in)                      :: filename
    LOGICAL, OPTIONAL                                    :: idbg
    INTEGER, INTENT(out)                                 :: ncid
    
    ! Local
    INTEGER                                              :: rcode
    LOGICAL                                              :: fexist
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! filename: name of the file to create
    
    fname = 'close_NCfile'
    
    dbg = set_dbg(idbg)
      
    rcode = NF90_CLOSE(ncid)
    IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)

    RETURN
    
  END SUBROUTINE close_NCfile  
    
  SUBROUTINE get_vars_ncid(ncid, Nvars, invarids, invarns)
  ! Subroutine to provide all the ids and names of the variables existing to an open netCDF
  !    previously knowing how many variables the file has via Nvars 
  !      rcode = nf90_inq_nvars(ncid, Nvars)
  
    USE netcdf

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER, INTENT(in)                                  :: ncid, Nvars
    INTEGER, DIMENSION(Nvars), INTENT(out)               :: invarids
    CHARACTER(len=S32), DIMENSION(Nvars), INTENT(out)    :: invarns
    
    ! Local
    INTEGER                                              :: iv
    INTEGER                                              :: rcode
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! ncid: netCDF open object
    ! Nvars: amount of variables
    ! invarids: Ids of the variables
    ! invarns: name of the variables
    fname = 'get_vars_ncid'
    
    invarids = -1
    invarns = ''
        
    DO iv=1, Nvars
      invarids(iv) = iv
      rcode = nf90_inquire_variable(ncid, invarids(iv), invarns(iv))
      IF (rcode /= NF90_NOERR) CALL handle_errf(rcode, fname)
    END DO
    
    RETURN

  END SUBROUTINE get_vars_ncid


#endif
END MODULE module_simradar_tools
