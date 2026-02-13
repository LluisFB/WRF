MODULE module_simradar
#ifdef SIMRADAR
! Module for the simulation of surface polimeric radar following simrad.AR, Galligani et al (CONGREMET XIII), 2018

  CONTAINS

  SUBROUTINE simradar_calc(                                                                           &
                      is_restart, is_output_time, nsteps,                                             & 
                      dt, xtime, curr_secs2,                                                          &
                      timedbg,                                                                        &
                      londbg, latdbg,                                                                 &
                      idbg, jdbg,                                                                     &

                      ! Input
                      lon, lat,                                                                       &
                      qv, qc, qr, qs, qi, qg, qh,                                                     &
                      mu, mub, th300,                                                                 &
                      ph, phb, rho,                                                                   &
                      nradar, nrnn, nrnncoord, nrcoord, nr4coord,                                     &
                      e_angle, e_azimuth, e_dist,                                                     &
                      ngrid_radar, ngrid3_radar,                                                      &
                      gatename, gatelon, gatelat, gatehgt, gateband                                   &
                      radarlons, radarlats, radardist, radarhgts, radarlhgts,                         &
                      azimuths,                                                                       &
                      ijdbg,                                                                          &
                      a, f, elipsn,                                                                   &
                      gatehdist, inradar, radnum,                                                     &
                      radarlons, radarlats, radardist, radarhgts, radarlhgts,                         &
                      radarlons_bnds, radarlat_bnds, radarhgt_bnds,                                   &
                      ngridradar, ijgridradar, ijdistgridradar,                                       &
                      lutgrid,                                                                        &
                      missing_value,                                                                  &

                      ! Output
                      zh_rain, zdr_rain, kdp_rain,                                                    &
                      zh_snow, zdr_snow, kdp_snow,                                                    &
                      zh_grau, zdr_grau, kdp_grau,                                                    &
                      zh_ice, zdr_ice, kdp_ice,                                                       &
                      zh_hail, zdr_hail, kdp_hail,                                                    &
                      zh_total, zdr_total, kdp_total,                                                 &
                      tv_rain, tv_snow, tv_grau,                                                      &
                      tv_ice, tv_hail,                                                                &
                      tv_wgt, vr,                                                                     &

                      ids,ide, jds,jde, kds,kde,                                                      &
                      ims,ime, jms,jme, kms,kme,                                                      &
                      ips,ipe, jps,jpe, kps,kpe,                                                      & ! patch  dims
                      i_start,i_end,j_start,j_end,kts,kte,num_tiles                                   &
                    )

!----------------------------------------------------------------------
! Main subroutine to compute simrad.AR variables

    USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
    USE module_configure 
    USE module_simradar_tools
    USE module_diagvar_simradar

    IMPLICIT NONE

!======================================================================
! Definitions
!-----------
!-- DT            time step (second)
!-- XTIME         forecast time
!-- curr_secs2    current time in seconds since simulation restart
!-- radt          frequency of radiative calls (minutes)
!
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- ips           start index for i in patch
!-- ipe           end index for i in patch
!-- jps           start index for j in patch
!-- jpe           end index for j in patch
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!-- num_soil_layers number of soil layers
!-- n_moist     number of water species
!-- m_qv, m_qc, m_qr, m_qs, m_qi, m_qg, m_qh: moist indices for each specie: vapour, cloud, rain, snow, ice, graupel, ice and hail
!
!-----------
! is_restart: whether if simulation is a restart
! is_output_time: whether if this iteration is an output time-step

! ------- INPUT vars -------
! lon: longitude
! lat: latitude
! qv: water vapour mixing ratio [kgkg-1]
! qc: cloud mixing ratio [kgkg-1]
! qr: rain mixing ratio [kgkg-1]
! qs: snow mixing ratio [kgkg-1]
! qi: ice mixing ratio [kgkg-1]
! qg: graupel mixing ratio [kgkg-1]
! qh: hail mixing ratio [kgkg-1]
! th300: potential air temperature -300. [K]
! mu: perturbation dry air mass in column [Pa]
! mub: base dry air mass in column [Pa]
! sina: local sine of map rotation [1]
! cosa: local cosine of map rotation [1]
! ph: 3D gepotential perturbation [m2s-2]
! phb: 3D gepotential base state [m2s-2]
! rho: air density [kgm-3]

!! Radar variables
! nradar: amount of radars
! gatename: name of the radar
! gatelon: longitude of the radar location
! gatelat: latitude of the radar location
! gatehgt: height of the radar location [m]
! distmax: horizontal range of measuremnt of the radar [km]
! heightmax: vertical range of measuremnt of the radar [km]
! radards: radial resolution of the radar space [m]
! radii_ddist: radii (km) at which arc resolution is the same as minimum ddist (same for all radars)
! em_ray: type of propagation of em ray
!    1: 4/3 Ref 
! eloption: type of ellipsoide to use
!    1: 'GRS80 / WGS84  (NAD83)'
!        a = 6378137.d0
!        f = 1.d0/298.25722210088d0
!    2: 'Clarke 1866    (NAD27)'
!        a = 6378206.4d0
!        f = 1.d0/294.9786982138d0
! e_angle: amount of angle beams (horizontal) in the radar space [computed a 2*pi*radii_ddist / (min(ddist))
! e_azimuth: amount of azimuth angles (vertical) in the radar space (same for all radars)
! azimuths: azimuth angles
!   IF azimuths(1) == -1; aximuths are internally computed from Nazimuth, by regular intervals 
!     30. / (e_azimuth-1)
! e_dist: amount of distances (horizontal) in the radar space [computed as 2.*max(distmax) / min(ddist)]
! lutfilen: name of the netCDF file with the look-up-tables for a given miscrophysics of the 
!   polymetries obntained via the T-matrix derived
! doppler_scheme: mtehod to use to retrieve the radial velocities after LUT
!   0: Reflectivity weighted terminal velocity
!   1:  Get PSD integrated backscattering cross section int = v S/S [not implemented]
! ngrid_radar: maximum amount of radar cells in a WRF grid
!! Compressed dimensions
! nrnn = nradar*ngrid_radar
! nrcoord = nradar*2
! nrnncoord = nradar*ngrid_radar*2
! nncoord = ngrid_radar*2
! nr4coord = nradar*4

! timedbg: frecuency in minutes of time-debugging
! londbg: longitude in degrees of grid-point debugging (-999., not used)
! latdbg: latitude in degrees of grid-point debugging (-999., not used)

! Output
!-------
! zh_rain: Reflectivity for Rain [dBZ]
! zdr_rain: ZDR for Rain [dBz]
! kdp_rain: KDP for Rain [degree/km]
! zh_snow: Reflectivity for Snow [dBZ]
! zdr_snow: ZDR for Snow [dBz]
! kdp_snow: KDP for Snow [degree/km]
! zh_grau: Reflectivity for Graupel [dBZ]
! zdr_grau: ZDR for Graupel [dBz]
! kdp_grau: KDP for Graupel [degree/km]
! zh_ice: Reflectivity for cloud Ice [dBZ]
! zdr_ice: ZDR for cloud Ice [dBz]
! kdp_ice: KDP for cloud Ice [degree/km]
! zh_hail: Reflectivity for Hail [dBZ]
! zdr_hail: ZDR for Hail [dBz]
! kdp_hail: KDP for Hail [degree/km]
! zh_total: Reflectivity for all species [dBZ]
! zdr_total: ZDR for all species [dBz]
! kdp_total: KDP for all species [degree/km]
! tv_rain: terminal velocity for Rain [ms-1]
! tv_snow: terminal velocity for Snow [ms-1]
! tv_grau: terminal velocity for Graupel [ms-1]
! tv_ice: terminal velocity for cloud Ice [ms-1]
! tv_hail: terminal velocity for Hail [ms-1]
! tv_wgt: weighted terminal velocity [ms-1]
! vr: Doppler velocity [ms-1]

!======================================================================

    INTEGER,      INTENT(in)                             ::                                           &
                                            ids,ide, jds,jde, kds,kde,                                &
                                            ims,ime, jms,jme, kms,kme,                                &
                                            ips,ipe, jps,jpe, kps,kpe,                                &
                                            kts,kte,                                                  &
                                            num_tiles,                                                &
                                            n_moist


    INTEGER, DIMENSION(num_tiles), INTENT(in)            :: i_start, i_end, j_start, j_end
    INTEGER, INTENT(in)                                  :: m_qv, m_qc, m_qr, m_qs, m_qi, m_qg, m_qh

    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: lon, lat
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: mu, mub, sina, cosa
    REAL, DIMENSION(kms:kme), INTENT(in)                 :: deta

    REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(in) :: qv, qc, qr, qs, qi, qg, qh, ph, phb, rho

    REAL, INTENT(in)                                     :: DT, XTIME
    REAL, INTENT(in)                                     :: curr_secs2
    LOGICAL, INTENT(in)                                  :: is_restart, is_output_time
    INTEGER, INTENT(in)                                  :: timedbg
    REAL, INTENT(in)                                     :: londbg, latdbg
    INTEGER, INTENT(in)                                  :: idbg, jdbg

    INTEGER, INTENT(inout)                               :: nsteps
  
! Specific Radar variables
    INTEGER, INTENT(in)                                  :: nradar, ngrid_radar
    INTEGER, INTENT(in)                                  :: e_azimuth, eloption
    INTEGER, INTENT(in)                                  :: e_angle, e_dist
    REAL, DIMENSION(nradar), INTENT(in)                  :: gatelon, gatelat, gatehgt, distmax,       &
      heightmax, radards
!    INTEGER, DIMENSION(nradar), INTENT(in)               :: gatetype, gateband, radinterp, em_ray
    INTEGER, DIMENSION(nradar), INTENT(in)               :: em_ray
    REAL, INTENT(in)                                     :: radii_ddist
    REAL, INTENT(in)                                     :: missing_value
    REAL, DIMENSION(e_azimuth), INTENT(inout)            :: azimuths
    REAL(kind=rr8)                                       :: a, f
    CHARACTER(len=30)                                    :: elips
    ! Compressed dimensions
    INTEGER, INTENT(in)                                  :: nrnn, nrnncoord, nrcoord, nr4coord 
    INTEGER, INTENT(in)                                  :: ngrid_radar, ngrid3_radar
    INTEGER, DIMENSION(2), INTENT(in)                    :: ijdbg
    REAL, DIMENSION(ims:ime, jms:jme, nradar),                                                        &
      INTENT(int)                                        :: gatehdist
    INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(in)     :: inradar
    INTEGER, DIMENSION(ims:ime, jms:jme, nradar),                                                     &
      INTENT(in)                                         :: radnum
    REAL, DIMENSION(e_angle, e_dist, nradar),INTENT(in)  :: radarlons, radarlats, radardist,          &
    REAL, DIMENSION(e_azimuth,e_dist,nradar),INTENT(in)  :: radarhgts, radarlhgts
    REAL,DIMENSION(e_angle,e_dist,nr4coord), INTENT(in)  :: radarlons_bnds, radarlat_bnds,            &
    REAL, DIMENSION(e_azimuth,e_dist,nr4coord), INTENT(in)                                            &
                                                         :: radarhgt_bnds
    INTEGER, DIMENSION(ims:ime,jms:jme, nradar),                                                      &
      INTENT(in)                                         :: ngridradar
    INTEGER, DIMENSION(ims:ime, jms:jme, nrnncoord),                                                  &
      INTENT(in)                                         :: ijgridradar
    REAL, DIMENSION(ims:ime, jms:jme, nrnn), INTENT(in)  :: ijdistgridradar
    
    ! LUT space
    TYPE(lut_vars), INTENT(in)                           :: lutgrid

    ! Output
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: zh_rain, zdr_rain, kdp_rain
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: zh_snow, zdr_snow, kdp_snow
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: zh_grau, zdr_grau, kdp_grau
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: zh_ice, zdr_ice, kdp_ice
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: zh_hail, zdr_hail, kdp_hail
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: zh_total, zdr_total, kdp_total
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: tv_rain, tv_snow, tv_grau
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: tv_ice, tv_hail
    REAL, DIMENSION(e_angle, e_azimuth, nrd)             :: tv_wgt, vr
    
    ! LOCAL  VAR
    INTEGER                                              :: i, j, k, m, its, ite, jts, jte, ij
    INTEGER                                              :: idp, jdp, dimz, dimz1, dz1
    INTEGER                                              :: i1, i_1, j1, j_1, k2, tdbg
    INTEGER, DIMENSION(num_tiles)                        :: i2, j2
    INTEGER                                              :: inrnn, inrnncoord, inncoord, inrcoord,    &
      inr4coord

    REAL                                                 :: xtimep
    REAL                                                 :: prv
    
    LOGICAL, EXTERNAL                                    :: wrf_dm_on_monitor
    CHARACTER(len=256)                                   :: msg
    LOGICAL                                              :: ijdbg, cdone
    INTEGER                                              :: dbg_level
 
    INTEGER                                              :: ir, ia, ih, id, iband
    INTEGER, DIMENSION(ngrid_radar, coord2)              :: ijgridradar_NOc
    REAL, DIMENSION(ngrid_radar, nradar)                 :: ijdistgridradar_NOc

    ! radar 3-dimensional
    INTEGER                                              :: ijkNtot
    INTEGER, DIMENSION(kms:kme)                          :: ijkngridradar
    INTEGER, DIMENSION(kms:kme, 3, ngrid3_radar)         :: ijkgridradar
    INTEGER, DIMENSION(kms:kme, 27, ngrid3_radar)        :: wghts

    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: zh_rain_NOc, zdr_rain_NOc, kdp_rain_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: zh_snow_NOc, zdr_snow_NOc, kdp_snow_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: zh_grau_NOc, zdr_grau_NOc, kdp_grau_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: zh_ice_NOc, zdr_ice_NOc, kdp_ice_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: zh_hail_NOc, zdr_hail_NOc, kdp_hail_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: zh_total_NOc, zdr_total_NOc, kdp_total_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: tv_rain_NOc, tv_snow_NOc, tv_grau_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: tv_ice_NOc, tv_hail_NOc
    REAL, DIMENSION(e_angle, e_azimuth, e_dist)          :: tv_wgt_NOc, vr_NOc

    fname = 'simradar_calc'

    ! Showing version of the module
    IF (curr_secs2 == 0. ) THEN
      WRITE(msg,*) '  simrad.AR module version: ' // TRIM(simradARversion) // ' ...'
      PRINT *, TRIM(msg)
      CALL wrf_debug(0,msg)
    END IF

    ! timedbg
    !  when time-debug in minutes
    tdbg = MOD(xtime, timedbg*1.)

    ! Getting the overall debug level
    CALL get_wrf_debug_level( dbg_level ) 
    cdone = .FALSE.

    ! L. Fita: NOTE about efficiency
    !   I am not sure if this is the more efficient way:
    !     1.- IF --> 2.- repeat DO
    !     1.- DO -->  2.- IF

    dimz = kme - kms + 1
    dimz1 = dimz + 1
    dz1 = dimz - 1
    mutot = mu + mub
    zg = (ph + phb)/9.81

    ! For test purposes
    IF ((londbg == -999.) .OR. (latdbg == -999.)) THEN
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , num_tiles
        i2(ij) = INT(0.5*(i_start(ij) + i_end(ij)))
        j2(ij) = INT(0.5*(j_start(ij) + j_end(ij)))
      END DO
      !$OMP END PARALLEL DO
    ELSE
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , num_tiles
        i2(ij) = idbg
        j2(ij) = jdbg
      END DO
      !$OMP END PARALLEL DO
    END IF
    k2 = 0.5*(kms+kme)

    ! Computing variables 
    ijdbg = .FALSE.

    !! Looping for all radars
    DO ir=1, nradar   
    
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , num_tiles
        DO i = i_start(ij), i_end(ij)
          DO j = j_start(ij), j_end(ij)
            IF (ij == 1 .AND. i == i2(ij) .AND. j == j2(ij)) ijdbg = .TRUE.
            i1 = MIN(i+1,i_end(ij))
            i_1 = MAX(i-1,i_start(ij))
            j1 = MIN(j+1,j_end(ij))
            j_1 = MAX(j-1,j_start(ij))
          
            IF (ij21Ncoords(i,j) == 0) CYCLE
        
            ii1 = i1-1
            ei1 = i1+1
            ij1 = j1-1
            ej1 = j1+1
            IF (i == ims) ii1 = i
            IF (i == ime) ei1 = i
            IF (j == jms) ij1 = j
            IF (j == jme) ej1 = j

            ! Decompressing
            DO n=1, ngridradar(i,j,ir)
              inrnn = (ir-1)*ngrid_radar + n
              ijdistgridradar_NOc(n) = ijdistgridradar(i,j,inrnn)
              DO l=1, 2
                inrnncoord = (ir-1)*ngrid_radar*2 + (n-1)*2 + l
                ijgridradar_NOc(n,l) = ijgridradar(i,j,inrnncoord)
              END DO
            END DO

            ! Computing 27 weights around each WRF i,j,k for all radar cells
            CALL weights_hcoordsR27around_polar_ellipsoide_1D(i, j, ii1, ei1, ij1, ej1, kms, kme, a,    &
              f, elipsn, lon(ii1:ei1,ij1:ej1), lat(ii1:ei1,ij1:ej1), zg(ii1:ei1,kms:kme,ij1:ej1),       &
              missing_value, ngridradar(i,j,ir), ngridradar(i,j,ir), 2, ijgridradar_NOc, 9, 1, e_dist,  &
              1, e_angle, 1, e_azimuth, radarlons(:,:,ir), radarlats(:,:,ir), radarhgts(:,:,ir),        &
              missing_value, ngrid3_radar, ijkngridradar, 3, ijkgridradar, ijkNtot, 27, wghts, dbg)
           
            loop_z: DO k=kms, kme
            IF (ijkngridradar(k) == 0) CYCLE

            printdbg = .FALSE.
            IF (dbg) THEN
              IF (MOD(i+j+k,10) == 0 .AND. ijkfound < 10) THEN
                IF (ijkfound < 10) THEN
                  printdbg = .TRUE.
                 ijkfound = ijkfound + 1
                END IF
              END IF
            END IF
  
            IF (printdbg) PRINT *,'  ' // TRIM(fname) // ": ", ijkfound, ' i,j,k ', i,j,k,            &
              ' amount of radar cells:', ijkngridradar(k)

            ! 27 values around reference space at i,j,k
            i123 = 1
            cvals = 0.
            IF (k==kms) THEN
              DO ii1=-1,1
                DO ii2=-1,1
                  DO ii3=0,1
                    cvals(i123,1) = qv(i+ii1,k+ii3,j+ii2)
                    cvals(i123,2) = qc(i+ii1,k+ii3,j+ii2)
                    cvals(i123,3) = qr(i+ii1,k+ii3,j+ii2)
                    cvals(i123,4) = qi(i+ii1,k+ii3,j+ii2)
                    cvals(i123,5) = qs(i+ii1,k+ii3,j+ii2)
                    cvals(i123,6) = qg(i+ii1,k+ii3,j+ii2)
                    cvals(i123,7) = ua(i+ii1,k+ii3,j+ii2)
                    cvals(i123,8) = va(i+ii1,k+ii3,j+ii2)
                    cvals(i123,9) = wa(i+ii1,k+ii3,j+ii2)
                    cvals(i123,10) = ta(i+ii1,k+ii3,j+ii2)
                    cvals(i123,11) = press(i+ii1,k+ii3,j+ii2)
                    i123 = i123 + 1
                  END DO
                END DO
              END DO
            ELSE IF (k1==kme) THEN
              DO ii1=-1,1
                DO ii2=-1,1
                  DO ii3=-1,0
                    cvals(i123,1) = qv(i+ii1,k+ii3,j+ii2)
                    cvals(i123,2) = qc(i+ii1,k+ii3,j+ii2)
                    cvals(i123,3) = qr(i+ii1,k+ii3,j+ii2)
                    cvals(i123,4) = qi(i+ii1,k+ii3,j+ii2)
                    cvals(i123,5) = qs(i+ii1,k+ii3,j+ii2)
                    cvals(i123,6) = qg(i+ii1,k+ii3,j+ii2)
                    cvals(i123,7) = ua(i+ii1,k+ii3,j+ii2)
                    cvals(i123,8) = va(i+ii1,k+ii3,j+ii2)
                    cvals(i123,9) = wa(i+ii1,k+ii3,j+ii2)
                    cvals(i123,10) = ta(i+ii1,k+ii3,j+ii2)
                    cvals(i123,11) = press(i+ii1,k+ii3,j+ii2)
                    i123 = i123 + 1
                  END DO
                END DO    
              END DO
            ELSE
              DO ii1=-1,1
                DO ii2=-1,1
                  DO ii3=-1,1
                    cvals(i123,1) = qv(i+ii1,k+ii3,j+ii2)
                    cvals(i123,2) = qc(i+ii1,k+ii3,j+ii2)
                    cvals(i123,3) = qr(i+ii1,k+ii3,j+ii2)
                    cvals(i123,4) = qi(i+ii1,k+ii3,j+ii2)
                    cvals(i123,5) = qs(i+ii1,k+ii3,j+ii2)
                    cvals(i123,6) = qg(i+ii1,k+ii3,j+ii2)
                    cvals(i123,7) = ua(i+ii1,k+ii3,j+ii2)
                    cvals(i123,8) = va(i+ii1,k+ii3,j+ii2)
                    cvals(i123,9) = wa(i+ii1,k+ii3,j+ii2)
                    cvals(i123,10) = ta(i+ii1,k+ii3,j+ii2)
                    cvals(i123,11) = press(i+ii1,k+ii3,j+ii2)
                    i123 = i123 + 1
                  END DO
                END DO
              END DO
            END IF
          
            ! Removing too small values...
            DO ii1=1, 27
              DO ii2=1, 11
                IF (cvals(ii1,ii2) < eps) cvals(ii1,ii2) = 0. 
              END DO
            END DO
  
            loop_ptsijk: DO ic= 1, ijkngridradar(k)
              !! Angle
              iy = ijkgridradar(k,ic,1)
              !! Dist
              ix = ijkgridradar(k,ic,2)
              !! Azimuth
              iz = ijkgridradar(k,ic,3)
              
              IF (printdbg) PRINT *,'    cell: ', ic, ' coords; ', ijkgridradar(k,ic,:), ' || lon: ', &
                radarlons(iy,ix,ir), ' lat: ', radarlats(iy,ix,ir), ' hgt: ', radarhgts(iz,ix,ir),    &
                ' amount of values: ', i123
              ! Need to filter Nan and missing values ?
              !ijkvals(i,j,k,ic) = 0.
              !DO ijk=1,acoords3
              !  ijkvals(i,j,k,ic) = ijkvals(i,j,k,ic) + cvals(ijk)*wghts(i,j,k,ic,ijk)
              !END DO
              iqv(iy,iz,ix) = SUM(cvals(:,1)*wghts(k,ic,:))
              iqc(iy,iz,ix) = SUM(cvals(:,2)*wghts(k,ic,:))
              iqr(iy,iz,ix) = SUM(cvals(:,3)*wghts(k,ic,:))
              iqi(iy,iz,ix) = SUM(cvals(:,4)*wghts(k,ic,:))
              iqs(iy,iz,ix) = SUM(cvals(:,5)*wghts(k,ic,:))
              iqg(iy,iz,ix) = SUM(cvals(:,6)*wghts(k,ic,:))
              iua(iy,iz,ix) = SUM(cvals(:,7)*wghts(k,ic,:))
              iva(iy,iz,ix) = SUM(cvals(:,8)*wghts(k,ic,:))
              iwa(iy,iz,ix) = SUM(cvals(:,9)*wghts(k,ic,:))
              ita(iy,iz,ix) = SUM(cvals(:,10)*wghts(k,ic,:))
              ipress(iy,iz,ix) = SUM(cvals(:,11)*wghts(k,ic,:))
           
              IF (printdbg ) PRINT *, '       27values: ', cvals(:,1), ' 27wghts: ', wghts(k,ic,:),   &
                ' int. value: ', iqv(iy,iz,ix)
  
            END DO loop_ptsijk
          END DO loop_z

    DO ia=ida, eda
      DO ih=idh, edh
        DO id=idd, edd
          printdbg = .FALSE.
          IF (dbg .AND. MOD(ia+ih+id, 10) == 0 .AND. ahdfound < 10) THEN
            ahdfound = ahdfound + 1
            printdbg = .TRUE.
          END IF

          ! Avoiding masked values
          IF (qv(ia,ih,id) > 0.9*fillValueR) CYCLE
          
          ! Filtering to 0.
          dqv = qv(ia,ih,id)*1.d0
          IF (qc(ia,ih,id) < qwater_constants%q_clou_min) THEN
            dqc = 0.d0
          ELSE
            dqc = qc(ia,ih,id)*1.d0
          END IF
          IF (qr(ia,ih,id) < qwater_constants%q_rain_min) THEN
            dqr = 0.d0
          ELSE
            dqr = qr(ia,ih,id)*1.d0
          END IF
          
          IF (qci(ia,ih,id) < qwater_constants%q_ice_min) THEN
            dqci = 0.d0
          ELSE
            dqci = qci(ia,ih,id)*1.d0
          END IF
          
          IF (qs(ia,ih,id) < qwater_constants%q_snow_min) THEN
            dqs = 0.d0
          ELSE
            dqs = qs(ia,ih,id)*1.d0
          END IF
          
          IF (qg(ia,ih,id) < qwater_constants%q_grau_min) THEN
            dqg = 0.d0
          ELSE
            dqg = qg(ia,ih,id)*1.d0
          END IF
          
          ! Avoiding wasting time
          IF (dqr*dqs*dqg == 0.d0) THEN
            Zh_r_out(ia,ih,id,1:3) = 0.
            Zdr_r_out(ia,ih,id,1:3) = 0.
            KDP_r_out(ia,ih,id,1:3) = 0.
            Zh_s_out(ia,ih,id,1:3) = 0.
            Zdr_s_out(ia,ih,id,1:3) = 0.
            KDP_s_out(ia,ih,id,1:3) = 0.
            Zh_g_out(ia,ih,id,1:3) = 0.
            Zdr_g_out(ia,ih,id,1:3) = 0.
            KDP_g_out(ia,ih,id,1:3) = 0.
            Zh_i_out(ia,ih,id,1:3) = 0.
            Zdr_i_out(ia,ih,id,1:3) = 0.
            KDP_i_out(ia,ih,id,1:3) = 0.
            vt_rain(ia,ih,id) = 0.
            vt_snow(ia,ih,id) = 0.
            vt_grau(ia,ih,id) = 0.
            wt(ia,ih,id,1:3) = 0.
            vr(ia,ih,id,1:3) = 0.
            CYCLE
          END IF
          du = u(ia,ih,id)*1.d0
          dv = v(ia,ih,id)*1.d0
          dw = w(ia,ih,id)*1.d0
          dt = t(ia,ih,id)*1.d0
          dp = p(ia,ih,id)*1.d0
          
          IF (printdbg) THEN
            PRINT *,'    ' // TRIM(fname) // ": ia,ih,id:", ia,ih,id
            PRINT *,'      qv: ', qv(ia,ih,id), ' qc: ', qc(ia,ih,id), ' qr: ', qr(ia,ih,id),         &
              ' qci: ', qci(ia,ih,id), ' qs: ', qs(ia,ih,id), ' qg: ', qg(ia,ih,id)
            PRINT *,'      ta: ', t(ia,ih,id), ' p: ', p(ia,ih,id), ' ua: ', u(ia,ih,id), ' va: ',    &
              v(ia,ih,id), ' wa: ', w(ia,ih,id)
          END IF
          
          ! Polarimetrics_fromLUT_1D(qx, iz, dimv, q_log_vec, Nbnds, Nelevs, Nfrq, Zh_x, Zdr_x,      &
          !    KDP_x, Zh_out, Zdr_out, KDP_out)
          CALL Polarimetrics_fromLUT_1D(dqr, ih, qwater_constants%dimv, q_rain_vec, Nbands, Nlutelevs,&
            Nlutfreqs, Zh_RAIN, Zdr_RAIN, KDP_RAIN, dZh_q_out, dZdr_q_out, dKDP_q_out)
          Zh_r_out(ia,ih,id,1:3) = REAL(dZh_q_out)
          Zdr_r_out(ia,ih,id,1:3) = REAL(dZdr_q_out)
          KDP_r_out(ia,ih,id,1:3) = REAL(dKDP_q_out)
          CALL Polarimetrics_fromLUT_1D(dqs, ih, qwater_constants%dimv, q_snow_vec, Nbands, Nlutelevs,&
            Nlutfreqs, Zh_SNOW, Zdr_SNOW, KDP_SNOW, dZh_q_out, dZdr_q_out, dKDP_q_out)
          Zh_s_out(ia,ih,id,1:3) = REAL(dZh_q_out)
          Zdr_s_out(ia,ih,id,1:3) = REAL(dZdr_q_out)
          KDP_s_out(ia,ih,id,1:3) = REAL(dKDP_q_out)
          CALL Polarimetrics_fromLUT_1D(dqg, ih, qwater_constants%dimv, q_grau_vec, Nbands, Nlutelevs,&
            Nlutfreqs, Zh_GRAU, Zdr_GRAU, KDP_GRAU, dZh_q_out, dZdr_q_out, dKDP_q_out)
          Zh_g_out(ia,ih,id,1:3) = REAL(dZh_q_out)
          Zdr_g_out(ia,ih,id,1:3) = REAL(dZdr_q_out)
          KDP_g_out(ia,ih,id,1:3) = REAL(dKDP_q_out)
          CALL Polarimetrics_fromLUT_1D(dqci, ih, qwater_constants%dimv, q_ice_vec, Nbands, Nlutelevs,&
            Nlutfreqs, Zh_ICE, Zdr_ICE, KDP_ICE, dZh_q_out, dZdr_q_out, dKDP_q_out)
          Zh_i_out(ia,ih,id,1:3) = REAL(dZh_q_out)
          Zdr_i_out(ia,ih,id,1:3) = REAL(dZdr_q_out)
          KDP_i_out(ia,ih,id,1:3) = REAL(dKDP_q_out)

          IF (printdbg) THEN
            DO ibnd=1, Nbands
              PRINT *, '    ibnd:', ibnd, ' RAIN Zh:', Zh_r_out(ia,ih,id,ibnd), ' Zdr: ',             &
                Zdr_r_out(ia,ih,id,ibnd), ' KDP: ', KDP_r_out(ia,ih,id,ibnd)
            END DO
          END IF
          
          ! Terminal velocities          
          CALL calc_tvelocities_1D(dqr, dqs, dqg, dp, dt, dvt_rain, dvt_snow, dvt_grau)
          vt_rain(ia,ih,id) = REAL(dvt_rain)
          vt_snow(ia,ih,id) = REAL(dvt_snow)
          vt_grau(ia,ih,id) = REAL(dvt_grau)
          IF (printdbg) PRINT *, '    terminal velocities rain:', vt_rain(ia,ih,id), ' snow: ',       &
            vt_snow(ia,ih,id), ' graupel: ', vt_grau(ia,ih,id)
          
          ! Calcualte radial velocities ... (not full spectrum)
          ! The average radial velocity vrad is the power-weighted sum of the projections of 
          ! U (eastward wind component), V (northward wind component), W (vertical wind component), 
          ! and vt, the hydrometeor terminal velocity, onto the axis of the radar beam defined by 
          ! angle and azimuth
          dwt = 0.D0
          dvr = 0.D0
          IF (doppler_scheme == 0) THEN
            doppler_scheme_name = 'massweightedmeanvt'
            
            !  Reflectivity weighted terminal velocity
            !    'Same approximation as wrf2radar'
            ! Radial velocity (calculated geometrically)
            cos_el = DCOS(razimuths1D(ih))
            sin_el = DSIN(razimuths1D(ih))
            sin_an = DSIN(rangles1D(ia))
            cos_an = DCOS(rangles1D(ia))
            DO ibnd=1, Nbands
              denom = Zh_r_out(ia,ih,id,ibnd) + Zh_s_out(ia,ih,id,ibnd) + Zh_g_out(ia,ih,id,ibnd) 
              IF (denom /= 0.D0) THEN
                dwt = (dvt_rain*Zh_r_out(ia,ih,id,ibnd) + dvt_snow*Zh_s_out(ia,ih,id,ibnd) +          &
                  dvt_grau*Zh_g_out(ia,ih,id,ibnd)) / denom
                wt(ia,ih,id,ibnd) = REAL(dwt)
              END IF
 
              ! Radial velocity
              dvr = du*cos_el*sin_an + dv*cos_el*cos_an + (dw-wt(ia,ih,id,ibnd))*sin_el
              vr(ia,ih,id,ibnd) = REAL(dvr)
            END DO
            
            IF (printdbg) THEN
              PRINT *, '    wt:', wt(ia,ih,id,:)
              PRINT *, '    radial terminal velocity:', vr(ia,ih,id,:)
            END IF
            
          ELSE IF (doppler_scheme == 1) THEN
            ! Get PSD integrated backscattering cross section int = v S/S
            CALL ErrMsg('Still need to add this functionality. Get the backscattering cross ' //      &
              'section for everygrid point', fname, -1)
          END IF
        END DO
      END DO
    END DO
          
          
          
            
          IF (dbg_level >= 75 .AND. ijdbg .AND. tdbg==0) THEN
            WRITE(msg,*) '  ' // TRIM(fname) // ' ',i,' ,',j,' : final checks _______' 
            CALL wrf_debug(75,msg)
          END IF
          
          
          
          ijdbg = .FALSE.

        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN
  
  END SUBROUTINE simradar_calc

#endif
END MODULE module_simradar
