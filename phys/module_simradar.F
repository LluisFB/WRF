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
!                      n_moist,                                                                        &
!                      m_qv, m_qc, m_qr, m_qs, m_qi, m_qg, m_qh,                                       &
!                      moist,                                                                          &
                      qv, qc, qr, qs, qi, qg, qh,                                                     &
                      lon, lat,                                                                       &
                      mu, mub, sina, cosa,                                                            &
                      u, v, w,                                                                        &
                      p, pb,                                                                          &
                      th300, ph, phb,                                                                 &
                      nradar, ngrid_radar,                                                            &
                      e_angle, e_azimuth, e_dist,                                                     &
                      coord2, coord3, coord4, around2, around3,                                       &
                      ngrid3_radar, nrd,                                                              &
                      nrnn, nrnncoord, nrcoord, nr4coord,                                             &
                      radars_angle, radars_azimuth, radars_dist,                                      &
                      ijradardbg,                                                                     &
                      radarlons, radarlats,                                                           &
                      radarhgts, radarlhgts,                                                          &
                      gatehdist, inradar, radnum,                                                     &
                      ngridradar, ijgridradar, ijdistgridradar,                                       &
                      lutgrid,                                                                        &

                      ! Output
                      !!! Optional. Modify Registry/registry.simradar and recompile
                      !!SIMRADARop!!  iqvm, iqcm, iqrm, iqim,                                           &
                      !!SIMRADARop!!  iqsm, iqgm, iqhm,                                                 &
                      !!SIMRADARop!!  iuam, ivam, iwam, itam, ipressm,                                  &
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
    !   NOTE radars are in polar space:
    !     dx: distance    dy: angle    dz: azimuth
    !     3Dvariables (dy,dz,dx),        lon,lats(dy,dx),        hgt(dx,dz)

    ! There is no domain-related differences in radar configuration
    USE module_configure, ONLY: model_config_rec
    USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
    USE module_configure 
    USE module_simradar_tools
    USE module_model_constants
!    USE module_diagvar_simradar

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
! mu: perturbation dry air mass in column [Pa]
! mub: base dry air mass in column [Pa]
! sina: local sine of map rotation [1]
! cosa: local cosine of map rotation [1]
! u: eastward wind speed [ms-1]
! v: northward wind speed [ms-1]
! wa: vertical wind speed [ms-1]
! p: 3D pressure perturbation [Pa]
! pb: 3D pressure base state [Pa]
! th300: potential air temperature -300. [K]
! ph: 3D gepotential perturbation [m2s-2]
! phb: 3D gepotential base state [m2s-2]

!! Radar variables
! nradar: amount of radars
! e_angle: amount of angle beams (horizontal) in the radar space
! e_azimuth: amount of azimuth angles (vertical) in the radar space (same for all radars)
! azimuths: azimuth angles
! e_dist: amount of distances (horizontal) in the radar space
! coord2: amount of 2-dimensional coordinates (2)
! coord3: amount of 3-dimensional coordinates (3)
! coord4: amount of coordinates for grid/cell 2-dimensional boundaries (4)
! around2: amoount of points of 3-dimensional surroundings (9) 
! around3: amoount of points of 3-dimensional surroundings (27) 
! ngrid_radar: maximum amount of radar cells in a WRF grid
! ngrid3_radar: maximum amount of 3-dimensional radar cells in a WRF grid
! inradar: pertinence of the grid cell to any radar (1: yes)
! radarnum: pertinence of the grid cell inside a given radar area [1]
! radars_dist: ground radar gate distances to their center [m]
! radars_angle: ground radar gate horizontal angles to emit beams [m]
! radars_azimuth: ground radar gate vertical angles to emit beams [m]
! radarlons: ground radar gate longitude [degree_east]
! radarlats: ground radar gate latitude [degree_north]
! radarhgts: ground radar gate height [m]
! radarlhgts: ground radar gate local height [m]
! radarhdist: ground radar horizontal distance of the grid point [m]
! gatehdist: distacnce of the grid points to the position of the radar [m]
! ngridradar: amount of radar cells inside the grid point [1]
! ijgridradar: coordinates of radar cells inside the grid point [-]
! ijdistgridradar: distance of radar cells inside the grid point to its center [m]

!! Compressed dimensions
! nrnn = nradar*ngrid_radar
! nrcoord = nradar*coord2
! nrnncoord = nradar*ngrid_radar*coord2
! nncoord = ngrid_radar*coord2
! nr4coord = nradar*coord4
! nrd: nradar*e_dist

! timedbg: frecuency in minutes of time-debugging
! londbg: longitude in degrees of grid-point debugging (-999., not used)
! latdbg: latitude in degrees of grid-point debugging (-999., not used)

! Output
!-------
!! Optional. !!SIMRADARop!!
!!   Modify Registry/registry.simradar and recompile
! iqvm: 3-dimensional radar interpolated variable of water vapour mixing ratio [kgkg-1]
! iqcm: 3-dimensional radar interpolated variable of cloud mixing ratio [kgkg-1]
! iqrm: 3-dimensional radar interpolated variable of rain mixing ratio [kgkg-1]
! iqim: 3-dimensional radar interpolated variable of ice mixing ratio [kgkg-1]
! iqsm: 3-dimensional radar interpolated variable of snow mixing ratio [kgkg-1]
! iqgm: 3-dimensional radar interpolated variable of graupel mixing ratio [kgkg-1]
! iqhm: 3-dimensional radar interpolated variable of hail mixing ratio [kgkg-1]
! iuam: 3-dimensional radar interpolated variable of eastward wind speed [ms-1]
! ivam: 3-dimensional radar interpolated variable of northward wind speed [ms-1]
! iwam: 3-dimensional radar interpolated variable of vertical wind speed [ms-1]
! itam: 3-dimensional radar interpolated variable of air temperature [K]
! ipressm: 3-dimensional radar interpolated variable of air pressure [Pa]
!! END Optional
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
                                            num_tiles
!                                            num_tiles,                                                &
!                                            n_moist


    INTEGER, DIMENSION(num_tiles), INTENT(in)            :: i_start, i_end, j_start, j_end
!    INTEGER, INTENT(in)                                  :: m_qv, m_qc, m_qr, m_qs, m_qi, m_qg, m_qh

    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: lon, lat
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: mu, mub, sina, cosa
    REAL, DIMENSION(ims:ime,kms:kme,jms:jme), OPTIONAL,                                               &
      INTENT(in)                                         :: qv, qc, qr, qs, qi, qg, qh
    REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(in) :: u, v, w, p, pb, th300, ph, phb

    REAL, INTENT(in)                                     :: DT, XTIME
    REAL, INTENT(in)                                     :: curr_secs2
    LOGICAL, INTENT(in)                                  :: is_restart, is_output_time
    INTEGER, INTENT(in)                                  :: timedbg
    REAL, INTENT(in)                                     :: londbg, latdbg
    INTEGER, INTENT(in)                                  :: idbg, jdbg

    INTEGER, INTENT(inout)                               :: nsteps
    INTEGER, INTENT(in)                                  :: coord2, coord3, coord4, around2, around3
  
! Specific Radar variables
    INTEGER, INTENT(in)                                  :: nradar, ngrid_radar
    INTEGER, INTENT(in)                                  :: e_azimuth
    INTEGER, INTENT(in)                                  :: e_angle, e_dist
!    REAL, DIMENSION(nradar), INTENT(in)                  :: gatelon, gatelat, gatehgt, distmax,       &
!      heightmax, radards
!    INTEGER, DIMENSION(nradar), INTENT(in)               :: gatetype, gateband, radinterp, em_ray
!    INTEGER, DIMENSION(nradar), INTENT(in)               :: em_ray
!    REAL, INTENT(in)                                     :: radii_ddist
!    REAL, INTENT(in)                                     :: missing_value
    REAL, DIMENSION(e_angle, nradar), INTENT(in)         :: radars_angle
    REAL, DIMENSION(e_azimuth, nradar), INTENT(in)       :: radars_azimuth
    REAL, DIMENSION(e_dist, nradar), INTENT(in)          :: radars_dist
    REAL(kind=rr8)                                       :: a, f
    CHARACTER(len=30)                                    :: elipsn
    ! Compressed dimensions
    INTEGER, INTENT(in)                                  :: nrnn, nrnncoord, nrcoord, nr4coord 
    INTEGER, INTENT(in)                                  :: ngrid3_radar, nrd
    INTEGER, DIMENSION(coord2), INTENT(in)               :: ijradardbg
    REAL, DIMENSION(ims:ime, jms:jme, nradar),                                                        &
      INTENT(in)                                         :: gatehdist
    INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(in)     :: inradar
    INTEGER, DIMENSION(ims:ime, jms:jme, nradar),                                                     &
      INTENT(in)                                         :: radnum
    REAL, DIMENSION(e_angle, e_dist, nradar), INTENT(in) :: radarlons, radarlats
    REAL, DIMENSION(e_azimuth,e_dist,nradar), INTENT(in) :: radarhgts, radarlhgts
!    REAL,DIMENSION(e_angle,e_dist,nr4coord), INTENT(in)  :: radarlons_bnds, radarlat_bnds,            &
!    REAL, DIMENSION(e_azimuth,e_dist,nr4coord), INTENT(in)                                            &
!                                                         :: radarhgt_bnds
    INTEGER, DIMENSION(ims:ime,jms:jme, nradar),                                                      &
      INTENT(in)                                         :: ngridradar
    INTEGER, DIMENSION(ims:ime, jms:jme, nrnncoord),                                                  &
      INTENT(in)                                         :: ijgridradar
    REAL, DIMENSION(ims:ime, jms:jme, nrnn), INTENT(in)  :: ijdistgridradar
    
    ! LUT space
    TYPE(lut_vars), INTENT(in)                           :: lutgrid

    ! Output 
    !! optional !! SIMRADARop. Uncomment, modify Regisry/registry.simradar and recompile
    !!REAL, DIMENSION(e_angle,e_azimiuth,nrd), INTENT(out)                                               &
    !!                                                     :: iqvm, iqcm, iqrm, iqim, iqsm, iqgm, iqhm,  &
    !!  iuam, ivam, iwam, itam, ipressm

    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: zh_rain, zdr_rain, kdp_rain
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: zh_snow, zdr_snow, kdp_snow
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: zh_grau, zdr_grau, kdp_grau
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: zh_ice, zdr_ice, kdp_ice
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: zh_hail, zdr_hail, kdp_hail
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: zh_total, zdr_total, kdp_total
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: tv_rain, tv_snow, tv_grau
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: tv_ice, tv_hail
    REAL, DIMENSION(e_angle, e_azimuth, nrd),INTENT(out) :: tv_wgt, vr
    
    ! LOCAL  VAR
    INTEGER                                              :: i, j, k, l, m, n, its, ite, jts, jte, ij
    INTEGER                                              :: idp, jdp, dimz, dimz1, dz1
    INTEGER                                              :: i1, i_1, j1, j_1, k1, k2, tdbg
    INTEGER                                              :: ii1, ei1, ij1, ej1
    INTEGER, DIMENSION(num_tiles)                        :: i2, j2
    INTEGER                                              :: inrnn, inrnncoord, inncoord, inrcoord,    &
      inr4coord, inrd

    REAL                                                 :: xtimep
    REAL                                                 :: prv
    
    LOGICAL, EXTERNAL                                    :: wrf_dm_on_monitor
    CHARACTER(len=256)                                   :: msg
    LOGICAL                                              :: ijdbg, cdone
    INTEGER                                              :: dbg_level
 
    INTEGER                                              :: ir, ia, ih, id, iband
    INTEGER                                              :: ix, iy, iz
    INTEGER, DIMENSION(ngrid_radar, coord2)              :: ijgridradar_NOc
    REAL, DIMENSION(ngrid_radar, nradar)                 :: ijdistgridradar_NOc

    ! radar 3-dimensional
    INTEGER                                              :: ijkNtot
    INTEGER, DIMENSION(kms:kme)                          :: ijkngridradar
    INTEGER, DIMENSION(kms:kme, coord3, ngrid3_radar)    :: ijkgridradar
    REAL, DIMENSION(kms:kme, around3, ngrid3_radar)      :: wghts
    
    INTEGER                                              :: ii2, ii3, ic
    INTEGER                                              :: ei2, ei3
    INTEGER                                              :: i123
    INTEGER                                              :: ijkfound
    REAL, DIMENSION(around3,12)                          :: cvals
    CHARACTER(len=1), DIMENSION(3)                       :: availrband

    REAL, DIMENSION(ims:ime,jms:jme)                     :: mutot
    REAL, DIMENSION(ims:ime,kms:kme,jms:jme)             :: zg, press, ta, unua, unva, ua, va
    
    REAL                                                 :: iqv, iqc, iqr, iqi, iqs, iqg, iqh, iua,   &
      iva, iwa, ita, ipress
    REAL(kind=rr8)                                       :: dqv, dqc, dqr, dqi, dqs, dqg, dqh, dua,   &
      dva, dwa, dta, dpress
    REAL(kind=rr8)                                       :: dZh_q, dZdr_q, dKDP_q
    REAL                                                 :: zh_rain_NOc, zdr_rain_NOc, kdp_rain_NOc
    REAL                                                 :: zh_snow_NOc, zdr_snow_NOc, kdp_snow_NOc
    REAL                                                 :: zh_grau_NOc, zdr_grau_NOc, kdp_grau_NOc
    REAL                                                 :: zh_ice_NOc, zdr_ice_NOc, kdp_ice_NOc
    REAL                                                 :: zh_hail_NOc, zdr_hail_NOc, kdp_hail_NOc
    REAL                                                 :: zh_total_NOc, zdr_total_NOc, kdp_total_NOc
    REAL(kind=rr8)                                       :: dtv_rain, dtv_snow, dtv_grau, dtv_ice,    &
      dtv_hail
    REAL                                                 :: tv_rain_NOc, tv_snow_NOc, tv_grau_NOc
    REAL                                                 :: tv_ice_NOc, tv_hail_NOc
    REAL(kind=rr8)                                       :: dangle, dazimuth, ddist
    REAL(kind=rr8)                                       :: cos_el, sin_el, cos_an, sin_an
    REAL(kind=rr8)                                       :: denom
    REAL(kind=rr8)                                       :: dwt, dvr
    REAL                                                 :: tv_wgt_NOc, vr_NOc, wt_NOc
    LOGICAL                                              :: dbg, printdbg
    CHARACTER(len=S32)                                   :: fname

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
    IF (dbg_level > 75) dbg = .TRUE.
    cdone = .FALSE.

    ! L. Fita: NOTE about efficiency
    !   I am not sure if this is the more efficient way:
    !     1.- IF --> 2.- repeat DO
    !     1.- DO -->  2.- IF

    dimz = kme - kms + 1
    dimz1 = dimz + 1
    dz1 = dimz - 1
    
    IF (ALL(inradar == 0)) THEN
      !!! Optional. Modify Registry/registry.simradar and recompile
      !!SIMRADARop!!iqvm = model_config_rec%missing_value
      !!SIMRADARop!!iqcm = model_config_rec%missing_value
      !!SIMRADARop!!iqrm = model_config_rec%missing_value
      !!SIMRADARop!!iqim = model_config_rec%missing_value
      !!SIMRADARop!!iqsm = model_config_rec%missing_value
      !!SIMRADARop!!iqgm = model_config_rec%missing_value
      !!SIMRADARop!!iqhm = model_config_rec%missing_value
      !!SIMRADARop!!iuam = model_config_rec%missing_value
      !!SIMRADARop!!ivam = model_config_rec%missing_value
      !!SIMRADARop!!iwam = model_config_rec%missing_value
      !!SIMRADARop!!itam = model_config_rec%missing_value
      !!SIMRADARop!!ipressm = model_config_rec%missing_value

      zh_rain = model_config_rec%missing_value
      zdr_rain = model_config_rec%missing_value
      kdp_rain = model_config_rec%missing_value
      zh_snow = model_config_rec%missing_value
      zdr_snow = model_config_rec%missing_value
      kdp_snow = model_config_rec%missing_value
      zh_grau = model_config_rec%missing_value
      zdr_grau = model_config_rec%missing_value
      kdp_grau = model_config_rec%missing_value
      zh_ice = model_config_rec%missing_value
      zdr_ice = model_config_rec%missing_value
      kdp_ice = model_config_rec%missing_value
      zh_hail = model_config_rec%missing_value
      zdr_hail = model_config_rec%missing_value
      kdp_hail = model_config_rec%missing_value
      zh_total = model_config_rec%missing_value
      zdr_total = model_config_rec%missing_value
      kdp_total = model_config_rec%missing_value
      tv_rain = model_config_rec%missing_value
      tv_snow = model_config_rec%missing_value
      tv_grau = model_config_rec%missing_value
      tv_ice = model_config_rec%missing_value
      tv_hail = model_config_rec%missing_value
      tv_wgt = model_config_rec%missing_value
      vr = model_config_rec%missing_value    
      RETURN
    ELSE        
      !! Getting standard atmospheric variables
      DO ij = 1 , num_tiles
        DO i = i_start(ij), i_end(ij)
          DO j = j_start(ij), j_end(ij)
            mutot(i,j) = mu(i,j) + mub(i,j)
            zg(i,:,j) = (ph(i,:,j) + phb(i,:,j))/9.81
            press(i,:,j) = p(i,:,j) + pb(i,:,j)
            ta(i,:,j) = (th300(i,:,j) + 300.)*(press(i,:,j)/p0)**(rcp)
            unua(i,:,j) = 0.5*(u(i+1,:,j) + u(i,:,j))
            unva(i,:,j) = 0.5*(u(i,:,j+1) + u(i,:,j))
            ua(i,:,j) = unua(i,:,j)*cosa(i,j) - unva(i,:,j)*sina(i,j)
            va(i,:,j) = unua(i,:,j)*cosa(i,j) + unva(i,:,j)*sina(i,j)
          END DO
        END DO
      END DO 
    END IF

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
    radars_loop: DO ir=1, nradar
    
      ijkfound = 0

      IF (TRIM(model_config_rec%radar_band(ir)) == 'S') THEN
        iband = 1
      ELSE IF (TRIM(model_config_rec%radar_band(ir)) == 'C') THEN
        iband = 2
      ELSE IF (TRIM(model_config_rec%radar_band(ir)) == 'X') THEN
        iband = 3
      ELSE
        availrband(1) = 'S'
        availrband(2) = 'C'
        availrband(3) = 'X'
        CALL ErrMsgSAvail("radar band '" // TRIM(model_config_rec%radar_band(ir)) // "' not ready !!",&
          3, availrband, fname)
      END IF       
    
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      tiles_loop: DO ij = 1 , num_tiles
        DO i = i_start(ij), i_end(ij)
          DO j = j_start(ij), j_end(ij)
            IF (ij == 1 .AND. i == i2(ij) .AND. j == j2(ij)) ijdbg = .TRUE.
            i1 = MIN(i+1,i_end(ij))
            i_1 = MAX(i-1,i_start(ij))
            j1 = MIN(j+1,j_end(ij))
            j_1 = MAX(j-1,j_start(ij))
          
            IF (ngridradar(i,j,ir) == 0) CYCLE
                    
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
              ijdistgridradar_NOc(n,ir) = ijdistgridradar(i,j,inrnn)
              DO l=1, coord2
                inrnncoord = (ir-1)*ngrid_radar*2 + (n-1)*2 + l
                ijgridradar_NOc(n,l) = ijgridradar(i,j,inrnncoord)
              END DO
            END DO

            ! Computing 27 weights around each WRF i,j,k for all radar cells
            CALL weights_hcoordsR27around_polar_ellipsoide_1D(i, j, ii1, ei1, ij1, ej1, kms, kme, a,  &
              f, elipsn, lon(ii1:ei1,ij1:ej1), lat(ii1:ei1,ij1:ej1), zg(ii1:ei1,kms:kme,ij1:ej1),     &
              model_config_rec%missing_value, ngrid_radar, ngridradar(i,j,ir), coord2,                &
              ijgridradar_NOc, around2, 1, e_dist, 1, e_angle, 1, e_azimuth, radarlons(:,:,ir),       &
              radarlats(:,:,ir), radarhgts(:,:,ir), model_config_rec%missing_value, ngrid3_radar,     &
              ijkngridradar, coord3, ijkgridradar, ijkNtot, around3, wghts, dbg)
           
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
  
              IF (printdbg) PRINT *,'  ' // TRIM(fname) // ": ", ijkfound, ' i,j,k ', i,j,k,          &
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
                      cvals(i123,7) = qh(i+ii1,k+ii3,j+ii2)
                      cvals(i123,8) = ua(i+ii1,k+ii3,j+ii2)
                      cvals(i123,9) = va(i+ii1,k+ii3,j+ii2)
                      cvals(i123,10) = w(i+ii1,k+ii3,j+ii2)
                      cvals(i123,11) = ta(i+ii1,k+ii3,j+ii2)
                      cvals(i123,12) = press(i+ii1,k+ii3,j+ii2)
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
                      cvals(i123,7) = qh(i+ii1,k+ii3,j+ii2)
                      cvals(i123,8) = ua(i+ii1,k+ii3,j+ii2)
                      cvals(i123,9) = va(i+ii1,k+ii3,j+ii2)
                      cvals(i123,10) = w(i+ii1,k+ii3,j+ii2)
                      cvals(i123,11) = ta(i+ii1,k+ii3,j+ii2)
                      cvals(i123,12) = press(i+ii1,k+ii3,j+ii2)
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
                      cvals(i123,7) = qh(i+ii1,k+ii3,j+ii2)
                      cvals(i123,8) = ua(i+ii1,k+ii3,j+ii2)
                      cvals(i123,9) = va(i+ii1,k+ii3,j+ii2)
                      cvals(i123,10) = w(i+ii1,k+ii3,j+ii2)
                      cvals(i123,11) = ta(i+ii1,k+ii3,j+ii2)
                      cvals(i123,12) = press(i+ii1,k+ii3,j+ii2)
                      i123 = i123 + 1
                    END DO
                  END DO
                END DO
              END IF
          
              ! Removing too small values...
              DO ii1=1, around3
                DO ii2=1, 12
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
              
                IF (printdbg) PRINT *,'    cell: ', ic, ' coords; ', ijkgridradar(k,ic,:),            &
                  ' || lon: ', radarlons(iy,ix,ir), ' lat: ', radarlats(iy,ix,ir), ' hgt: ',          &
                  radarhgts(iz,ix,ir), ' amount of values: ', i123
                ! Need to filter Nan and missing values ?
                !ijkvals(i,j,k,ic) = 0.
                !DO ijk=1,acoords3
                !  ijkvals(i,j,k,ic) = ijkvals(i,j,k,ic) + cvals(ijk)*wghts(i,j,k,ic,ijk)
                !END DO
              
                ! Let's be memory constrained (might be ofered for output checks?)
                iqv = SUM(cvals(:,1)*wghts(k,ic,:))
                iqc = SUM(cvals(:,2)*wghts(k,ic,:))
                iqr = SUM(cvals(:,3)*wghts(k,ic,:))
                iqi = SUM(cvals(:,4)*wghts(k,ic,:))
                iqs = SUM(cvals(:,5)*wghts(k,ic,:))
                iqg = SUM(cvals(:,6)*wghts(k,ic,:))
                iqh = SUM(cvals(:,7)*wghts(k,ic,:))
                iua = SUM(cvals(:,8)*wghts(k,ic,:))
                iva = SUM(cvals(:,9)*wghts(k,ic,:))
                iwa = SUM(cvals(:,10)*wghts(k,ic,:))
                ita = SUM(cvals(:,11)*wghts(k,ic,:))
                ipress = SUM(cvals(:,12)*wghts(k,ic,:))
          
                IF (printdbg ) PRINT *, '       27values: ', cvals(:,1), ' 27wghts: ', wghts(k,ic,:), &
                  ' int. value: ', iqv
                
                !! Computing radar diagnostics
                ! Avoiding masked values
                IF (iqv > 0.9*fillValueR) CYCLE
            
                ! Filtering to 0.
                dqc = iqc*1.d0
                IF (iqc .lt. lutgrid%q_clou_min) THEN
                  dqc = 0.d0
                ELSE
                  dqc = iqc*1.d0
                END IF

                IF (iqr .lt. lutgrid%q_rain_min) THEN
                  dqr = 0.d0
                ELSE
                  dqr = iqr*1.d0
                END IF
      
                IF (iqi .lt. lutgrid%q_ice_min) THEN
                  dqi = 0.d0
                ELSE
                  dqi = iqi*1.d0
                END IF
          
                IF (iqs .lt. lutgrid%q_snow_min) THEN
                  dqs = 0.d0
                ELSE
                  dqs = iqs*1.d0
                END IF
          
                IF (iqg .lt. lutgrid%q_grau_min) THEN
                  dqg = 0.d0
                ELSE
                  dqg = iqg*1.d0
                END IF
                    
                IF (iqh .lt. lutgrid%q_hail_min) THEN
                  dqh = 0.d0
                ELSE
                  dqh = iqh*1.d0
                END IF

                ! Avoiding wasting time
                IF (dqr+dqs+dqg+dqh == 0.d0) THEN
                  zh_rain_NOc = 0.
                  zdr_rain_NOc = 0.
                  kdp_rain_NOc = 0.               
                  zh_snow_NOc = 0.
                  zdr_snow_NOc = 0.
                  kdp_snow_NOc = 0.
                  zh_grau_NOc = 0.
                  zdr_grau_NOc = 0.
                  kdp_grau_NOc = 0.
                  zh_ice_NOc = 0.
                  zdr_ice_NOc = 0.
                  kdp_ice_NOc = 0.
                  zh_hail_NOc = 0.
                  zdr_hail_NOc = 0.
                  kdp_hail_NOc = 0.
                  zh_total_NOc = 0.
                  zdr_total_NOc = 0.
                  kdp_total_NOc = 0.
                  tv_rain_NOc = 0.
                  tv_snow_NOc = 0.
                  tv_grau_NOc = 0.
                  wt_NOc = 0.
                  vr_NOc = 0.
                  CYCLE
                END IF
                dua = iua*1.d0
                dva = iva*1.d0
                dwa = iwa*1.d0
                dta = ita*1.d0
                dpress = ipress*1.d0
              
                IF (printdbg) THEN
                  PRINT *,'    ' // TRIM(fname) // " radar space: ia,ih,id:", iy,iz,ix
                  PRINT *,'      qv: ', iqv, ' qc: ', iqc, ' qr: ', iqr, ' qci: ', iqi, ' qs: ',      &
                    iqs, ' qg: ', iqg, ' qh: ', iqh
                  PRINT *,'      ta: ',ita, ' press: ', ipress, ' ua: ', iua, ' va: ', iva, ' wa: ',  &
                    iwa
                END IF
              
                ! Polarimetrics_fromLUT_1D(qx, iz, dimv, q_log_vec, Nbnds, Nelevs, Nfrq, Zh_x, Zdr_x,  &
                !    KDP_x, Zh_out, Zdr_out, KDP_out)
                ! Rain
                CALL Polarimetrics_fromLUT_1D(dqr, iz, iband, lutgrid%dimv, lutgrid%q_rain_vec,        &
                  lutgrid%Nlutqxbins, lutgrid%Zh_RAIN(iband,iz,:), lutgrid%Zdr_RAIN(iband,iz,:),       &
                  lutgrid%KDP_RAIN(iband,iz,:), dZh_q, dZdr_q, dKDP_q)
                zh_rain_NOc = REAL(dZh_q)
                zdr_rain_NOc = REAL(dZdr_q)
                kdp_rain_NOc = REAL(dKDP_q)

                ! Snow
                CALL Polarimetrics_fromLUT_1D(dqs, iz, iband, lutgrid%dimv, lutgrid%q_snow_vec,       &
                  lutgrid%Nlutqxbins, lutgrid%Zh_SNOW(iband,iz,:), lutgrid%Zdr_SNOW(iband,iz,:),      &
                  lutgrid%KDP_SNOW(iband,iz,:), dZh_q, dZdr_q, dKDP_q)
                zh_snow_NOc = REAL(dZh_q)
                zdr_snow_NOc = REAL(dZdr_q)
                kdp_snow_NOc = REAL(dKDP_q)

                ! Graupel
                CALL Polarimetrics_fromLUT_1D(dqg, iz, iband, lutgrid%dimv, lutgrid%q_grau_vec,       &
                  lutgrid%Nlutqxbins, lutgrid%Zh_GRAU(iband,iz,:), lutgrid%Zdr_GRAU(iband,iz,:),      &
                  lutgrid%KDP_GRAU(iband,iz,:), dZh_q, dZdr_q, dKDP_q)
                zh_grau_NOc = REAL(dZh_q)
                zdr_grau_NOc = REAL(dZdr_q)
                kdp_grau_NOc = REAL(dKDP_q)

                ! Cloud ice
                CALL Polarimetrics_fromLUT_1D(dqi, iz, iband, lutgrid%dimv, lutgrid%q_ice_vec,        &
                  lutgrid%Nlutqxbins, lutgrid%Zh_ICE(iband,iz,:), lutgrid%Zdr_ICE(iband,iz,:),        &
                  lutgrid%KDP_ICE(iband,iz,:), dZh_q, dZdr_q, dKDP_q)
                zh_ice_NOc = REAL(dZh_q)
                zdr_ice_NOc = REAL(dZdr_q)
                kdp_ice_NOc = REAL(dKDP_q)

                ! Hail
                CALL Polarimetrics_fromLUT_1D(dqh, iz, iband, lutgrid%dimv, lutgrid%q_hail_vec,       &
                  lutgrid%Nlutqxbins, lutgrid%Zh_HAIL(iband,iz,:), lutgrid%Zdr_HAIL(iband,iz,:),      &
                  lutgrid%KDP_HAIL(iband,iz,:), dZh_q, dZdr_q, dKDP_q)
                zh_hail_NOc = REAL(dZh_q)
                zdr_hail_NOc = REAL(dZdr_q)
                kdp_hail_NOc = REAL(dKDP_q)

                IF (printdbg) THEN
                  PRINT *, '    band:', TRIM(model_config_rec%radar_band(ir)), ' RAIN Zh:',           &
                    zh_rain_NOc, ' Zdr: ', zdr_rain_NOc, ' KDP: ', kdp_rain_NOc
                END IF
            
                ! Terminal velocities          
                CALL calc_tvelocities_1D(dqr, dqs, dqg, dqh, dpress, dta, dtv_rain, dtv_snow,         &
                  dtv_grau, dtv_hail)
                tv_rain_NOc = REAL(dtv_rain)
                tv_snow_NOc = REAL(dtv_snow)
                tv_grau_NOc = REAL(dtv_grau)
                tv_hail_NOc = REAL(dtv_hail)
                IF (printdbg) PRINT *, '    terminal velocities rain:', tv_rain_NOc, ' snow: ',       &
                  tv_snow_NOc, ' graupel: ', tv_grau_NOc, ' hail: ', tv_hail_NOc
              
                ! Calcualte radial velocities ... (not full spectrum)
                ! The average radial velocity vrad is the power-weighted sum of the projections of 
                ! U (eastward wind component), V (northward wind component), W (vertical wind component), 
                ! and vt, the hydrometeor terminal velocity, onto the axis of the radar beam defined by 
                ! angle and azimuth
                dwt = 0.D0
                dvr = 0.D0
                IF (model_config_rec%doppler_scheme == 0) THEN
                  !doppler_scheme_name = 'massweightedmeanvt'
                
                  !  Reflectivity weighted terminal velocity
                  !    'Same approximation as wrf2radar'
                  ! Radial velocity (calculated geometrically)
                  dangle = radars_angle(iy,ir)*1.D0 / dRad
                  dazimuth = radars_azimuth(iz,ir)*1.D0 / dRad
                  cos_el = DCOS(dazimuth)
                  sin_el = DSIN(dazimuth)
                  sin_an = DSIN(dangle)
                  cos_an = DCOS(dangle)
                  denom = zh_rain_NOc + zh_snow_NOc + zh_grau_NOc + zh_hail_NOc
                  IF (denom /= 0.D0) THEN
                    dwt = (dtv_rain*zh_rain_NOc + dtv_snow*zh_snow_NOc + dtv_grau*zh_grau_NOc +       &
                      dtv_hail*zh_hail_NOc) / denom
                    wt_NOc = REAL(dwt)
                  END IF
     
                  ! Radial velocity
                  dvr = dua*cos_el*sin_an + dva*cos_el*cos_an + (dwa-wt_NOc)*sin_el
                  vr_NOc = REAL(dvr)
                  
                  IF (printdbg) THEN
                    PRINT *, '    wt:', wt_NOc
                    PRINT *, '    radial terminal velocity:', vr_NOc
                  END IF
                
                ELSE IF (model_config_rec%doppler_scheme == 1) THEN
                  ! Get PSD integrated backscattering cross section int = v S/S
                  CALL ErrMsg('Still need to add this functionality. Get the backscattering cross '// &
                    'section for everygrid point', fname, -1)
                END IF
              
                ! Compressing
                inrd = (ir-1)*e_dist + ix
              
                !! Optional
                !! SIMRADARop !!
                !!iqvm(iy,iz,inrd) = iqv
                !!iqcm(iy,iz,inrd) = iqc
                !!iqrm(iy,iz,inrd) = iqr
                !!iqim(iy,iz,inrd) = iqi
                !!iqsm(iy,iz,inrd) = iqs
                !!iqgm(iy,iz,inrd) = iqg
                !!iqhm(iy,iz,inrd) = iqh
                !!iuam(iy,iz,inrd) = iua
                !!ivam(iy,iz,inrd) = iva
                !!iwam(iy,iz,inrd) = iwa
                !!itam(iy,iz,inrd) = ita
                !!ipressm(iy,iz,inrd) = ipress

                zh_rain(iy,iz,inrd) = zh_rain_NOc
                zdr_rain(iy,iz,inrd) = zdr_rain_NOc
                kdp_rain(iy,iz,inrd) = kdp_rain_NOc
              
                zh_snow(iy,iz,inrd) = zh_snow_NOc
                zdr_snow(iy,iz,inrd) = zdr_snow_NOc
                kdp_snow(iy,iz,inrd) = kdp_snow_NOc
  
                zh_grau(iy,iz,inrd) = zh_grau_NOc
                zdr_grau(iy,iz,inrd) = zdr_grau_NOc
                kdp_grau(iy,iz,inrd) = kdp_grau_NOc

                zh_ice(iy,iz,inrd) = zh_ice_NOc
                zdr_ice(iy,iz,inrd) = zdr_ice_NOc
                kdp_ice(iy,iz,inrd) = kdp_ice_NOc

                zh_hail(iy,iz,inrd) = zh_hail_NOc
                zdr_hail(iy,iz,inrd) = zdr_hail_NOc
                kdp_hail(iy,iz,inrd) = kdp_hail_NOc

                zh_total(iy,iz,inrd) = zh_rain_NOc + zh_snow_NOc + zh_grau_NOc + zh_ice_NOc +         &
                  zh_hail_NOc
                zdr_total(iy,iz,inrd) = zdr_rain_NOc + zdr_snow_NOc + zdr_grau_NOc + zdr_ice_NOc +    &
                  zdr_hail_NOc
                kdp_total(iy,iz,inrd) = kdp_rain_NOc + kdp_snow_NOc + kdp_grau_NOc + kdp_ice_NOc +    &
                  kdp_hail_NOc

                tv_rain(iy,iz,inrd) = tv_rain_NOc
                tv_snow(iy,iz,inrd) = tv_snow_NOc
                tv_grau(iy,iz,inrd) = tv_grau_NOc
                tv_ice(iy,iz,inrd) = tv_ice_NOc
                tv_hail(iy,iz,inrd) = tv_hail_NOc

                tv_wgt(iy,iz,inrd) = wt_NOc
                vr(iy,iz,inrd) = vr_NOc

              END DO loop_ptsijk
            END DO loop_z      
          
            IF (dbg_level >= 75 .AND. ijdbg .AND. tdbg==0) THEN
              WRITE(msg,*) '  ' // TRIM(fname) // ' ',i,' ,',j,' : final checks _______' 
              CALL wrf_debug(75,msg)
            END IF
          
            ijdbg = .FALSE.

          END DO
        END DO
      END DO tiles_loop
      !$OMP END PARALLEL DO

    END DO radars_loop

    RETURN
  
  END SUBROUTINE simradar_calc

#endif
END MODULE module_simradar
