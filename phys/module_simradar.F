MODULE module_simradar
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
                      mu, mub,                                                                        &
                      ph, phb, rho,                                                                   &
                      nradar,                                                                         &
                      e_radlon, e_radlat,                                                             &
                      e_radhgt, e_radnz,                                                              &
                      gatelon, gatelat, gatehgt, gatehwght,                                           &
                      inradar, radnum,                                                                &
                      radards, radinterp, missing_value,                                              &

                     ! Output
                      vr,                                                                             &

                     !! S-band
                     aihtot_s,                                                                        &
                     aih_rain_s, aih_snow_s,                                                          &
                     aih_grau_s, aih_hail_s,                                                          &
                     kdptot_s,                                                                        &
                     kdp_rain_s, kdp_snow_s,                                                          &
                     kdp_grau_s, kdp_hail_s,                                                          &
                     zdrtot_s,                                                                        &
                     zdr_rain_s, zdr_snow_s,                                                          &
                     zdr_grau_s, zdr_hail_s,                                                          &
                     zhtot_s,                                                                         &
                     zh_rain_s, zh_snow_s,                                                            &
                     zh_grau_s, zh_hail_s,                                                            &

                     !! C-band
                     aihtot_c,                                                                        &
                     aih_rain_c, aih_cnow_c,                                                          &
                     aih_grau_c, aih_hail_c,                                                          &
                     kdptot_c,                                                                        &
                     kdp_rain_c, kdp_cnow_c,                                                          &
                     kdp_grau_c, kdp_hail_c,                                                          &
                     zdrtot_c,                                                                        &
                     zdr_rain_c, zdr_cnow_c,                                                          &
                     zdr_grau_c, zdr_hail_c,                                                          &
                     zhtot_c,                                                                         &
                     zh_rain_c, zh_cnow_c,                                                            &
                     zh_grau_c, zh_hail_c,                                                            &

                     !! X-band
                     aihtot_x,                                                                        &
                     aih_rain_x, aih_xnow_x,                                                          &
                     aih_grau_x, aih_hail_x,                                                          &
                     kdptot_x,                                                                        &
                     kdp_rain_x, kdp_xnow_x,                                                          &
                     kdp_grau_x, kdp_hail_x,                                                          &
                     zdrtot_x,                                                                        &
                     zdr_rain_x, zdr_xnow_x,                                                          &
                     zdr_grau_x, zdr_hail_x,                                                          &
                     zhtot_x,                                                                         &
                     zh_rain_x, zh_xnow_x,                                                            &
                     zh_grau_x, zh_hail_x,                                                            &

                     ids,ide, jds,jde, kds,kde,                                                       &
                     ims,ime, jms,jme, kms,kme,                                                       &
                     ips,ipe, jps,jpe, kps,kpe,                                                       & ! patch  dims
                     i_start,i_end,j_start,j_end,kts,kte,num_tiles                                    &
                    )

!----------------------------------------------------------------------
! Main subroutine to compute CORDEX required variables

  USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
  USE module_configure 
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
! mu: perturbation dry air mass in column [Pa]
! mub: base dry air mass in column [Pa]
! sina: local sine of map rotation [1]
! cosa: local cosine of map rotation [1]
! ph: 3D gepotential perturbation [m2s-2]
! phb: 3D gepotential base state [m2s-2]
! rho: air density [kgm-3]

! nrad: amount of surface radar
! e_radlon: radar grid points on the x-axis
! e_radlat: radar grid points on the y-axis
! e_radhgt: radar grid points on the z-axis
! e_radnz: amount and height extent of radars (nradar)*(e_radhgt)

! radards: nominal horizontal resolution of radar values in m
! radinterp: method of interpolation to the radar space (1: lineal)
! missing_value: missing value
! gatelon: ground radar gate longitude [degree_east]
! gatelat: ground radar gate latitude [degree_north]
! gatehgt: ground radar gate height [m]
! gatehwght: ground radar horizontal weights for interpolation [1]
! inradar: pertinence of the grid cell to any radar (1: yes) [1]
! radarnum: pertinence of the grid cell inside a given radar area [1]

! timedbg: frecuency in minutes of time-debugging
! londbg: longitude in degrees of grid-point debugging (-999., not used)
! latdbg: latitude in degrees of grid-point debugging (-999., not used)

! Output
!-------
! vr: ground radar Doppler Velocity [ms-1]

!! S-band
! aihtot_s: S-band ground radar total attenuation [dB/km]
! aih_rain_s: S-band ground radar rain attenuation [dB/km]
! aih_snow_s: S-band ground radar snow attenuation [dB/km]
! aih_grau_s: S-band ground radar graupel attenuation [dB/km]
! aih_hail_s: S-band ground radar hail attenuation [dB/km]
! kdptot_s: S-band ground radar total specific differential phase [deg/km]
! kdp_rain_s: S-band ground radar rain specific differential phase [deg/km]
! kdp_snow_s: S-band ground radar snow specific differential phase [deg/km]
! kdp_grau_s: S-band ground radar graupel specific differential phase [deg/km]
! kdp_hail_s: S-band ground radar hail specific differential phase [deg/km]
! zdrtot_s: S-band ground radar total differential reflectivity [dBz]
! zdr_rain_s: S-band ground radar rain differential reflectivity [dBz]
! zdr_snow_s: S-band ground radar snow differential reflectivity [dBz]
! zdr_grau_s: S-band ground radar graupel differential reflectivity [dBz]
! zdr_hail_s: S-band ground radar hail differential reflectivity [dBz]
! zhtot_s: S-band ground radar total horizontal reflectivity [dBz]
! zh_rain_s: S-band ground radar rain horizontal reflectivity [dBz]
! zh_snow_s: S-band ground radar snow horizontal reflectivity [dBz]
! zh_grau_s: S-band ground radar graupel horizontal reflectivity [dBz]
! zh_hail_s: S-band ground radar hail horizontal reflectivity [dBz]

!! C-band
! aihtot_c: C-band ground radar total attenuation [dB/km]
! aih_rain_c: C-band ground radar rain attenuation [dB/km]
! aih_cnow_c: C-band ground radar snow attenuation [dB/km]
! aih_grau_c: C-band ground radar graupel attenuation [dB/km]
! aih_hail_c: C-band ground radar hail attenuation [dB/km]
! kdptot_c: C-band ground radar total specific differential phase [deg/km]
! kdp_rain_c: C-band ground radar rain specific differential phase [deg/km]
! kdp_cnow_c: C-band ground radar snow specific differential phase [deg/km]
! kdp_grau_c: C-band ground radar graupel specific differential phase [deg/km]
! kdp_hail_c: C-band ground radar hail specific differential phase [deg/km]
! zdrtot_c: C-band ground radar total differential reflectivity [dBz]
! zdr_rain_c: C-band ground radar rain differential reflectivity [dBz]
! zdr_cnow_c: C-band ground radar snow differential reflectivity [dBz]
! zdr_grau_c: C-band ground radar graupel differential reflectivity [dBz]
! zdr_hail_c: C-band ground radar hail differential reflectivity [dBz]
! zhtot_c: C-band ground radar total horizontal reflectivity [dBz]
! zh_rain_c: C-band ground radar rain horizontal reflectivity [dBz]
! zh_cnow_c: C-band ground radar snow horizontal reflectivity [dBz]
! zh_grau_c: C-band ground radar graupel horizontal reflectivity [dBz]
! zh_hail_c: C-band ground radar hail horizontal reflectivity [dBz]

!! X-band
! aihtot_x: X-band ground radar total attenuation [dB/km]
! aih_rain_x: X-band ground radar rain attenuation [dB/km]
! aih_xnow_x: X-band ground radar snow attenuation [dB/km]
! aih_grau_x: X-band ground radar graupel attenuation [dB/km]
! aih_hail_x: X-band ground radar hail attenuation [dB/km]
! kdptot_x: X-band ground radar total specific differential phase [deg/km]
! kdp_rain_x: X-band ground radar rain specific differential phase [deg/km]
! kdp_xnow_x: X-band ground radar snow specific differential phase [deg/km]
! kdp_grau_x: X-band ground radar graupel specific differential phase [deg/km]
! kdp_hail_x: X-band ground radar hail specific differential phase [deg/km]
! zdrtot_x: X-band ground radar total differential reflectivity [dBz]
! zdr_rain_x: X-band ground radar rain differential reflectivity [dBz]
! zdr_xnow_x: X-band ground radar snow differential reflectivity [dBz]
! zdr_grau_x: X-band ground radar graupel differential reflectivity [dBz]
! zdr_hail_x: X-band ground radar hail differential reflectivity [dBz]
! zhtot_x: X-band ground radar total horizontal reflectivity [dBz]
! zh_rain_x: X-band ground radar rain horizontal reflectivity [dBz]
! zh_xnow_x: X-band ground radar snow horizontal reflectivity [dBz]
! zh_grau_x: X-band ground radar graupel horizontal reflectivity [dBz]
! zh_hail_x: X-band ground radar hail horizontal reflectivity [dBz]

!======================================================================

  INTEGER,      INTENT(in)                               ::                                           &
                                            ids,ide, jds,jde, kds,kde,                                &
                                            ims,ime, jms,jme, kms,kme,                                &
                                            ips,ipe, jps,jpe, kps,kpe,                                &
                                            kts,kte,                                                  &
                                            num_tiles,                                                &
                                            n_moist


  INTEGER, DIMENSION(num_tiles), INTENT(in)              :: i_start, i_end, j_start, j_end
  INTEGER, INTENT(in)                                    :: m_qv, m_qc, m_qr, m_qs, m_qi, m_qg, m_qh

  REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)           :: lon, lat
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)           :: mu, mub, sina, cosa
  REAL, DIMENSION(kms:kme), INTENT(in)                   :: deta

  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(in)   :: qv, qc, qr, qs, qi, qg, qh, ph, phb, rho

  REAL, INTENT(in)                                       :: DT, XTIME
  REAL, INTENT(in)                                       :: curr_secs2
  LOGICAL, INTENT(in)                                    :: is_restart, is_output_time
  INTEGER, INTENT(in)                                    :: timedbg
  REAL, INTENT(in)                                       :: londbg, latdbg
  INTEGER, INTENT(in)                                    :: idbg, jdbg

  INTEGER, INTENT(inout)                                 :: nsteps
  
! Specific Radar variables
  INTEGER, INTENT(in)                                    :: nrad, e_radlon, e_radlat, e_radhgt, e_radnz
  INTEGER, INTENT(in)                                    :: radinterp
  REAL, INTENT(in)                                       :: radards, missing_value
  REAL, DIMENSION(nrad), INTENT(in)                      :: gatelon, gatelat, gatehgt
  REAL, DIMENSION(e_radnz,ims:ime,jms:jme), INTENT(in)   :: gatehwght
  INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(in)        :: inradar, radarnum
  
  REAL, DIMENSION

END MODULE module_simradar

