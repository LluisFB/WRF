MODULE module_simradar_tools
! Module wiuth the tools for simrad.AR
!   L. Fita, CIMA-IFAECI, November 2025

!!!!!!! Subroutines / functions
! is_inside: Subroutine to check if a tile gets a radar
! DIRCT1: Solution of the geodetic direct problem after T. Vincenty
! position_from_point: Subroutine tro compute the coordinates at a given distance and angle from a 
!   given point
! radar_inside: Subroutine to determine if a radar fells inside the tile
! real_distances: Subroutine to compute the diferential of coordinattes of a distance from a given 
!   point

  IMPLICIT NONE
  
  ! Integer 4 bytes 
  INTEGER, PARAMETER                                     :: ii4 = SELECTED_INT_KIND(4)
  ! REal 16 bytes 
  INTEGER, PARAMETER                                     :: rr8 = SELECTED_REAL_KIND(8)
  
  REAL(kind=rr8), PARAMETER                              :: pi = 4.d0*DATAN(1.d0)
  REAL(kind=rr8), PARAMETER                              :: rad = 180.d0/pi
  
  CHARACTER(len=33)                                      :: errormsg = 'ERROR -- error -- ERROR -- error'
  
  CONTAINS  
  
  SUBROUTINE position_from_point(lon_pt, lat_pt, a, f, elips, dist, angle, lon_dist, lat_dist)
  ! Subroutine tro compute the coordinates at a given distance and angle from a given point
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
    CHARACTER(len=50)                                    :: fname

    !!!!!!! Variables
    ! lon_pt, lat_pt: coordinates of the point
    ! a: semimajor axis equatorial (in meters)
    ! f: flattening
    ! elips: name of the ellipsis being used
    ! dist: distance [m]
    ! angle: angle
    fname = 'position_from_point'
      
    glat1 = lat_pt / rad
    glon1 = lon_pt / rad
    faz = angle / rad
  
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
   
    lat_dist = glat2*rad
    lon_dist = glon2*rad
  
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
    
    CHARACTER(len=50)                                    :: fname

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
    BAZ=DATAN2(SA,BAZ)+PI

    RETURN

  END SUBROUTINE DIRCT1
  
  SUBROUTINE radar_inside (ims, ime, jms, jme, lon, lat, bnds_lon, bnds_lat, rlon, rlat, rdist,       &
    missing_value, gatehdist, inradar, radnum, inside)
  ! Subroutine to determine if a radar fells inside the tile
  
    IMPLICIT NONE
  
    INTEGER, INTENT(in)                                  :: ims, ime, jms, jme
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: lon, lat
    REAL, DIMENSION(ims:ime,jms:jme,4), INTENT(in)       :: bnds_lon, bnds_lat
    REAL, INTENT(in)                                     :: rlon, rlat, rdist, missing_value
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(out)        :: gatehdist
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(out)     :: inradar, radnum
    LOGICAL, INTENT(out)                                 :: inside
  
    ! Local
    INTEGER                                              :: i, j
    INTEGER                                              :: ii, ie, ji, je
    INTEGER                                              :: ibnd
    CHARACTER(len=50)                                    :: fname
    
    !!!!!!! Variables
    ! lon, lat: 2D array longitude, latitude matrices at the center 
    ! bnds_lon, bnds_lat: 2D array longitude, latitude matrices of the vertices
    ! rlonv, rlat: location of the radar
    ! rdist: maximum distance to the radar
    ! gatehdist: distance to the radar in the model grid
    ! inradar: whether the model grid lays inside the measurement area of any radar
    ! radnum: pertinence of the model grid to a given radar
    fname = 'radar_inside'
  
    ! Is this radar inside the range of the tile?
    inside = .FALSE.
    
    gatehdist = missing_value
    inradar = -1
    radnum = -1
      
    ! Do we need to be efficient?
    !! S border
    ii = ims
    ie = ime
    ji = jms
    je = jms
    CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside)

    ! Going to the borders
    DO ibnd=1, 4
      CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,    &
        rlat, rdist, inside)  
    END DO
    
    !! W border
    IF (.NOT.inside) THEN
      ii = ims
      ie = ims
      ji = jms
      je = jme
      CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside)

      ! Going to the borders
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside)  
      END DO
    END IF
    
    !! N border
    IF (.NOT.inside) THEN
      ii = ims
      ie = ime
      ji = jme
      je = jme
      CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside)

      ! Going to the borders
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside)  
      END DO
    END IF
        
    !! E border
    IF (.NOT.inside) THEN
      ii = ime
      ie = ime
      ji = jms
      je = jme
      CALL is_inside(ii, ie, ji, je, lon(ii:ie,ji:je), lat(ii:ie,ji:je), rlon, rlat, rdist, inside)

      ! Going to the borders
      DO ibnd=1, 4
        CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,  &
          rlat, rdist, inside)  
      END DO
    END IF

    ! Filling other values with a single loop
    IF (inside) THEN
      DO i=ims, ime
        DO j=jms, jme     
          CALL real_distances(ims, ime, jms, jme, lon, lat, rlon, rlat, gatehdist(i,j))
          IF (gatehdist(i,j) <= rdist) THEN
            inradar(i,j) = 1
            radnum(i,j) = 1        
          ELSE
            gatehdist(i,j) = missing_value
          END IF
        END DO
      END DO
    END IF

    RETURN

  END SUBROUTINE radar_inside
  
  SUBROUTINE is_inside(ime, ims, jms, jme, lon, lat, rlon, rlat, rdist, inside)
  ! Subroutine to check if a tile gets a radar
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: ims, ime, jms, jme
    REAL, DIMENSION(ims:ime,jms:jme), INTENT(in)         :: lon, lat
    REAL, INTENT(in)                                     :: rlon, rlat, rdist
    LOGICAL, INTENT(inout)                               :: inside
  
    ! Local
    REAL                                                 :: mindist
    REAL, DIMENSION(:,:), ALLOCATABLE                    :: distborder  
    CHARACTER(len=50)                                    :: fname
    
    !!!!!!! Variables
    ! lon, lat: 2D array longitude, latitude matrices at the center 
    ! rlonv, rlat: location of the radar
    ! rdist: maximum distance to the radar
    ! gatehdist: distance to the radar in the model grid
    fname = 'is_inside'
      
    IF (ALLOCATED(distborder)) DEALLOCATE(distborder)
    ALLOCATE(distborder(ims:ime,jms:jme))
    CALL real_distances(ims, ime, jms, jme,                                                           &
      lon(ims:ime,jms:jme), lat(ims:ime,jms:jme),                                                     &
      rlon, rlat,                                                                                     &
      distborder)
    mindist = MINVAL(distborder)
    IF (mindist < rdist) inside = .TRUE.

    IF (ALLOCATED(distborder)) DEALLOCATE(distborder)
    
    RETURN

  END SUBROUTINE is_inside
  
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
    REAL, PARAMETER                                      :: R = 6371.009 
    REAL, DIMENSION(ims:ime,jms:jme)                     :: dlat, dlon, latm
    CHARACTER(len=50)                                    :: fname
    
    !!!!!!! Variables
    ! lon, lat: 2D array longitude, latitude matrices
    ! lonv, latv: point to compute the distance
    fname = 'real_distance'
  
    dist = 0.
    
    dlat = lat - latv
    dlon = lon - lonv
    latm = (lat + latv) / 2.
    
    dist = 2.*R * SQRT( (sin(dlat /2.)*cos(dlon/2.))**2 + (cos(latm)*sin(dlon/2.))**2 )
   
    RETURN
 
  END SUBROUTINE real_distances

END MODULE module_simradar_tools
