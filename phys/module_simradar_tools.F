MODULE module_simradar_tools
! Module wiuth the tools for simrad.AR
!   L. Fita, CIMA-IFAECI, November 2025

  USE module_wrf_error

!!!!!!! Subroutines / functions
! distance_between_points: Subroutine tro compute the coordinates at a given distance and angle from a 
!   given point
! DIRCT1: Solution of the geodetic direct problem after T. Vincenty
! inverse: Subroutine to compute distance between 2 points
! is_inside: Subroutine to check if a tile gets a radar
! pertinence_polar: Subroutine to compute the 2-dimensional pertinence of a lonlat polar radial space 
!   (r) into another regular lonlat (WRF)
! position_from_point: Subroutine tro compute the coordinates at a given distance and angle from a 
!   given point
! read_lut: Subroutine to read the T-MATRIX look-up-table (LUT) from an specific microphysics from a 
!   netCDF file
! radar_inside: Subroutine to determine if a radar fells inside the tile
! real_distances_ellipsoide: Subroutine to compute the diferential of coordinattes of a distance from 
!   a given point on an ellipsoide


! basic: (from PyNCplot)
! ErrMsg: Subroutine to stop execution and provide an error message
! ErrMsg[I/R/S]Avail: Subroutine to stop execution by non exsitent Integer/Real/String value from 
!   available ones and provide an error message
! [I/R]toS: Function to transform an Integer/Real to String
! set_dbg: Function to setup debug value
! Sto[I/R/L]: Function to transform a String to an Integer/Real/Boolean

! netCDF: (from PyNCplot)
! get_dimlength_ncid: Subroutine to get the length of a dimension from an existing to an open netCDF
! get_Ndims_ncid: Subroutine to get the amount of dimensions from an open netCDF file
! get_Nvars_ncid: Subroutine to get the amount of variables from an open netCDF file
! get_vars_ncid: Subroutine to provide all the ids and names of the variables existing to an open 
!   netCDF
! get_varD3D_ncid: Subroutine to get a 3D double Real variable from a netCDF file unit passing 
!   starting and ending
! get_varS0D_ncid: Subroutine to get a 0D String variable from a netCDF file unit passing starting 
!   and ending
! isin_ncunit: Function to tell if a given variable is inside a netcdf file unit
! isin_dim_ncid: Function to tell if a given dimension is inside a netcdf file unit
! open_NCfile: Subroutine to open a basic netCDF file

  IMPLICIT NONE
  
  ! Integer 4 bytes 
  INTEGER, PARAMETER                                     :: ii4 = SELECTED_INT_KIND(4)
  ! REal 16 bytes 
  INTEGER, PARAMETER                                     :: rr8 = SELECTED_REAL_KIND(8)
  
  REAL(kind=rr8), PARAMETER                              :: pi = 4.d0*DATAN(1.d0)
  REAL(kind=rr8), PARAMETER                              :: rad = 180.d0/pi

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
  
  CHARACTER(len=33)                                      :: errormsg = 'ERROR -- error -- ERROR -- error'
  
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
    
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_rain_vec
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_snow_vec
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_grau_vec
    REAL(kind=rr8), DIMENSION(:), ALLOCATABLE            :: q_ice_vec
    
  END TYPE lut_vars

  CONTAINS
     
  SUBROUTINE CALL read_lut(filen, lutgrid, idbg)
  ! Subroutine to read the T-MATRIX look-up-table (LUT) from an specific microphysics from a netCDF 
  !   file
  
    IMPLICIT NONE
  
    CHARACTER(len=S256), INTENT(in)                      :: filen
    LOGICAL, OPTIONAL                                    :: idbg
    TYPE(lut_vars), INTENT(out)                          :: lutgrid
    
    ! Local
    INTEGER                                              :: ncid, ierr
    CHARACTER(len=S16)                                   :: S16a
    CHARACTER(len=S32)                                   :: S32a
    CHARACTER(len=S256)                                  :: S256a
    LOGICAL                                              :: dbg
    CHARACTER(len=Sl)                                    :: fname
    
    !!!!!!! Variables
    ! lutfilen: name of the netCDF file with the values
    ! lutgrid: derived type where to kept the variables
    
    fname = 'read_lut'
    
    dbg = set_dbg(idbg)
    
    CALL open_NCfile(filen, ncid, dbg)

    CALL get_Ndims_ncid(lutid, Nlutdims)
    CALL get_Nvars_ncid(lutid, Nlutvars)
    IF (dbg) PRINT *,'  ' // TRIM(fname) // ": getting from file '" // TRIM(lutfilen) //              &
      "' Ndims: ", Nlutdims, ' Nvars: ', Nlutvars

    S32a = 'WRF_mp_physics'
    CALL get_varS0D_ncid(ncid, 256, S32a, S256a, dbg)
    S16a = TRIM(S256a)
    lutgrid%mp_scheme = S16a
    
    ! LUT space
    S32a = 'elevation'
    CALL get_dimlength_ncid(ncid, S32a, lutgrid%Nlutelevs)
    S32a = 'radar_band'
    CALL get_dimlength_ncid(ncid, S32a, lutgrid%Nbands)
    S32a = 'qxbin'
    CALL get_dimlength_ncid(ncid, S32a, lutgrid%Nlutqxbins)
    
    !! Variables
    !!! RAIN
    IF (ALLOCATED(lutgrid%Zh_RAIN)) DEALLOCATE(lutgrid%Zh_RAIN)
    ALLOCATE(lutgrid%Zh_RAIN(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_RAIN', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zh_RAIN'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zh_RAIN, dbg)
    IF (ALLOCATED(lutgrid%Zdr_RAIN)) DEALLOCATE(lutgrid%Zdr_RAIN)
    ALLOCATE(lutgrid%Zdr_RAIN(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_RAIN', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zdr_RAIN'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zdr_RAIN, dbg)
    IF (ALLOCATED(lutgrid%LDR_RAIN)) DEALLOCATE(lutgrid%LDR_RAIN)
    ALLOCATE(lutgrid%LDR_RAIN(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_RAIN', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'LDR_RAIN'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%LDR_RAIN, dbg)
    IF (ALLOCATED(lutgrid%Aih_RAIN)) DEALLOCATE(lutgrid%Aih_RAIN)
    ALLOCATE(lutgrid%Aih_RAIN(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_RAIN', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aih_RAIN'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aih_RAIN, dbg)
    IF (ALLOCATED(lutgrid%Aiv_RAIN)) DEALLOCATE(lutgrid%Aiv_RAIN)
    ALLOCATE(lutgrid%Aiv_RAIN(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_RAIN', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aiv_RAIN'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aiv_RAIN, dbg)
    IF (ALLOCATED(lutgrid%KDP_RAIN)) DEALLOCATE(lutgrid%KDP_RAIN)
    ALLOCATE(lutgrid%KDP_RAIN(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_RAIN', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'KDP_RAIN'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%KDP_RAIN, dbg)

    !!! SNOW
    IF (ALLOCATED(lutgrid%Zh_SNOW)) DEALLOCATE(lutgrid%Zh_SNOW)
    ALLOCATE(lutgrid%Zh_SNOW(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_SNOW', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zh_SNOW'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zh_SNOW, dbg)
    IF (ALLOCATED(lutgrid%Zdr_SNOW)) DEALLOCATE(lutgrid%Zdr_SNOW)
    ALLOCATE(lutgrid%Zdr_SNOW(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_SNOW', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zdr_SNOW'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zdr_SNOW, dbg)
    IF (ALLOCATED(lutgrid%LDR_SNOW)) DEALLOCATE(lutgrid%LDR_SNOW)
    ALLOCATE(lutgrid%LDR_SNOW(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_SNOW', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'LDR_SNOW'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%LDR_SNOW, dbg)
    IF (ALLOCATED(lutgrid%Aih_SNOW)) DEALLOCATE(lutgrid%Aih_SNOW)
    ALLOCATE(lutgrid%Aih_SNOW(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_SNOW', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aih_SNOW'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aih_SNOW, dbg)
    IF (ALLOCATED(lutgrid%Aiv_SNOW)) DEALLOCATE(lutgrid%Aiv_SNOW)
    ALLOCATE(lutgrid%Aiv_SNOW(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_SNOW', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aiv_SNOW'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aiv_SNOW, dbg)
    IF (ALLOCATED(lutgrid%KDP_SNOW)) DEALLOCATE(lutgrid%KDP_SNOW)
    ALLOCATE(lutgrid%KDP_SNOW(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_SNOW', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'KDP_SNOW'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%KDP_SNOW, dbg)

    !!! GRAU
    IF (ALLOCATED(lutgrid%Zh_GRAU)) DEALLOCATE(lutgrid%Zh_GRAU)
    ALLOCATE(lutgrid%Zh_GRAU(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_GRAU', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zh_GRAU'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zh_GRAU, dbg)
    IF (ALLOCATED(lutgrid%Zdr_GRAU)) DEALLOCATE(lutgrid%Zdr_GRAU)
    ALLOCATE(lutgrid%Zdr_GRAU(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_GRAU', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zdr_GRAU'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zdr_GRAU, dbg)
    IF (ALLOCATED(lutgrid%LDR_GRAU)) DEALLOCATE(lutgrid%LDR_GRAU)
    ALLOCATE(lutgrid%LDR_GRAU(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_GRAU', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'LDR_GRAU'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%LDR_GRAU, dbg)
    IF (ALLOCATED(lutgrid%Aih_GRAU)) DEALLOCATE(lutgrid%Aih_GRAU)
    ALLOCATE(lutgrid%Aih_GRAU(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_GRAU', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aih_GRAU'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aih_GRAU, dbg)
    IF (ALLOCATED(lutgrid%Aiv_GRAU)) DEALLOCATE(lutgrid%Aiv_GRAU)
    ALLOCATE(lutgrid%Aiv_GRAU(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_GRAU', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aiv_GRAU'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aiv_GRAU, dbg)
    IF (ALLOCATED(lutgrid%KDP_GRAU)) DEALLOCATE(lutgrid%KDP_GRAU)
    ALLOCATE(lutgrid%KDP_GRAU(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_GRAU', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'KDP_GRAU'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%KDP_GRAU, dbg)

    !!! ICE
    IF (ALLOCATED(lutgrid%Zh_ICE)) DEALLOCATE(lutgrid%Zh_ICE)
    ALLOCATE(lutgrid%Zh_ICE(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zh_ICE', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zh_ICE'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zh_ICE, dbg)
    IF (ALLOCATED(lutgrid%Zdr_ICE)) DEALLOCATE(lutgrid%Zdr_ICE)
    ALLOCATE(lutgrid%Zdr_ICE(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Zdr_ICE', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Zdr_ICE'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Zdr_ICE, dbg)
    IF (ALLOCATED(lutgrid%LDR_ICE)) DEALLOCATE(lutgrid%LDR_ICE)
    ALLOCATE(lutgrid%LDR_ICE(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'LDR_ICE', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'LDR_ICE'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%LDR_ICE, dbg)
    IF (ALLOCATED(lutgrid%Aih_ICE)) DEALLOCATE(lutgrid%Aih_ICE)
    ALLOCATE(lutgrid%Aih_ICE(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aih_ICE', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aih_ICE'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aih_ICE, dbg)
    IF (ALLOCATED(lutgrid%Aiv_ICE)) DEALLOCATE(lutgrid%Aiv_ICE)
    ALLOCATE(lutgrid%Aiv_ICE(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'Aiv_ICE', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'Aiv_ICE'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%Aiv_ICE, dbg)
    IF (ALLOCATED(lutgrid%KDP_ICE)) DEALLOCATE(lutgrid%KDP_ICE)
    ALLOCATE(lutgrid%KDP_ICE(Nbands,Nlutelevs,Nlutqxbins), STAT=ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'KDP_ICE', 3, (/Nbands,Nlutelevs,Nlutqxbins/), fname, ierr)
    S32a = 'KDP_ICE'
    CALL get_varD3D_ncid(ncid, 1, 1, 1, Nbands, Nlutelevs, Nlutqxbins, S32a, lutgrid%KDP_ICE, dbg)

!    S32a = 'q_rain_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_rain_min, dbg)
!    S32a = 'q_rain_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_rain_max, dbg)
!    S32a = 'qn_rain_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_rain_min, dbg)
!    S32a = 'qn_rain_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_rain_max, dbg)
!    S32a = 'q_snow_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_snow_min, dbg)
!    S32a = 'q_snow_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_snow_max, dbg)
!    S32a = 'qn_snow_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_snow_min, dbg)
!    S32a = 'qn_snow_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_snow_max, dbg)
!    S32a = 'q_grau_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_grau_min, dbg)
!    S32a = 'q_grau_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_grau_max, dbg)
!    S32a = 'qn_grau_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_grau_min, dbg)
!    S32a = 'qn_grau_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_grau_max, dbg)
!    S32a = 'q_ice_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_ice_min, dbg)
!    S32a = 'q_ice_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_ice_max, dbg)
!    S32a = 'qn_ice_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_ice_min, dbg)
!    S32a = 'qn_ice_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%qn_ice_max, dbg)
!    S32a = 'q_clou_min'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_clou_min, dbg)
!    S32a = 'q_clou_max'
!    CALL get_varR0D_ncid(ncid, S32a, lutgrid%q_clou_max, dbg)
!    S32a = 'dim_ice'
!    CALL get_varN0D_ncid(ncid, S32a, lutgrid%dim_ice, dbg)
!    S32a = 'dimv'
!    CALL get_varN0D_ncid(ncid, S32a, lutgrid%dimv, dbg)

    IF (ALLOCATED(lutgrid%elevations)) DEALLOCATE(lutgrid%elevations)
    ALLOCATE(lutgrid%elevations(lutgrid%Nlutelevs), ierr)
    IF (ierr /= 0) CALL ErrMsgAlloc('', 'elevations', 1, (/lutgrid%Nlutelevs/), fname, ierr)

    S32a = 'elevation'
    CALL get_varR1D_ncid(ncid, 1, Nlutelevs, S32a, lutgrid%elevations, dbg)

    CALL close_NCfile(filen, ncid, dbg)

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
    CHARACTER(len=50)                                      :: fname

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
    REAL(kind=rr8)                                       :: dsig, ep2, eps
    REAL(kind=rr8)                                       :: prev, tem1, tem2, temp, test
    REAL(kind=rr8)                                       :: sinal, sinlam, sinsig
    REAL(kind=rr8)                                       :: beta1, sinu1, cosu1
    REAL(kind=rr8)                                       :: beta2, sinu2, cosu2
    LOGICAL                                              :: converged, keeplooping
    CHARACTER(len=Sl)                                    :: fname
       
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
        sinal = cosu1*cosu2*sinlam/dsign(eps,sinsig)
      ELSE
        sinal = cosu1*cosu2*sinlam/sinsig
      ENDIF

      cosal2 = -sinal**2+1.d0
      !*** v18   (rapp part II)
      !***** costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig
      IF(dabs(cosal2).lt.eps) THEN
        !*** avoid div 0
        costm = -2.d0*(sinu1*sinu2/dsign(eps,cosal2))+cossig
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

        IF (kindv .eq. 1) THEN
          lam = l+d*sinal
          IF (DABS(lam-test) .lt. tol) THEN
            converged = .TRUE.
            EXIT
            
          END IF

          IF (DABS(lam) .gt. dpict) THEN
          
            keeplooping = .FALSE.

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
              costm= -2.d0*(sinu1*sinu2/dsign(eps,cosal2))+cossig
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
            costm = -2.d0*(sinu1*sinu2/dsign(eps,cosal2)) + cossig
          ELSE
            costm = -2.d0*(sinu1*sinu2/cosal2)+cossig
          ENDIF

          costm2=costm*costm

        ENDIF
    
      END DO

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
    CHARACTER(len=50)                                    :: fname

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
    CHARACTER(len=50)                                    :: fname
    
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
    DO ibnd=1, 4
      CALL is_inside(ii, ie, ji, je, bnds_lon(ii:ie,ji:je,ibnd), bnds_lat(ii:ie,ji:je,ibnd), rlon,    &
        rlat, rdist, inside, dbg)  
    END DO

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
      IF (dbg) PRINT *,'    ' // TRIM(fname) // ": found ", SUM(radnum, mask=mask), " inside"
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
    CHARACTER(len=50)                                    :: fname
    
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

    IF (dbg) PRINT *, '   ' //TRIM(fname)//' for ims,ime;jms,jme ',ims,' , ',ime,' ; ',jms,' , ',jme, &
      ' min,max lon ', MINVAL(lon), ',', MAXVAL(lon), ' min,max lat ', MINVAL(lat), ',', MAXVAL(lat), &
      ' mindist: ', mindist, ' inside ', inside

    IF (ALLOCATED(distborder)) DEALLOCATE(distborder)
    
    RETURN

  END SUBROUTINE is_inside

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
    CHARACTER(len=50)                                    :: fname
    
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
      "ERROR -- error -- ERROR -- error:   " // TRIM(funcn) // ": " // TRIM(msg) //                   &
      ' not ready !!  Available ones: ' // (TRIM(ItoS(avail)) // ', ', i=1, Navail)
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
      "ERROR -- error -- ERROR -- error:   " // TRIM(funcn) // ": " // TRIM(msg) //                   &
      ' not ready !!  Available ones: ' // "'" // (TRIM(avail)) // "', ", i=1, Navail)
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
    emsg = 'ERORR -- error -- ERROR -- error'

    IF (errN /= 0) THEN
      WRITE (diag_message , * ) &
      "ERROR -- error -- ERROR -- error:   " // TRIM(funcn) // ": " // TRIM(msg) //                   &
      ' Problem allocating : ' // "'" // TRIM(arrayn)) // "' (", (TRIM(ItoS(shapev(i))) // ', ',i=1,  &
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
    CHARACTER(LEN=50)                                    :: emsg

    !!!!!!! Variables
    ! msg: message to print with the error
    ! funcn: name of the funciton, section to localize the error
    ! errN: number of the error provided for a given FORTRAN function

    emsg = 'ERORR -- error -- ERROR -- error'

    IF (errN /= 0) THEN
      PRINT *,TRiM(emsg)
      PRINT *,'  ' // TRIM(funcn) // ': ' // TRIM(msg) // ' !!'
      PRINT *,'    error number:', errN
      STOP

      WRITE (diag_message , * ) &
      "ERROR -- error -- ERROR -- error:   " // TRIM(funcn) // ": " // TRIM(msg) // ' !!'
      CALL wrf_error_fatal ( diag_message )
    END IF

    RETURN

  END SUBROUTINE ErrMsg

  CHARACTER(LEN=50) FUNCTION RtoS(Rval)
  ! RtoS: Function to transform a real to String

    IMPLICIT NONE

    REAL(r_k), INTENT(IN)                                :: Rval

    WRITE(RtoS,'(F21.10)')Rval

  END FUNCTION RtoS

  CHARACTER(LEN=50) FUNCTION ItoS(Ival)
  ! ItoS: Function to transform an integer to String

    IMPLICIT NONE

    INTEGER, INTENT(IN)                                  :: Ival

    ! Local
    CHARACTER(LEN=50)                                    :: itoS0

    WRITE(ItoS0,'(I50)')Ival
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
    CHARACTER(len=Sl)                                    :: fname
    
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
    CHARACTER(len=Sl)                                    :: fname
    
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
    CHARACTER(len=Sl)                                    :: fname
    
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
    CHARACTER(len=50)                                    :: fname

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
    CHARACTER(len=Sl)                                    :: fname

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
    CHARACTER(len=Sl)                                    :: fname
    
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
    CHARACTER(len=Sl)                                    :: fname

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
    CHARACTER(len=Sl)                                    :: fname

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
  
END MODULE module_simradar_tools
