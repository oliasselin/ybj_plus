module IO_psi
  use parameters
  use mpi
  use fft
  use netcdf
  implicit none
  
!  integer :: idink,idoutk,idkx,idky,idkz,idkri,idktm,istat
!  integer :: idzxk,idzyk,idzzk,idttk,idtimek 

  integer :: idin,idout !Files ID
  integer :: idx,idy,idz,idt !Dimension IDs
  integer :: idpsi,idtime    !Variable IDs

CONTAINS
  
subroutine ncdumpout(psik,psir,ts)
  !! Create netcdf files and write the fields for restart (output)
  !! Only psir is dumped in a netcdf file
  implicit none

  double complex, dimension(iktx,ikty,n3h1) ::   psik
  double precision, dimension(n1d,n2d,n3h1)   :: psir

  double precision,    dimension(iktx,ikty,n3h1) :: psi_save
  real, dimension(n1,n2,n3h0) :: psi_clean             !suboptimal: de-haloed, de-extra x-dim values version of psi, real (not double)

  real,    intent(in)   :: ts
  integer, dimension(3) :: ncdims
  integer, dimension(3) :: nccount,ncstart

  !Move to r-space!
  psi_save = psik !Save to avoid inverse fft
  call fft_c2r(psik,psir,n3h1)

  !---------------------------------------------------!
  !Suboptimal: get rid of halos and extra x-dim values!

  do izh0=1,n3h0
     izh1=izh0+1
     do iy=1,n2
        do ix=1,n1

           psi_clean(ix,iy,izh0) = real(psir(ix,iy,izh1))

        end do
     end do
  end do
  !Recover the k-space psi without fft back
  psik=psi_save
  !---------------------------------------------------!

!---------------------------------------------------------------------
!     DEFINING THE OUTPUT NETCDF FILE 

!!! prep restart file (output) (define the variables and dims)
  if (mype.eq.0)  print*,'Yo! dumping in netcdf restart file'
  
  ! Create the netCDF file. The nf90_clobber parameter tells netCDF,
  call check( nf90_create("psi_fout.ncf",ior(NF90_NETCDF4,NF90_MPIIO),idout, comm = MPI_COMM_WORLD,info = MPI_INFO_NULL) )
  
  ! Define the dimension: kx, ky, kz. time. RI used as another dim to distinguish between real and imag parts
  call check( nf90_def_dim(idout, "x",int(n1,4),idx) )
  call check( nf90_def_dim(idout, "y",int(n2,4),idy) )
  call check( nf90_def_dim(idout, "z",int(n3,4),idz) )
  call check( nf90_def_dim(idout, "t",        1,idt) )
     
  ! ncdimsk array used to pass the dimension IDs to the variables
!  ncdimsk = (/idkx,idky,idkz,idkri/)
  ncdims = (/idx,idy,idz/)
  
  ! Define the variables which are the fields that going to be stored
  call check( nf90_def_var(idout,"psi" ,NF90_FLOAT,ncdims,idpsi ) )
  call check( nf90_def_var(idout,"time",NF90_FLOAT,idt   ,idtime) )
     
  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(idout) )

!---------------------------------------------------------------------
  !     WRITING VARIABLES IN NETCDF FILE
  
  ! how many to count in each dimension when writing files
  nccount(1) = n1
  nccount(2) = n2
  nccount(3) = n3h0

  ! where to start on the output file
  ncstart(1) = 1
  ncstart(2) = 1
  
!  call check( nf90_put_var(idout, idpsi, real(psir)) )

  ! wrtie time (time is written only one time so just root process is used
  if (mype.eq.0) then 
     call check( nf90_put_var(idout, idtime, ts) )
  endif

  ! in the z-direction mype=0 is in the first place
  ncstart(3) = mype*n3h0+1
  call check( nf90_put_var(idout, idpsi, psi_clean,ncstart,nccount))    

  
!---------------------------------------------------------------------
!     CLOSING the  NETCDF FILE 
  call check (nf90_close(idout))



contains
  subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped!!!"
    end if
  end subroutine check
  
end subroutine ncdumpout



subroutine ncreadin(psik,psir,ts)
! Read the netcdf data from an input file (psi.in.ncf)                                                                
  implicit none

  double complex, dimension(iktx,ikty,n3h1) ::   psik
  double precision, dimension(n1d,n2d,n3h1)   :: psir

  real, dimension(n1,n2,n3h0) :: psi_clean   !Directly from the read file

  integer :: n1in,n2in,n3in         !dimensions read in the file. Must match n1,n2,n3 for now

!  integer, dimension(3) :: ncdims
  integer, dimension(3) :: nccount,ncstart


!  complex, intent(out),    dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
!  real, dimension(iktx,ikty,iktzp) :: wr,wi
  real,    intent(out) :: ts

!  integer, dimension(5) :: ncstartr,ncstarti,nccount

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access                                                            
  print*,'Yo! reading from netcdf restart file, BTW im mype',mype
  call check( nf90_open('psi_pout.ncf', NF90_NOWRITE,idin) )
  print*,'Yo! just opened netcdf restart file, BTW im mype',mype

  ! Get the dimensions IDs based on their name                                                                                   
  call check( nf90_inq_dimid(idin, "x", idx) )
  call check( nf90_inq_dimid(idin, "y", idy) )
  call check( nf90_inq_dimid(idin, "z", idz) )
!  call check( nf90_inq_dimid(idin, "t", idt) )

  ! Get the dimension length and check if the grid resolution matches                                                   
  call check( nf90_inquire_dimension(idin,idx,len=n1in))
  call check( nf90_inquire_dimension(idin,idy,len=n2in))
  call check( nf90_inquire_dimension(idin,idz,len=n3in))

  if (n1in.ne.n1 .or. n2in.ne.n2 .or. n3in.ne.n3) then
    print*,'Sorry, do not know how to change resolution.'
    stop
  endif

  ! Get the variables IDs                                                                                                                                                      
!  call check( nf90_inq_varid(idin,"time",idtime))
  call check( nf90_inq_varid(idin, "psi",idpsi) )

  ! prepare for reading variables                                                                                                                                              
  ncstart    = 1
  nccount(1) = n1
  nccount(2) = n2
  nccount(3) = n3h0
  ncstart(3)= int(mype*n3h0+1)

  !Get variables
!  call check(  nf90_get_var(idin,idtime,ts))
  call check( nf90_get_var(idin,idpsi,psi_clean,ncstart,nccount))
  
  !Back to psik!
  psir = 0.D0
  do izh0=1,n3h0
     izh1=izh0+1
     do iy=1,n2
        do ix=1,n1

          psir(ix,iy,izh1) =  psi_clean(ix,iy,izh0) 

        end do
     end do
  end do
  call fft_r2c(psir,psik,n3h1)
  call generate_halo_q(psik)
  !------------!

  call check (nf90_close(idin))

contains
  subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status)),'mype=',mype
       stop "Stopped!!!"
    end if
  end subroutine check

end subroutine ncreadin




end module IO_psi
