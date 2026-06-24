module types
  implicit none
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer, parameter     :: i1b=1 !integer 1 bit
  integer, parameter     :: i2b=2 !integer 2 bits
  integer, parameter     :: i4b=4 !integer 4 bits
  integer, parameter     :: i8b=8 !integer 8 bits
  integer, parameter     :: ik=i8b !main integer type
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer, parameter     :: sp=4  !single precision real
  integer, parameter     :: dp=8  !double precision real
  integer, parameter     :: qp=8 !"quad" accumulator kind (double here)
  integer, parameter     :: rk=dp !main real type
  integer, parameter     :: rks=dp!real type for arrays
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer, parameter     :: csp=4  !single precision complex
  integer, parameter     :: cdp=8  !double precision real
  integer, parameter     :: ck=cdp !main complex type
  integer, parameter     :: cks=rks!complex type for arrays
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real(rk), parameter    :: PI=3.14159265358979323846264338327950&
       &2884197169399375105820974944592307816406286208998628034825&
       &342117068_rk
  integer(ik) :: mpisize=1
end module types

module input_output
  use types
  use hdf5
  implicit none

contains

  subroutine read_hdf5_file(n1,n2,n3,u1,u2,u3,b1,b2,b3,filename,mhd)
    implicit none
    integer(ik), intent(IN)                           :: n1,n2,n3
    real(rks), dimension(1:n1,1:n2,1:n3), intent(OUT) :: u1,u2,u3,b1,b2,b3
    character(*), intent(IN)                          :: filename
    logical, intent(IN)                               :: mhd
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read velocity (/u) and, for MHD, magnetic (/b) vector fields from an
    ! ALIAKMON-style HDF5 file. Each dataset is stored as (3, n1, n2, n3).
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(hid_t)                             :: file_id
    integer                                    :: error
    real(rks), dimension(:,:,:,:), allocatable :: vec

    print '(3a)', 'Reading ', trim(filename), ' ...'

    call h5open_f(error)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

    allocate(vec(1:3,1:n1,1:n2,1:n3))

    call read_vector('/u', vec)
    u1 = vec(1,:,:,:)
    u2 = vec(2,:,:,:)
    u3 = vec(3,:,:,:)

    if(mhd) then
       call read_vector('/b', vec)
       b1 = vec(1,:,:,:)
       b2 = vec(2,:,:,:)
       b3 = vec(3,:,:,:)
    else
       b1 = 0.0_rks
       b2 = 0.0_rks
       b3 = 0.0_rks
    end if

    deallocate(vec)
    call h5fclose_f(file_id, error)
    call h5close_f(error)

    return

  contains

    subroutine read_vector(dataset_name, dataset)
      implicit none
      character(*), intent(IN)                              :: dataset_name
      real(rks), dimension(1:3,1:n1,1:n2,1:n3), intent(OUT) :: dataset
      integer(hid_t)                   :: dset_id
      integer(hsize_t), dimension(1:4) :: dims
      integer                          :: err

      dims = [3_hsize_t, int(n1,hsize_t), int(n2,hsize_t), int(n3,hsize_t)]
      call h5dopen_f(file_id, dataset_name, dset_id, err)
      if(rks==sp) then
         call h5dread_f(dset_id, H5T_NATIVE_REAL, dataset, dims, err)
      else
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dims, err)
      end if
      call h5dclose_f(dset_id, err)

      return
    end subroutine read_vector

  end subroutine read_hdf5_file

end module input_output

module random
  use types
  implicit none
  !random number generator
  integer(i4b) :: jseed,ifrst
  data JSEED,IFRST/123456789,0/



contains

  subroutine SRAND(ISEED)
    implicit none
    integer(ik) :: iseed
    !
    !  This subroutine sets the integer seed to be used with the
    !  companion RAND function to the value of ISEED.  A flag is
    !  set to indicate that the sequence of pseudo-random numbers
    !  for the specified seed should start from the beginning.
    !
    !
    JSEED = int(ISEED,i4b)
    IFRST = 0
    !
    return
  end subroutine SRAND

  real(rk) function RAND()
    implicit none
    !
    !  This function returns a pseudo-random number for each invocation.
    !  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
    !  standard number generator whose Pascal code appears in the article:
    !
    !     Park, Steven K. and Miller, Keith W., "Random Number Generators:
    !     Good Ones are Hard to Find", Communications of the ACM,
    !     October, 1988.
    !
    integer(i4b), parameter :: MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,&
         &MOMDMP=2836
    !
    integer(i4b) ::  HVLUE, LVLUE, TESTV
    integer(i4b), save :: NEXTN

    !
    if (IFRST .eq. 0) then
       NEXTN = JSEED
       IFRST = 1
    endif
    !
    HVLUE = NEXTN / MOBYMP
    LVLUE = mod(NEXTN, MOBYMP)
    TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
    if (TESTV .gt. 0) then
       NEXTN = TESTV
    else
       NEXTN = TESTV + MODLUS
    endif
    RAND = real(NEXTN)/real(MODLUS)
    !
    return
  end function RAND


end module random

module structure
  !use, intrinsic :: ieee_arithmetic
  use types
  implicit none
  logical :: MHD=.true.
contains


  subroutine structure_functions(n1,n2,n3,u1,u2,u3,b1,b2,b3,maxincr,maxord,maxpoints)
    use random
    implicit none
    integer(ik), intent(IN) :: n1,n2,n3,maxord,maxpoints
    integer(ik), intent(INOUT) :: maxincr
    integer :: nfilestrfun=0
    real(rks), dimension(1:n1,1:n2,1:n3), intent(IN) :: u1,u2,u3,b1,b2,b3
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik) :: i,j,k,n,ord
    real(rk) :: dx
    real(qp) :: tot
    real(rk), dimension(1:3) :: vel1,vel2,lvec,mag1,mag2,du,db
    real(rk) :: duL,duT,dbL,dbT,dzpL,dzmL,bL0
    real(rk) :: pv,pvt,pb,pbt,pzp,pzm
    real(rk) :: scale=1.0_rk
    character(256) :: fmt,fname
    integer :: lstrfun_file=21,trstrfun_file=22,psstrfun_file=25,&
         &pdstrfun_file=26,fourthirds_file=701,blstrfun_file=702,&
         &btrstrfun_file=703
    integer(ik) :: i1,i2,j1,j2,k1,k2
    real(rk) :: dist
    real(rk) :: ff
    real(rk), dimension(:),allocatable :: strfuncsv,strfuncsvt,strfuncsb,&
         &strfuncsbt,strfuncszp,strfuncszm
    integer(ik) :: nn1,nn2,nn3
    integer(i2b), dimension(:), allocatable :: iii,jjj,kkk
    integer(ik), dimension(:), allocatable  :: mm
    integer(i8b) :: niii,m,npoints
    integer(ik) :: inc,npointsi,irand,cputime

    nn1=n1
    nn2=n2
    nn3=n3
    !maxincr=n1/2

    write(fname,'(a,i5.5,a)') 'lstrfun.',nfilestrfun,'.dat'
    open(lstrfun_file,file=trim(fname),action='write',form='formatted')
    write(fname,'(a,i5.5,a)') 'trstrfun.',nfilestrfun,'.dat'
    open(trstrfun_file,file=trim(fname),action='write',form='formatted')

    if(MHD) then
       write(fname,'(a,i5.5,a)') 'fourthirds.',nfilestrfun,'.dat'
       open(fourthirds_file,file=trim(fname),action='write',&
            &form='formatted')
       write(fname,'(a,i5.5,a)') 'blstrfun.',nfilestrfun,'.dat'
       open(blstrfun_file,file=trim(fname),action='write',&
            &form='formatted')
       write(fname,'(a,i5.5,a)') 'btrstrfun.',nfilestrfun,'.dat'
       open(btrstrfun_file,file=trim(fname),action='write',&
            &form='formatted')
       write(fname,'(a,i5.5,a)') 'psstrfun.',nfilestrfun,'.dat'
       open(psstrfun_file,file=trim(fname),action='write',&
            &form='formatted')
       write(fname,'(a,i5.5,a)') 'pdstrfun.',nfilestrfun,'.dat'
       open(pdstrfun_file,file=trim(fname),action='write',&
            &form='formatted')
    end if
    write(fmt,'(a,i0,a)')  '(',maxord+1,'e35.14)'




    n=max(n1,n2,n3)
    tot=0

    allocate(strfuncsv(maxord),strfuncsvt(maxord))
    allocate(strfuncsb(maxord),strfuncsbt(maxord))
    allocate(strfuncszp(maxord),strfuncszm(maxord))


    niii=int(4./3.*3.14159*n**2)
    allocate(iii(niii))
    allocate(jjj(niii))
    allocate(kkk(niii))
    allocate(mm(niii))
    call system_clock(cputime)
    call srand(cputime)
    do inc=1,maxincr

       npointsi=0
       do i=-inc-1,inc+1
          do j=-inc-1,inc+1
             do k=-inc-1,inc+1
                if(inc**2 <= i**2+j**2+k**2.and.i**2+j**2+k**2<(inc+1)**2) then
                   npointsi=npointsi+1
                   if(npointsi>niii) then
                      print *, 'error: shell point count exceeds niii =', niii
                      stop 1
                   end if
                   iii(npointsi)=i
                   jjj(npointsi)=j
                   kkk(npointsi)=k
                end if
             end do
          end do
       end do

       mm=0
       if(npointsi>maxpoints) then
          npoints=maxpoints
          do i=1,maxpoints
             ! pick a uniform index in [1, npointsi], distinct from earlier picks
             do
                irand=int(rand()*npointsi)+1
                if(.not.any(irand==mm(1:i-1))) exit
             end do
             mm(i)=irand
          end do
       else
          npoints=npointsi
          do m=1,npointsi
             mm(m)=m
          end do
       end if
       
       tot=0.0
       strfuncsv=0.  ; strfuncsvt=0.
       strfuncsb=0.  ; strfuncsbt=0.
       strfuncszp=0. ; strfuncszm=0.
       ff=0.
       print *,'increment: ', inc
       !$omp parallel do private(i1,j1,m,lvec,dist,vel1,vel2,mag1,mag2,du,db,&
       !$omp& duL,duT,dbL,dbT,dzpL,dzmL,bL0,pv,pvt,pb,pbt,pzp,pzm,ord,i2,j2,k2) &
       !$omp& reduction(+:tot,strfuncsv,strfuncsvt,strfuncsb,strfuncsbt,&
       !$omp& strfuncszp,strfuncszm,ff)
       do k1=1,nn3
          do j1=1,nn2
             do i1=1,nn1
                do m=1,npoints
                   
                   i2=i1+iii(mm(m))
                   j2=j1+jjj(mm(m))
                   k2=k1+kkk(mm(m))

                   lvec(1)=i2-i1
                   lvec(2)=j2-j1
                   lvec(3)=k2-k1
                   dist=sqrt(dot_product(lvec(:),lvec(:)))

                   lvec(:)=lvec(:)/dist
                   tot=tot+1                      


                   vel1(1)=u1(i1,j1,k1)
                   vel1(2)=u2(i1,j1,k1)
                   vel1(3)=u3(i1,j1,k1)

                   i2=per(i2,nn1)
                   j2=per(j2,nn2)
                   k2=per(k2,nn3)
                   vel2(1)=u1(i2,j2,k2)
                   vel2(2)=u2(i2,j2,k2)
                   vel2(3)=u3(i2,j2,k2)
                   ! signed velocity increment, split into longitudinal and
                   ! transverse parts relative to the unit separation lvec
                   du(:)=vel2(:)-vel1(:)
                   duL=dot_product(du(:),lvec(:))
                   duT=sqrt(max(dot_product(du(:),du(:))-duL*duL,0.0_rk))

                   if(MHD) then
                      mag1(1)=b1(i1,j1,k1)
                      mag1(2)=b2(i1,j1,k1)
                      mag1(3)=b3(i1,j1,k1)
                      mag2(1)=b1(i2,j2,k2)
                      mag2(2)=b2(i2,j2,k2)
                      mag2(3)=b3(i2,j2,k2)
                      db(:)=mag2(:)-mag1(:)
                      dbL=dot_product(db(:),lvec(:))
                      dbT=sqrt(max(dot_product(db(:),db(:))-dbL*dbL,0.0_rk))
                      ! Elsasser z+/- = u +/- b : longitudinal increments
                      dzpL=duL+dbL
                      dzmL=duL-dbL
                   end if

                   ! accumulate order-p sums via running products
                   pv=1.0_rk ; pvt=1.0_rk
                   pb=1.0_rk ; pbt=1.0_rk ; pzp=1.0_rk ; pzm=1.0_rk
                   do ord=1,maxord
                      pv =pv *abs(duL) ; strfuncsv(ord) =strfuncsv(ord) +pv
                      pvt=pvt*duT      ; strfuncsvt(ord)=strfuncsvt(ord)+pvt
                      if(MHD) then
                         pb =pb *abs(dbL) ; strfuncsb(ord) =strfuncsb(ord) +pb
                         pbt=pbt*dbT      ; strfuncsbt(ord)=strfuncsbt(ord)+pbt
                         pzp=pzp*abs(dzpL); strfuncszp(ord)=strfuncszp(ord)+pzp
                         pzm=pzm*abs(dzmL); strfuncszm(ord)=strfuncszm(ord)+pzm
                      end if
                   end do

                   ! four-thirds (Politano-Pouquet) flux: one signed value per
                   ! point-pair, independent of the structure-function order.
                   if(MHD) then
                      bL0=dot_product(mag1(:),lvec(:))
                      ff=ff+scale*(duL**3 - 6.0_rk*bL0*bL0*duL)
                   end if

                end do
             end do
          end do
       end do
       !$omp end parallel do


       dx=2*PI/(n-1)

       write(lstrfun_file,fmt)  inc*dx,abs(strfuncsv(:))/tot
       write(trstrfun_file,fmt) inc*dx,abs(strfuncsvt(:))/tot

       if(MHD) then
          write(blstrfun_file,fmt)  inc*dx,abs(strfuncsb(:))/tot
          write(btrstrfun_file,fmt) inc*dx,abs(strfuncsbt(:))/tot
          write(psstrfun_file,fmt)  inc*dx,abs(strfuncszp(:))/tot
          write(pdstrfun_file,fmt)  inc*dx,abs(strfuncszm(:))/tot
          write(fourthirds_file,'(2e30.15)') inc*dx,ff/tot
       end if

    end do




    return


  end subroutine structure_functions

  pure elemental function per(ij,n)
    integer(ik) :: per
    integer(ik), intent(IN) :: ij,n
    if(ij<1) then
       per=ij+n
    else if(ij>n) then
       per=ij-n
    else
       per=ij
    end if

    return

  end function per

end module structure

program strfun
  use types
  use structure
  use input_output
  implicit none
  integer(ik) :: n,n1,n2,n3
  real(rks), dimension(:,:,:), allocatable :: u1,u2,u3,b1,b2,b3
  character(256) :: buf
  integer(ik) :: nord,nmaxincr,maxpoints
  integer(ik) :: nfile

  call get_command_argument(1,buf)
  read(buf,*) n
  call get_command_argument(2,buf)
  read(buf,*) mpisize
  call get_command_argument(3,buf)
  read(buf,*) nord
  call get_command_argument(4,buf)
  read(buf,*) maxpoints
  call get_command_argument(5,buf)
  read(buf,*) nfile
  
  n1=n
  n2=n
  n3=n
  allocate(u1(n,n,n3),u2(n,n,n3),u3(n,n,n3))
  allocate(b1(n,n,n3),b2(n,n,n3),b3(n,n,n3))
  MHD=.true.

  write(buf,'(a,i6.6,a)') 'output.', nfile, '.h5'
  call read_hdf5_file(n1,n2,n3,u1,u2,u3,b1,b2,b3,trim(buf),MHD)





  print *, 'Calculating structure functions...'

  nmaxincr=n/2
  call structure_functions(n,n,n3,u1,u2,u3,b1,b2,b3,nmaxincr,nord,maxpoints)

  print *, 'Done.'


  stop

end program strfun
