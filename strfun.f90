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
end module types

module input_output
  use types
  use hdf5
  implicit none

contains

subroutine check_h5(error, context)
   implicit none
   integer, intent(IN) :: error
   character(*), intent(IN) :: context
   if(error < 0) then
      print '(a,i0,2a)', 'FATAL: HDF5 error ', error, ' during ', trim(context)
      stop 1
   end if
end subroutine check_h5

subroutine probe_hdf5_fields(filename, mhd, rad)
   ! Open the input file read-only and report which optional fields are
   ! present: magnetic field (/b) enables MHD, and the radiation pair
   ! (/G and /q) enables RAD.
   implicit none
   character(*), intent(IN) :: filename
   logical, intent(OUT)     :: mhd, rad
   integer(hid_t) :: file_id
   integer        :: error
   logical        :: has_b, has_g, has_q

   call h5open_f(error)
   call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
   call check_h5(error, 'opening '//trim(filename))

   call h5lexists_f(file_id, '/b', has_b, error)
   call h5lexists_f(file_id, '/G', has_g, error)
   call h5lexists_f(file_id, '/q', has_q, error)

   call h5fclose_f(file_id, error)
   call h5close_f(error)

   mhd = has_b
   rad = has_g .and. has_q
   if(has_g .neqv. has_q) then
      print '(a)', 'Warning: only one of /G, /q present; radiation disabled.'
   end if
   print '(a,l1,a,l1)', 'Detected fields -> MHD: ', mhd, '  RAD: ', rad

   return
end subroutine probe_hdf5_fields

subroutine write_strfun_to_hdf5(hdf5_fname, strfun_name, maxincr, maxord, strfun_data, tot)
   use types
   implicit none
   character(*), intent(IN) :: hdf5_fname
   character(*), intent(IN) :: strfun_name
   integer(ik), intent(IN) :: maxincr, maxord
   real(rk), dimension(1:maxincr, 1:maxord), intent(IN) :: strfun_data
   real(rk), dimension(1:maxincr), intent(IN) :: tot
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(rk), dimension(1:maxincr) :: out
   character(64) :: name_buf
   integer(ik) :: ord

   do ord=1,maxord
      out(:) = strfun_data(:,ord)/tot(:)
      write(name_buf,'(a,i0)') strfun_name, ord
      call add_hdf5_1d(hdf5_fname, name_buf, out)
   end do
   
   end subroutine write_strfun_to_hdf5

   subroutine add_hdf5_1d(filename, dsetname, data_array)
    use hdf5
    implicit none

    ! arguments
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dsetname
    real(kind=8), dimension(:), intent(in), contiguous :: data_array

    ! parameters & dimensions
    integer, parameter :: rank = 1
    integer(hsize_t), dimension(1) :: dims
    
    ! identifiers & error flag
    integer(hid_t) :: file_id, dataspace_id, dataset_id
    integer :: error
    logical, save :: first_time = .true.
    logical :: file_exists
    logical :: dset_exists

    ! get the size of the passed array dynamically
    dims(1) = size(data_array)

    ! 1. initialize hdf5 library
    call h5open_f(error)

    ! 2. first call: create/truncate the file (warn if it already exists);
    !    later calls: open read-write and append, replacing existing datasets
    if(first_time) then
       inquire(file=filename, exist=file_exists)
       if(file_exists) then
          print *, 'Warning: file ', trim(filename), ' already exists. Overwriting.'
       end if
       call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)
       call check_h5(error, 'creating '//trim(filename))
       first_time = .false.
    else
        call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)
        call check_h5(error, 'opening '//trim(filename))
        
        ! check if a dataset with this name already exists
        call h5lexists_f(file_id, dsetname, dset_exists, error)
        if (dset_exists) then
            ! delete old link to prevent "name already exists" crash
            call h5ldelete_f(file_id, dsetname, error)
        end if
    end if

    ! 3. create the 1d dataspace based on array size
    call h5screate_simple_f(rank, dims, dataspace_id, error)

    ! 4. create the dataset (using native double type)
    call h5dcreate_f(file_id, dsetname, h5t_native_double, dataspace_id, &
                     dataset_id, error)
    call check_h5(error, 'creating dataset '//trim(dsetname))

    ! 5. write the array data (using native double type)
    call h5dwrite_f(dataset_id, h5t_native_double, data_array, dims, error)
    call check_h5(error, 'writing dataset '//trim(dsetname))

    ! 6. clean up identifiers
    call h5dclose_f(dataset_id, error)
    call h5sclose_f(dataspace_id, error)
    call h5fclose_f(file_id, error)

    ! 7. close hdf5 library
    call h5close_f(error)
end subroutine add_hdf5_1d


  subroutine read_hdf5_file(n1,n2,n3,u1,u2,u3,b1,b2,b3,g,q1,q2,q3,filename,mhd,rad)
    implicit none
    integer(ik), intent(IN)                           :: n1,n2,n3
    real(rks), dimension(1:n1,1:n2,1:n3), intent(OUT) :: u1,u2,u3
    real(rks), dimension(:,:,:), intent(OUT)          :: b1,b2,b3
    real(rks), dimension(:,:,:), intent(OUT)          :: g,q1,q2,q3
    character(*), intent(IN)                          :: filename
    logical, intent(IN)                               :: mhd,rad
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read velocity (/u) and, optionally, magnetic (/b) and radiative-flux
    ! (/q) vector fields plus the scalar incident radiation (/G) from an
    ! ALIAKMON-style HDF5 file. Vector datasets are stored as (3,n1,n2,n3),
    ! scalar datasets as (n1,n2,n3).
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(hid_t)                             :: file_id
    integer                                    :: error
    real(rks), dimension(:,:,:,:), allocatable :: vec

    print '(3a)', 'Reading ', trim(filename), ' ...'

    call h5open_f(error)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
    call check_h5(error, 'opening '//trim(filename))

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

    if(rad) then
       call read_scalar('/G', g)
       call read_vector('/q', vec)
       q1 = vec(1,:,:,:)
       q2 = vec(2,:,:,:)
       q3 = vec(3,:,:,:)
    else
       g  = 0.0_rks
       q1 = 0.0_rks
       q2 = 0.0_rks
       q3 = 0.0_rks
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
      call check_h5(err, 'opening dataset '//trim(dataset_name))
      if(rks==sp) then
         call h5dread_f(dset_id, H5T_NATIVE_REAL, dataset, dims, err)
      else
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dims, err)
      end if
      call check_h5(err, 'reading dataset '//trim(dataset_name))
      call h5dclose_f(dset_id, err)

      return
    end subroutine read_vector

    subroutine read_scalar(dataset_name, dataset)
      implicit none
      character(*), intent(IN)                          :: dataset_name
      real(rks), dimension(1:n1,1:n2,1:n3), intent(OUT) :: dataset
      integer(hid_t)                   :: dset_id
      integer(hsize_t), dimension(1:3) :: dims
      integer                          :: err

      dims = [int(n1,hsize_t), int(n2,hsize_t), int(n3,hsize_t)]
      call h5dopen_f(file_id, dataset_name, dset_id, err)
      call check_h5(err, 'opening dataset '//trim(dataset_name))
      if(rks==sp) then
         call h5dread_f(dset_id, H5T_NATIVE_REAL, dataset, dims, err)
      else
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dims, err)
      end if
      call check_h5(err, 'reading dataset '//trim(dataset_name))
      call h5dclose_f(dset_id, err)

      return
    end subroutine read_scalar

  end subroutine read_hdf5_file

end module input_output

module random
  use types
  implicit none
  !random number generator
  integer(i4b) :: jseed,ifrst
  data JSEED,IFRST/123456789,0/



contains

  subroutine pm_srand(ISEED)
    implicit none
    integer(ik) :: iseed
    !
    !  This subroutine sets the integer seed to be used with the
    !  companion pm_rand function to the value of ISEED.  A flag is
    !  set to indicate that the sequence of pseudo-random numbers
    !  for the specified seed should start from the beginning.
    !
    !
    JSEED = int(ISEED,i4b)
    IFRST = 0
    !
    return
  end subroutine pm_srand

  real(rk) function pm_rand()
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
    pm_rand = real(NEXTN)/real(MODLUS)
    !
    return
  end function pm_rand


end module random

module structure
  !use, intrinsic :: ieee_arithmetic
  use types
  implicit none
  logical :: MHD=.true.
  logical :: RAD=.false.
contains


  subroutine structure_functions(n1,n2,n3,u1,u2,u3,b1,b2,b3,g,q1,q2,q3,maxincr,maxord,maxpoints,nfile)
    use types
    use random
    use input_output
    implicit none
    integer(ik), intent(IN) :: n1,n2,n3,maxord,maxpoints,nfile
    integer(ik), intent(IN) :: maxincr
    real(rks), dimension(1:n1,1:n2,1:n3), intent(IN) :: u1,u2,u3
    real(rks), dimension(:,:,:), intent(IN)          :: b1,b2,b3
    real(rks), dimension(:,:,:), intent(IN)          :: g,q1,q2,q3
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik) :: i,j,k,n,ord
    real(qp) :: tot
    real(rk), dimension(1:3) :: vel1,vel2,lvec,mag1,mag2,du,db,qr1,qr2,dq
    real(rk) :: duL,duT,dbL,dbT,dzpL,dzmL,bL0
    real(rk) :: dG,dqL,dqT
    real(rk) :: pv,pvt,pb,pbt,pzp,pzm
    real(rk) :: pg,pql,pqt
    real(rk) :: scale=1.0_rk
    integer(ik) :: i1,i2,j1,j2,k1,k2
    real(rk) :: dist
    real(rk) :: ff
    real(rk), dimension(:,:),allocatable :: strfuncsv,strfuncsvt,strfuncsb,&
         &strfuncsbt,strfuncszp,strfuncszm
    real(rk), dimension(:,:),allocatable :: strfuncsg,strfuncsql,strfuncsqt
    real(rk), dimension(:), allocatable :: dx
    real(rk), dimension(:), allocatable :: totinc
    real(rk), dimension(:), allocatable :: ffinc
    integer(ik) :: nn1,nn2,nn3
    integer(i2b), dimension(:), allocatable :: iii,jjj,kkk
    integer(ik), dimension(:), allocatable  :: mm
    integer(i8b) :: niii,m,npoints
    integer(ik) :: inc,npointsi,irand,cputime
    character(256) :: hdf5_fname

    nn1=n1
    nn2=n2
    nn3=n3
    n=max(n1,n2,n3)
    tot=0

    allocate(dx(maxincr))
    allocate(totinc(maxincr))
    allocate(ffinc(maxincr))
    allocate(strfuncsv(maxincr,maxord),strfuncsvt(maxincr,maxord))
    allocate(strfuncsb(maxincr,maxord),strfuncsbt(maxincr,maxord))
    allocate(strfuncszp(maxincr,maxord),strfuncszm(maxincr,maxord))
    allocate(strfuncsg(maxincr,maxord),strfuncsql(maxincr,maxord),&
         &strfuncsqt(maxincr,maxord))


    ! generous upper bound on the lattice-point count in the largest shell
    ! (inc = maxincr); i8b avoids overflow and the runtime check below is the
    ! real backstop should this ever be too small
    niii=int(4.0_rk/3.0_rk*PI*(real(maxincr+1,rk)**3-real(maxincr,rk)**3),i8b) &
         +int(8.0_rk*PI*real(maxincr,rk),i8b)+64_i8b
    allocate(iii(niii))
    allocate(jjj(niii))
    allocate(kkk(niii))
    allocate(mm(niii))
    call system_clock(cputime)
    call pm_srand(cputime)
    do inc=1,maxincr
       dx(inc) = real(inc, rk) * 2.0 * PI/ real(n, rk)
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
                   iii(npointsi)=int(i,i2b)
                   jjj(npointsi)=int(j,i2b)
                   kkk(npointsi)=int(k,i2b)
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
                irand=int(pm_rand()*npointsi)+1
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
       strfuncsv(inc, :)=0.  ; strfuncsvt(inc, :)=0.
       strfuncsb(inc, :)=0.  ; strfuncsbt(inc, :)=0.
       strfuncszp(inc, :)=0. ; strfuncszm(inc, :)=0.
       strfuncsg(inc, :)=0.  ; strfuncsql(inc, :)=0. ; strfuncsqt(inc, :)=0.
       ff=0.
       print *,'increment: ', inc
       !$omp parallel do private(i1,j1,m,lvec,dist,vel1,vel2,mag1,mag2,du,db,&
       !$omp& duL,duT,dbL,dbT,dzpL,dzmL,bL0,pv,pvt,pb,pbt,pzp,pzm,ord,i2,j2,k2,&
       !$omp& qr1,qr2,dq,dG,dqL,dqT,pg,pql,pqt) &
       !$omp& reduction(+:tot,strfuncsv,strfuncsvt,strfuncsb,strfuncsbt,&
       !$omp& strfuncszp,strfuncszm,strfuncsg,strfuncsql,strfuncsqt,ff)
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

                   if(RAD) then
                      ! scalar incident-radiation increment
                      dG=g(i2,j2,k2)-g(i1,j1,k1)
                      ! radiative-flux vector increment, split into
                      ! longitudinal and transverse parts relative to lvec
                      qr1(1)=q1(i1,j1,k1)
                      qr1(2)=q2(i1,j1,k1)
                      qr1(3)=q3(i1,j1,k1)
                      qr2(1)=q1(i2,j2,k2)
                      qr2(2)=q2(i2,j2,k2)
                      qr2(3)=q3(i2,j2,k2)
                      dq(:)=qr2(:)-qr1(:)
                      dqL=dot_product(dq(:),lvec(:))
                      dqT=sqrt(max(dot_product(dq(:),dq(:))-dqL*dqL,0.0_rk))
                   end if

                   ! accumulate order-p sums via running products
                   pv=1.0_rk ; pvt=1.0_rk
                   pb=1.0_rk ; pbt=1.0_rk ; pzp=1.0_rk ; pzm=1.0_rk
                   pg=1.0_rk ; pql=1.0_rk ; pqt=1.0_rk
                   do ord=1,maxord
                      pv =pv *abs(duL) ; strfuncsv(inc,ord) =strfuncsv(inc,ord) +pv
                      pvt=pvt*duT      ; strfuncsvt(inc,ord)=strfuncsvt(inc,ord)+pvt
                      if(MHD) then
                         pb =pb *abs(dbL) ; strfuncsb(inc,ord) =strfuncsb(inc,ord) +pb
                         pbt=pbt*dbT      ; strfuncsbt(inc,ord)=strfuncsbt(inc,ord)+pbt
                         pzp=pzp*abs(dzpL); strfuncszp(inc,ord)=strfuncszp(inc,ord)+pzp
                         pzm=pzm*abs(dzmL); strfuncszm(inc,ord)=strfuncszm(inc,ord)+pzm
                      end if
                      if(RAD) then
                         pg =pg *abs(dG)  ; strfuncsg(inc,ord) =strfuncsg(inc,ord) +pg
                         pql=pql*abs(dqL) ; strfuncsql(inc,ord)=strfuncsql(inc,ord)+pql
                         pqt=pqt*dqT      ; strfuncsqt(inc,ord)=strfuncsqt(inc,ord)+pqt
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

       ! pair count for this increment: used to normalize this row only
       totinc(inc) = tot
       ffinc(inc) = ff

    end do

    write(hdf5_fname,'(a,i6.6,a)') 'strfun.',nfile,'.h5'
    call add_hdf5_1d(hdf5_fname,'Dx',dx)
    call write_strfun_to_hdf5(hdf5_fname, 'Du_l', maxincr, maxord, strfuncsv, totinc)
    call write_strfun_to_hdf5(hdf5_fname, 'Du_t', maxincr, maxord, strfuncsvt, totinc)
    if(MHD) then
       call write_strfun_to_hdf5(hdf5_fname, 'Db_l', maxincr, maxord, strfuncsb, totinc)
       call write_strfun_to_hdf5(hdf5_fname, 'Db_t', maxincr, maxord, strfuncsbt, totinc)
       call write_strfun_to_hdf5(hdf5_fname, 'Dzp', maxincr, maxord, strfuncszp, totinc)
       call write_strfun_to_hdf5(hdf5_fname, 'Dzm', maxincr, maxord, strfuncszm, totinc)
       ! Politano-Pouquet 4/3 flux: one signed value per increment
       call add_hdf5_1d(hdf5_fname, 'Fourthirds', ffinc/totinc)
    end if
    if(RAD) then
       call write_strfun_to_hdf5(hdf5_fname, 'DG',   maxincr, maxord, strfuncsg,  totinc)
       call write_strfun_to_hdf5(hdf5_fname, 'Dq_l', maxincr, maxord, strfuncsql, totinc)
       call write_strfun_to_hdf5(hdf5_fname, 'Dq_t', maxincr, maxord, strfuncsqt, totinc)
    end if

    return

  end subroutine structure_functions

  pure elemental function per(ij,n)
    ! periodic wrap of a 1-based index into [1, n] for any integer ij
    integer(ik) :: per
    integer(ik), intent(IN) :: ij,n
    per = modulo(ij-1, n) + 1
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
  real(rks), dimension(:,:,:), allocatable :: g,q1,q2,q3
  character(256) :: buf, fname
  integer(ik) :: nord,nmaxincr,maxpoints,nfile
  integer :: ip1,ip2

  if(command_argument_count() /= 5) then
   n = 128
   nord = 3
   maxpoints = 100
   nmaxincr = 64
   fname = 'output.999999.h5'
  else
   call get_command_argument(1,buf)
   read(buf,*) n
   call get_command_argument(2,buf)
   read(buf,*) nord
   call get_command_argument(3,buf)
   read(buf,*) maxpoints
   call get_command_argument(4,buf)
   read(buf,*) nmaxincr
   ! read the filename directly (list-directed input would stop at a '/')
   call get_command_argument(5,fname)
  end if

  ! input files follow the output.<nfile>.h5 convention; parse the snapshot
  ! number nfile out of the filename so the output can be tagged with it
  ip2 = index(trim(fname), '.h5', back=.true.)
  ip1 = index(fname(1:max(ip2-1,1)), '.', back=.true.)
  if(ip1>0 .and. ip2>ip1+1) then
     read(fname(ip1+1:ip2-1),*) nfile
  else
     print *, 'Warning: could not parse nfile from ', trim(fname), '; using 0'
     nfile = 0
  end if

  n1=n
  n2=n
  n3=n
  ! detect which optional fields (/b, /G, /q) are present, then allocate
  ! and read only what is actually in the file
  call probe_hdf5_fields(trim(fname), MHD, RAD)
  allocate(u1(n1,n2,n3),u2(n1,n2,n3),u3(n1,n2,n3))
  if(MHD) then
     allocate(b1(n1,n2,n3),b2(n1,n2,n3),b3(n1,n2,n3))
  else
     ! no magnetic field needed; allocate empty so the dummy args stay valid
     allocate(b1(0,0,0),b2(0,0,0),b3(0,0,0))
  end if
  if(RAD) then
     allocate(g(n1,n2,n3),q1(n1,n2,n3),q2(n1,n2,n3),q3(n1,n2,n3))
  else
     ! no radiation fields needed; allocate empty so the dummy args stay valid
     allocate(g(0,0,0),q1(0,0,0),q2(0,0,0),q3(0,0,0))
  end if

  call read_hdf5_file(n1,n2,n3,u1,u2,u3,b1,b2,b3,g,q1,q2,q3,trim(fname),MHD,RAD)

  call structure_functions(n,n,n3,u1,u2,u3,b1,b2,b3,g,q1,q2,q3,nmaxincr,nord,maxpoints,nfile)

  print *, 'Done.'


  stop

end program strfun
