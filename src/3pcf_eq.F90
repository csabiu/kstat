program ThreePCF
  use kdtree2_module
  use kdtree2_precision_module
  implicit none
#ifdef MPI
	include 'mpif.h'
#endif

  real(kdkind), dimension(:,:), allocatable :: my_array, my_array2,bins,tbins
  real(kdkind), allocatable ::  v1(:),v2(:),p1(:),p2(:),p3(:)

  type(kdtree2), pointer :: tree
  real(kdkind), allocatable :: Zddd(:,:,:),Zddr(:,:,:),Zdrr(:,:,:),Zrrr(:,:,:),crap(:,:,:), wgt1(:)

  real(kdkind) :: ndd,ndr,nrr,r1min,r1max,r2min,r2max
  integer :: k, i, j, l,d,chunk,nbins,ind
  character :: file1*2000, file2*2000, ranfile*2000, outfile*2000
  type(kdtree2_result),allocatable :: resultsb(:), results2(:)
  integer   ::  Ndata,Nrand,nn1,nn2,nnpairs,mm1,mm2,nmu,mu1,mu2
  real(kdkind)  :: aux,odr,theta,odt,mid1,mid2, theta2,theta3, pivalby2
  real(kdkind) :: rmin, rmax, dr, vec_dist
  logical :: wgt, logbins, RSD, loadran,saveran
  integer :: myid , ntasks , ierr
  real(kdkind), parameter :: pival=3.14159265358979323846d0
#ifdef MPI 
  integer , dimension( MPI_STATUS_SIZE ) :: status
#endif
  integer , parameter :: master=0,msgtag1=11,msgtag2=12,msgtag3=13,msgtag4=14

! ! I n i t i a l i z a t i o n of the MPI environment
myid=0
ntasks=1

#ifdef MPI 
 call MPI_INIT( ierr )
 call MPI_COMM_RANK( MPI_COMM_WORLD , myid, ierr )
 call MPI_COMM_SIZE( MPI_COMM_WORLD , ntasks , ierr )
#endif

!--- set some default parameters ---
  call default_params()

  
  if(myid==master) print*,'reading options'
  call parseOptions()

 allocate(bins(nbins+2,2))
 allocate(tbins(nbins,2))

if(logbins) then
  dr=(log10(rmax)-log10(rmin))/nbins
else
  dr=(rmax-rmin)/nbins
endif

  do i=1,nbins
    if(logbins) then
      bins(i,1)=10.0**(log10(rmin) + dr*(float(i)-1.0))
      bins(i,2)=10.0**(log10(rmin) + (dr*i))
    else
      bins(i,1)=rmin + dr*(float(i)-1.0)
      bins(i,2)=rmin + dr*i
    endif
    if(myid==0) print*,bins(i,:)
  enddo

   if(myid==master)  print*,'allocated bins'

  call count_files()

  Allocate(my_array(d,Ndata+Nrand))
  Allocate(my_array2(d,Ndata+Nrand))

  if ( wgt ) then
    Allocate(wgt1(Ndata+Nrand))
  endif

  allocate(v1(d))
  allocate(v2(d))
  allocate(p1(d))
  allocate(p2(d))
  allocate(p3(d))

  call read_files()

print*,' building tree on node', myid
tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)     ! this is how you create a tree.
print*,' built tree on node', myid

  call allocate_arrays()


!odr=1./(pival/(2.*nmu))
odr=1./(pival/(1.*nmu))
pivalby2=pival/2.

print*,'beginning data loop on thread:',myid

#ifdef MPI
chunk=floor(float(Ndata)/ntasks)+1
do i=max(myid*chunk,1),min(((myid+1)*chunk)-1,Ndata),1
#else
do i=1,Ndata,1
#endif

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins,2)*bins(nbins,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

    do ind=nbins,1,-1       
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(ind,1)*bins(ind,1))      

      do j=nn1+1,nn2,1
        v1=my_array(:,resultsb(j)%idx)

        do k=nn1+1,nn2,1
          v2=my_array(:,resultsb(k)%idx)        

          call triplet_bin( )

          if(vec_dist<bins(ind,2) .and. vec_dist>bins(ind,1)) then

          if(resultsb(j)%idx <=Ndata) then 
          
            if(resultsb(k)%idx <=Ndata) then

              if ( wgt ) then
                Zddd(ind,mu1,mu2)=Zddd(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zddd(ind,mu1,mu2)=Zddd(ind,mu1,mu2)+1.d0
              endif

            else

              if ( wgt ) then
                Zddr(ind,mu1,mu2)=Zddr(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zddr(ind,mu1,mu2)=Zddr(ind,mu1,mu2)+1.d0
              endif
            endif

          else 

            if(resultsb(k)%idx <=Ndata) then

              if ( wgt ) then
                Zddr(ind,mu1,mu2)=Zddr(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zddr(ind,mu1,mu2)=Zddr(ind,mu1,mu2)+1.d0
              endif

            else

              if ( wgt ) then
                Zdrr(ind,mu1,mu2)=Zdrr(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zdrr(ind,mu1,mu2)=Zdrr(ind,mu1,mu2)+1.d0
              endif

            endif
          endif
        endif
        
      enddo
   enddo
nn2=nn1
enddo

enddo

print*,'finished data loop in thread:', myid

if(loadran) then
  if(myid==master) print*,'Random counts exist, skipping some calculation'
  goto 57
endif

!RR loop
print*,'beginning random loop on thread:',myid

#ifdef MPI
chunk=floor(float(Nrand)/ntasks)+1
do i=Ndata+max(myid*chunk,1),Ndata+min(((myid+1)*chunk)-1,Nrand),1
#else
do i=Ndata+1,Ndata+Nrand,1
#endif


   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins,2)*bins(nbins,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

    do ind=nbins,1,-1       
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(ind,1)*bins(ind,1))      

      do j=nn1+1,nn2,1
        v1=my_array(:,resultsb(j)%idx)

        do k=nn1+1,nn2,1
          v2=my_array(:,resultsb(k)%idx)        
          
          call triplet_bin( )

          if(vec_dist<bins(ind,2) .and. vec_dist>bins(ind,1)) then

         if(resultsb(j)%idx <=Ndata) then 
          if(resultsb(k)%idx <=Ndata) then 
            if ( wgt ) then
               Zddr(ind,mu1,mu2)=Zddr(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zddr(ind,mu1,mu2)=Zddr(ind,mu1,mu2)+1.d0
            endif
         else
            if ( wgt ) then
               Zdrr(ind,mu1,mu2)=Zdrr(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zdrr(ind,mu1,mu2)=Zdrr(ind,mu1,mu2)+1.d0
            endif
          endif
            
         else 

          if(resultsb(k)%idx <=Ndata) then 
           if ( wgt ) then
               Zdrr(ind,mu1,mu2)=Zdrr(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zdrr(ind,mu1,mu2)=Zdrr(ind,mu1,mu2)+1.d0
            endif
         else
            if ( wgt ) then
               Zrrr(ind,mu1,mu2)=Zrrr(ind,mu1,mu2)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zrrr(ind,mu1,mu2)=Zrrr(ind,mu1,mu2)+1.d0
            endif
         endif
       endif
       endif
      enddo
   enddo
nn2=nn1
enddo

enddo
print*,'finished RR loop on thread:', myid

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

57 continue

#ifdef MPI
call mpi_collect()
#endif

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

if(myid==master) then

call normalise_counts()
call load_save_randoms()

open(11,file=outfile,status='unknown')
do i=1,nbins
   do j=1,nmu
     do k=1,nmu   
        print*,bins(i,1),bins(i,2), (j-1)/odr, j/odr, (k-1)/odr, k/odr,(Zddd(i,j,k)-Zddr(i,j,k)+Zdrr(i,j,k)-Zrrr(i,j,k))/Zrrr(i,j,k)
        write(11,'(12(e14.7,1x))') bins(i,1),bins(i,2), (j-1)/odr,j/odr, (k-1)/odr, k/odr,Zddd(i,j,k), Zddr(i,j,k), Zdrr(i,j,k), &
        Zrrr(i,j,k), &
        (Zddd(i,j,k)-Zddr(i,j,k)+Zdrr(i,j,k)-Zrrr(i,j,k))/Zrrr(i,j,k), &
        (Zddd(i,j,k)/Zrrr(i,j,k)) - 1.0
    enddo
   enddo  
enddo
close(11)

endif

! deallocate the memory for the data.
call kdtree2_destroy(tree)  

call deallocate_arrays()

if(myid==master) then
  print*, "Exit... stage left!"
  stop
else
  stop
endif


contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parseOptions ()
     USE extension, ONLY : getArgument, nArguments
      implicit none

      logical lexist
      integer(kind=4)    :: i,n
      character(len=6)   :: opt
      character(len=2000) :: arg

      n = nArguments()

      if (n < 8 .and. myid==master) then
        print*,' '
        print*,'Not enough input parameters. Please read the following help info'
        print*,' '
        call print_help
        stop
      end if

      do i = 1,n,2
         call getArgument(i,opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
!            call fatal_error
         end if
         call getArgument(i+1,arg)
         select case (opt)
         case ('-gal')
            file1 = trim(arg)
         case ('-ran')
            file2 = trim(arg)
         case ('-out')
            outfile = trim(arg)
         case ('-rmin')
            read (arg,*) rmin
         case ('-rmax')
            read (arg,*) rmax
         case ('-nbins')
            read (arg,*) nbins
         case ('-nmu')
            read (arg,*) nmu
         case ('-RSD')
            if(trim(arg)=='.true.') then
              RSD=.true.
              if(myid==master) print*,'Anisotropic analysis'
            elseif(trim(arg)=='.false.') then
              RSD=.false.
            else
                stop "Error using -RSD flag. Either .true. or .false."
            endif
         case ('-wgt')
            if(trim(arg)=='.true.') then
              wgt=.true.
              if(myid==master) print*,'Using weighted points.'
            elseif(trim(arg)=='.false.') then
              wgt=.false.
            else
                stop "Error using -wgt flag. Either .true. or .false."
            endif
        case('-RR')        
            ranfile = trim(arg)
            inquire (file=ranfile,exist=lexist)
            if(lexist) then
               write(*,*) " Random count file exists!"
               loadran=.true.
            else
               write(*,*) " Random count file does not exist. It will be created!"
               saveran=.true.
            end if     

         case ('-help')
            call print_help
      stop

         case default
            print '("unknown option ",a6," ignored")', opt
            stop
         end select

      end do

   end subroutine parseOptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_help
print*,'PURPOSE: Code for calculating the 3PCF of a 3D point set. '
print*,'     '
print*,'     '
print*,'CALLING SEQUENCE:'
print*,'      3pcf [-gal gal_file] [-ran ran_file] [-out out_file] [-wgt WGT]'
print*,'           [-rmin Rmin] [-rmax Rmax] [-nbins Nbins] '
print*,'           [-RSD RSD] [-nmu Nmu] [-RR RR_Counts]'
print*,' '
print*,'      eg: 3pcf -gal dr9.gal -ran dr9.ran -r1min 10.0 -r1max 12.0 -r2min 5.0 -r2max 6.0 -nbins 10'
print*,' '
print*,'INPUTS:'
print*,'       gal_file/ranfile = strings containing the name of a 3/4 column point file'
print*,'                 lines 1-3: XYZ comoving, optional 4: weight  '
print*,'   '
print*,'       out_file = filename for output correlation func. '
print*,'   '
print*,'       Rmin = Min. radial seperation'
print*,'   '
print*,'       Rmax = Max. radial seperation'
print*,' '
print*,'OPTIONAL:'
print*,'       Nbins = Number of radial bins to calculate'
print*,' '
print*,'       WGT = logical value to define if weights are to be included '
print*,'              In this case, the input gal/ran files should be 4 columns '
print*,' '
print*,'       RSD = logical value to request anisotropic analysis.'
print*,'              In this case the number of angular bins Nmu should be set'
print*,' '
print*,'    RR_Counts = string containing the name of a file regarding RR data counts.  '
print*,'                If file exists, then the RR counting will be skipped (saveing CPU time)'
print*,'                If it does not exist, the file will be created for future use.'
print*,'                Consider using this option if you are doing Monte Carlo runs with same random catalogue.' 
print*,''
   end subroutine print_help
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine count_files ()
    implicit none

  Ndata=0
  open(11,file=file1,status='unknown')
12 continue
  read(11,*,err=12,end=13)aux
  Ndata=Ndata+1
  goto 12
13 close(11)
  if(myid==0) print*,'Preparing to read ',Ndata, 'data points'

  Nrand=0
  open(11,file=file2,status='unknown')
21 continue
  read(11,*,err=21,end=22)aux
  Nrand=Nrand+1
  goto 21
22 close(11)
  if(myid==0) print*,'Preparing to read ',Nrand, 'random points'

end subroutine count_files

  subroutine read_files ()
    implicit none

  if(myid==0)  print*, 'opening ',file1
  open(11,file=file1,status='unknown')
  do i=1,Ndata
    if ( wgt ) then
     read(11,*,end=14)my_array(1:3,i),wgt1(i)
    else
     read(11,*,end=14)my_array(1:3,i)
   endif

  enddo
14 close(11)
  if(myid==0)  print*,'Finished reading data file 1'  

  if(myid==0)  print*, 'opening ',file2
  open(11,file=file2,status='unknown')
  do i=Ndata+1,Ndata+Nrand
    if ( wgt ) then
     read(11,*,end=23)my_array(1:3,i),wgt1(i)
   else
     read(11,*,end=23)my_array(1:3,i)
   endif

  enddo
23 close(11)
  if(myid==0)  print*,'Finished reading data file 2'  

if(wgt) then
!    wgt1=wgt1/100000.0000000000
!    wgt1(1:Ndata)=wgt1(1:Ndata)/sum(wgt1(1:Ndata))
!    wgt1(Ndata+1:Ndata+Nrand)=wgt1(Ndata+1:Ndata+Nrand)/sum(wgt1(Ndata+1:Ndata+Nrand))
endif
  

end subroutine read_files

subroutine allocate_arrays ()
  implicit none
  allocate(resultsb(Ndata+Nrand))
  allocate(results2(Ndata+Nrand))

  allocate(Zddd(nbins,nmu,nmu))
  allocate(Zddr(nbins,nmu,nmu))
  allocate(Zdrr(nbins,nmu,nmu))
  allocate(Zrrr(nbins,nmu,nmu))
  allocate(crap(nbins,nmu,nmu))
  
  crap=0.0d0
  Zddd=0.0d0
  Zddr=0.0d0
  Zdrr=0.0d0
  Zrrr=0.0d0

end subroutine allocate_arrays

subroutine mpi_collect()
#ifdef MPI
  if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zddd, crap, nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zddd=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zddr, crap, nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zddr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdrr, crap, nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zrrr, crap, nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zrrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

end subroutine

subroutine normalise_counts()
  if ( wgt ) then
   ndd=sum(wgt1(1:Ndata))
   nrr=sum(wgt1(Ndata+1:Ndata+Nrand))
   Zddd=Zddd/(ndd*ndd*ndd)
   Zrrr=Zrrr/(nrr*nrr*nrr)
   if (loadran) then
       Zddr=1.5_kdkind*Zddr/(3*ndd*ndd*nrr)
       Zdrr=3.0_kdkind*Zdrr/(3*ndd*nrr*nrr)
    else
       Zddr=Zddr/(3*ndd*ndd*nrr)
       Zdrr=Zdrr/(3*ndd*nrr*nrr)
   endif
else
   ndd=float(Ndata)
   nrr=float(Nrand)
   Zddd=Zddd/(ndd*ndd*ndd)
   Zddr=Zddr/(3*ndd*ndd*nrr)
   Zdrr=Zdrr/(3*ndd*nrr*nrr)
   Zrrr=Zrrr/(nrr*nrr*nrr)
endif


end subroutine

SUBROUTINE dist(v1,v2,d)
implicit none
real(kdkind):: v1(3),v2(3), d

d=sqrt((v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2)

return
end SUBROUTINE dist

subroutine deallocate_arrays()
implicit none

if (allocated(my_array)) deallocate(my_array)
if (allocated(my_array2)) deallocate(my_array2)
if (allocated(Zddd)) deallocate(Zddd)
if (allocated(Zddr)) deallocate(Zddr)
if (allocated(Zdrr)) deallocate(Zdrr)
if (allocated(Zrrr)) deallocate(Zrrr)
if (allocated(crap)) deallocate(crap)

if (allocated(p1)) deallocate(p1)
if (allocated(p2)) deallocate(p2)
if (allocated(p3)) deallocate(p3)

if (allocated(wgt1)) deallocate(wgt1)
if (allocated(bins)) deallocate(bins)
if (allocated(tbins)) deallocate(tbins)
if (allocated(v1)) deallocate(v1)
if (allocated(v2)) deallocate(v2)

if (allocated(resultsb)) deallocate(resultsb)
if (allocated(results2)) deallocate(results2)

end subroutine deallocate_arrays

subroutine load_save_randoms()
if(saveran) then
  open(11,file=ranfile,status='unknown',form='unformatted')
    write(11)Zrrr
  close(11)
endif

if(loadran) then
  open(11,file=ranfile,status='unknown',form='unformatted')
    read(11)Zrrr
  close(11)
endif

end subroutine load_save_randoms

subroutine default_params()

  d=3
  wgt=.false.
  logbins=.false.
  nbins=0
  rmin=0.0
  rmax=0.0
  outfile='result.txt'
  RSD=.false.
  nmu=1
  mu1=1
  mu2=1
  
  loadran=.false.
  saveran=.false.

end subroutine default_params

subroutine triplet_bin( )

call dist(v1,v2,vec_dist)

if(RSD) then
        p1=my_array(:,i)
        p2=v1-p1
        p3=v2-p1

       theta2 = ACOS(DOT_PRODUCT(p2,p1)/(SQRT(DOT_PRODUCT(p2,p2))*SQRT(DOT_PRODUCT(p1,p1))))
       theta3 = ACOS(DOT_PRODUCT(p3,p1)/(SQRT(DOT_PRODUCT(p3,p3))*SQRT(DOT_PRODUCT(p1,p1))))

!       if(theta2>pivalby2) theta2=pival-theta2
!       if(theta3>pivalby2) theta3=pival-theta3

       mu1=floor(theta2*odr)+1
       mu2=floor(theta3*odr)+1
       if(mu1>nmu) mu1=nmu
       if(mu2>nmu) mu2=nmu
       if(mu1<=0) mu1=1
       if(mu2<=0) mu2=1
                
endif

end subroutine 

end program  ThreePCF
