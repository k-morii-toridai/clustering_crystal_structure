module neib
  use m_ftox
  !     ! for given plat,pdiff rmax, return dislist
  integer,protected:: ng
  real(8),allocatable,protected:: dislist(:)
contains
  subroutine neibourpair(ix,iy,alat,plat,rbas21,rmax,atomid,ipr)
    !! For given plat and rbas21=ratom2-ratom1 (cartesian),
    !! we give distance lists, within rmax.
    !! -----------------------------------------------------------------------------------      
    !! Show possible distance of lattice as  |rbas21- (n1*plat1 + n2*plat2 + n3*plat3)|< rmax
    implicit none
    real(8):: plat(3,3),bohr,xv,alat,pi=4*atan(1d0),qlat(3,3),rmax,pwgmax,qdiff(3) &
         ,dum,platx(3,3),pwgmin,gmax0,rvec(3),distance,rbas21(3),rmax2,dis2
    integer,allocatable:: igv2x(:,:),kv_iv(:)
    integer:: ig,napw,nn1,nn2,nn3,mshlst,nmin(3),nmax(3),n1,n2,n3,ndatx,ipr,ix,iy
    character*24:: atomid(*)
    !c      print *, plat(1:3,1),plat(1:3,2),plat(1:3,3),rbas21(1:3),rmax
    !c      PLAT(1:3,1)=[0.5, 0.5, 1.0]
    !c      PLAT(1:3,2)=[0.5, 1.0, 0.5]
    !c      PLAT(1:3,3)=[1.0, 0.5, 0.5]
    !c      rbas1=[0,0,0]       !cartesian atom1
    !c      rbas2=[.5,.5,0.0]    !cartesian atom2
    !
    rmax2=rmax**2
    call gvlstnx(plat(:,1),plat(:,2),plat(:,3),rbas21,rmax, nmin(1),nmax(1))
    call gvlstnx(plat(:,2),plat(:,3),plat(:,1),rbas21,rmax, nmin(2),nmax(2))
    call gvlstnx(plat(:,3),plat(:,1),plat(:,2),rbas21,rmax, nmin(3),nmax(3))
    !! For given rbas21(1:3) (cartesian),
    !c$$$      print *,'rbas21 =',rbas21
    !c$$$      print *,'plat1 =',plat(:,1)
    !c$$$      print *,'plat2 =',plat(:,2)
    !c$$$      print *,'plat3 =',plat(:,3)
!    print *,'nmin1 mnax1 =',nmin(1),nmax(1),'nmin2 mnax2 =',nmin(2),nmax(2),'nmin3 mnax3 =',nmin(3),nmax(3)
    ndatx=(nmax(1)-nmin(1)+1)*(nmax(2)-nmin(2)+1)*(nmax(3)-nmin(3)+1)
    if(allocated(dislist)) deallocate(dislist)
    allocate(dislist(ndatx))
    ig=0
    do n1 = nmin(1),nmax(1)
       do n2 = nmin(2),nmax(2)
          do n3 = nmin(3),nmax(3)
             rvec = rbas21 + (plat(:,1)*n1 + plat(:,2)*n2 + plat(:,3)*n3)
             dis2=sum(rvec*rvec)
             if( dis2<= rmax2) then
                write(ipr,"(' ',2i6,f12.6,' ',3f11.6,' ',3i4,'  ',a,' ',a)")&
                     ix,iy,alat*sqrt(dis2),alat*rvec,n1,n2,n3,trim(atomid(ix)),trim(atomid(iy))
                ig=ig+1
                dislist(ig)=sqrt(dis2)
             endif
          enddo
       enddo
    enddo
    ng=ig
    !stop 'vvvvvvvvvvvvvvvvvvvvv'
  end subroutine neibourpair
end module neib
!=============================================================================================
subroutine gvlstnx(q0,q1,q2,qp,rmax, nmin,nmax)  !! Multiples of primitive vectors.
  !! we look for allow n; n is nmin<n< nmax, satisfying  |qp+n1*q1+n2*q2+n3|< rmax
  !C     ----------------------------------------------------------------------
  !Ci   q0    : first  primitive vector
  !Ci   q1,q2 : other two primitive vectors.
  !Ci   qp    :  added to G vectors: sphere is centered G=qp.
  !Ci  rmax  : cutoff Radius
  !Co Outputs
  !Co   nmin  :search for lattice vectors limited to (nmin..nmax)*q0
  !Co   nmax  :search for lattice vectors limited to (nmin..nmax)*q0
  !! NOTE: n1=constant makes a plane spanned by q1 and q2. We calculate distance to the plane 
  !! (n1,n2,n3) on the plane has a distane larger than the distance.
  implicit none
  integer nmin,nmax,nt1,nt2
  real(8):: qperp(3),ddot,qqperp,q0(3),q1(3),q2(3),qp(3),rmax,distance1,distanceqp
  !C ... qperp = q1 x q2 / |q1 x q2| ; qqperp = q . qperp
  qperp(1)  = q1(2)*q2(3) - q1(3)*q2(2) !external product qprep=q1 x q2
  qperp(2)  = q1(3)*q2(1) - q1(1)*q2(3) !this is normal to the plane for given n1
  qperp(3)  = q1(1)*q2(2) - q1(2)*q2(1)
  qperp = qperp/sqrt(sum(qperp*qperp))  !normalize vector. sum(a,b) is innner product
  distanceqp = sum(qp*qperp) !signed distance to qp from origin
  distance1=   sum(q0*qperp) !signed distance to plane for n1=1 from origin
  !! we search n1 for |distance1*n1 - distanceqp|<rmax
  !! fortran int gives an integer closer to zero for given float number.
  nmax= floor  ( max((rmax-distanceqp)/distance1,(-rmax-distanceqp)/distance1) )
  nmin= ceiling( min((rmax-distanceqp)/distance1,(-rmax-distanceqp)/distance1) )
end subroutine gvlstnx

!!----------------------------------------------------------      
program neib1
  !     ! cell
  use neib,only: ng, dislist, Neibourpair
  use m_ftox
  implicit none
  integer::u,v,i,npair,ix,id1,id2,irr,iu, natom,system,lkeyw,iarg,iy,ntype,n,nnn
  integer,allocatable::nsite(:)
  real(8):: plat(3,3),rbas21(3),rmax
  character*120:: aaa,bbb,filename,tmpline,fout,strn,cifid,keyw,atypeall,postype
  character*24:: atomid(10000),filetype,atype(200)
  real(8):: pos(3,10000),alat
!
  if( iargc()/=2) then
     write(*,*)'usage: Supply POSCARname and rmax. '
     write(*,*)'     : For example, >neib1 poscar.1 5.0'
     stop 
  endif
  iarg=1
  call getarg(iarg,strn)
  read(strn,*) cifid
  iarg=2
  call getarg(iarg,strn)
  read(strn,*) rmax
  
  filetype=cifid
  postype='direct'
  if(.false.) then !trim(filetype)/='POSCAR') then
     fout=trim(cifid)//'.cell'
     open(newunit=iu, file=trim(fout))
     irr=system('cif2cell '//trim(cifid)//'.cif'//' > '//trim(fout))
     do 
        read(iu,"(a)") tmpline
        !print *,'mmmmm'//trim(tmpline)//'vvv'
        if(trim(tmpline)=='Bravais lattice vectors :') then
           read(iu,*) plat(1:3,1)
           read(iu,*) plat(1:3,2)
           read(iu,*) plat(1:3,3)
        elseif(trim(tmpline)=='All sites, (lattice coordinates):') then
           read(iu,*)
           ix=0
           do
              ix=ix+1
              read(iu,"(a)",end=1001,err=1001) tmpline
              if(trim(tmpline)=='') goto 1001
              read(tmpline,*) atomid(ix), pos(:,ix)
              !print *,'ix=',ix,trim(atomid(ix)), pos(:,ix)
              if(trim(atomid(ix))=='') then
                 ix=ix-1
                 goto 1001
              endif
           enddo
        endif
     enddo
     ix=ix-1
     natom=ix
  else !if(trim(filetype)=='poscar') then !only for Si case
     !write(*,ftox) 'Read POSCAR!'
     fout=trim(cifid) !'POSCAR'
     open(newunit=iu, file=trim(fout))
     read(iu,*)
     read(iu,*) alat
     read(iu,*) plat(1:3,1)
     read(iu,*) plat(1:3,2)
     read(iu,*) plat(1:3,3)
     read(iu,"(a)") atypeall
     atypeall=adjustl(atypeall)
     !write(*,*)'atypeall=',trim(atypeall)
     i=0
     do
        i=i+1
        n=scan(atypeall,' ')
        !write(*,*)'atypeall=',trim(atypeall(1:n)),'###',n
        read(atypeall(1:n),"(a)") atype(i)
        atypeall = adjustl(atypeall(n+1:))
        !write(*,*)'out=',i,trim(atype(i))
        if(len_trim(atypeall)==0) exit
     enddo
     ntype=i
     allocate(nsite(ntype))
     read(iu,*) nsite
     natom = sum(nsite)
     write(*,ftox)nsite,' !nsite'
     do i=1,ntype
        write(*,*)'readin type',i,trim(atype(i)),nsite(i)
     enddo
     nnn=0
     do i=1,ntype
        atomid(nnn+1:nnn+nsite(i))=atype(i)
        nnn=nnn+nsite(i)
     enddo
     read(iu,*) postype
     call lower_case(postype,postype)
     write(*,ftox)trim(postype),' !postype'
     do ix=1,natom
        read(iu,*) pos(:,ix)
     enddo
     write(*,ftox) alat, ' ! unit (angstrom)'
  endif    
1001 continue
  do ix=1,natom
     write(*,ftox) ix,trim(atomid(ix)),ftof(pos(:,ix))
  enddo
!  
  !write(*,ftox)'rmax plat natom:'
  write(*,ftox) ftof(rmax),' !rmax (angstrom)'
  write(*,ftox) ftof(plat(1:3,1)), ' !plat(:,1)'
  write(*,ftox) ftof(plat(1:3,2)), ' !plat(:,2)'
  write(*,ftox) ftof(plat(1:3,3)), ' !plat(:,3)'
  write(*,ftox) natom, ' !natom'
  if(trim(postype)=='cartesian') then
     write(*,ftox) 'Atom position: Cartesian coordinate'
     do ix=1,natom
        write(*,ftox) ix,trim(atomid(ix)),ftof(pos(:,ix))
     enddo
  else    
!     do ix=1,natom
!        write(*,ftox) ftof(pos(:,ix)),ix,trim(atomid(ix)),'  !Atom position: fractional(direct)' ! / cartetian'
!     enddo
     do ix=1,natom
        write(*,ftox) ftof(matmul(plat,pos(:,ix))),'pos atomid=',ix,trim(atomid(ix)),'  !Atom position: cartetian'
     enddo
  endif
!  open(newunit=u, file=trim(filename)//'_inp.dat', action='read')
!   read(u,*) plat(1:3,1)
!  read(u,*) plat(1:3,2)
!  read(u,*) plat(1:3,3)
!  read(u,*) rmax  

  !write(*,ftox) natom*natom, ' !natom*natom'
  !do ix=1,natom
  !   do iy=1,natom
  !      if(trim(postype)/='cartesian') rbas21= matmul(plat,pos(:,iy)- pos(:,ix)) !fractional coodinate
  !      if(trim(postype)=='cartesian') rbas21= pos(:,iy)- pos(:,ix) !cartesian coodinate
  !      write(*,ftox)ix,iy,ftof(rbas21,6),' ! ix iy r2-r1 (cartesian)'
  !   enddo
  !enddo
  
  open(newunit=v, file=trim(cifid)//'.nnlist')
  do ix=1,natom
     do iy=1,natom
        if(trim(postype)/='cartesian') rbas21= matmul(plat,pos(:,iy)- pos(:,ix)) !fractional coodinate
        if(trim(postype)=='cartesian') rbas21= pos(:,iy)- pos(:,ix) !cartesian coodinate
        !write(v,ftox)ix,iy,ftof(rbas21,6),' ',ftof(rmax),'! ix iy r2-r1 rmax'
        call neibourpair(ix,iy,alat,plat,rbas21,rmax/alat,atomid,v)
     enddo
  enddo
  close(u)
  close(v)
  write(*,ftox)'OK! Wrote ',trim(cifid)//'.nnlist'
end program neib1

subroutine lower_case(uword,lword)
character (len=*) uword,lword
integer i,ic,nlen
nlen = len_trim(uword)
do i=1,nlen
ic = ichar(uword(i:i))
if (ic >= 65 .and. ic <= 90) lword(i:i) = char(ic+32)
end do
end subroutine
