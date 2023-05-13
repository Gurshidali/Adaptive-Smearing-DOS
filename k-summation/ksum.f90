module global
        implicit none
        integer :: i,j,Norbs,NRPTS,ir,info,io,iv
        integer :: k,ik, NKPTS,NKTOT,NMESH
        integer :: Nfreq,iw,sm_flag
        integer, allocatable :: irvec(:,:),ndegen(:)

        character*80 :: header

        real*8 :: r1,r2,r3,r4,rdotk,pi
        real*8 :: w0,w,dw,eta,deltak,bfac,A0,A1
        real*8, allocatable:: ikvec(:,:),wtk(:),ndeg(:),sig2(:,:),xw(:)
        real*8,allocatable :: deigv_k(:,:,:)
        real*8,allocatable :: eigv(:),dos(:,:),eigv_k(:,:)

        complex*16 :: z1,fac,ii,zero
        complex*16 :: num,den, determinant
        complex*16,allocatable :: Gf(:,:,:),ctmp(:,:),ctmp2(:,:),cofacm(:,:)
        complex*16,allocatable :: Gf1(:,:,:)
        complex*16, allocatable :: ham_r(:,:,:),ham_k(:,:,:),t_eta_k(:,:,:),eigvec_k(:,:,:),ceigvec(:,:)
        complex*16, allocatable :: dham_k(:,:,:,:), eigvec(:,:), eta_k(:,:,:), eta_k1(:,:,:)
        logical :: write_hk,readkpts
        external determinant
        ! MPI variables                                 
        integer :: ierr, rank, nprocs 
        integer :: jstart, jend, N, nm 
end module global
!---------------------------------------------------------
!---------------------------------------------------------
Program read_ham
        use mpi
        use global

        call MPI_INIT(ierr)                                   
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        ! constants
        pi=4.d0*datan(1.d0)
        zero=dcmplx(0.d0,0.d0)
        ii=dcmplx(0.d0,1.d0)

        if (rank == 0) then
        write(6,*) 'Default broadening factor is 2.0'
        write(6,'(a40)',advance='no') 'Enter broadening factor for dos: '
        read(5,*) bfac
        write(6,*) 'Gaussian(1),Hermite(2) or Lorentzian(3) smearing? '
        write(6,'(a40)',advance='no') 'Enter flag for smearing'
        read(5,*) sm_flag
        write(6,*) 'Suggested NMESH is 25'
        write(6,'(a40)',advance='no') 'Enter NMESH for computing the k-grid'
        read(5,*) NMESH
        end if

        call MPI_BCAST(bfac, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(sm_flag, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(NMESH, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)       

        dw=0.01d0
        Nfreq=2000
        eta=1.D-4
        w0=0.d0 !Fermi energy
!        w0=17.5698d0
        !Energy window is -dw*Nfreq + w0, dw*Nfreq + w0


	! ham_r(NRPTS,Norbs,Norbs) is the real-space tight-binding Hamiltonian that is
        ! being read from H0_tb.dat


        open(unit=15,file='H0_tb.dat',status='unknown')
        read(15,*) header
        read(15,*) header
        read(15,*) header
        read(15,*) header
        read(15,*) Norbs
        read(15,*) NRPTS

        allocate(ham_r(NRPTS,Norbs,Norbs),stat=info)
        allocate(irvec(NRPTS,3),stat=info)
        allocate(ndegen(NRPTS),stat=info)
        allocate(ndeg(NRPTS),stat=info)

        read(15,*) (ndegen(ir), ir=1, NRPTS)
        do ir=1,NRPTS
          read(15,*) irvec(ir,:)      !reads real space vectors one by one...
          do io=1,Norbs**2            !reads matrix elements of each vectors...
            read(15,*) i,j,r1,r2
            z1=dcmplx(r1,r2)
            ham_r(ir,i,j)=z1
          end do
        end do
        close(15)



        ! Normalize ndegen
        r1=0.d0
        do ir=1,NRPTS
          r1=r1+dfloat(ndegen(ir))
        end do
        ndeg=dfloat(ndegen)/r1
          

	!Fourier transform to k-space and get ham_k(NKPTS,Norbs,Norbs)

        readkpts=.false. ! Keep this false if you need agreement 
                         ! with Wannier90 DoS
        if(readkpts) then
          !Read in the K-points from kpts.dat

          open(unit=15,file='kpts.dat',status='unknown')
          read(15,*) NKPTS
          allocate(ikvec(NKPTS,3),stat=info)
          allocate(wtk(NKPTS),stat=info) 

        else

          NKPTS=NMESH**3
          allocate(ikvec(NKPTS,3),stat=info)
          allocate(wtk(NKPTS),stat=info)

          open(unit=8,file='wtk.dat', status='unknown')
          ik=1          
          wtk=1.d0
          do i=1,NMESH
            do j=1,NMESH
              do k=1,NMESH
                ikvec(ik,1)=dfloat(i-1)/dfloat(NMESH)  !3D uniform k-mesh generating...
                ikvec(ik,2)=dfloat(j-1)/dfloat(NMESH) 
                ikvec(ik,3)=dfloat(k-1)/dfloat(NMESH)
                write(8,*) ik, ikvec(ik,:), wtk(ik)
                ik=ik+1
              end do
            end do
          end do
          NKTOT=NKPTS
        end if

        close(8)

        allocate(ham_k(NKPTS,Norbs,Norbs),stat=info)
        allocate(dham_k(NKPTS,3,Norbs,Norbs),stat=info)
        allocate(ctmp(Norbs,Norbs))


        if(readkpts) then
          do ik=1,NKPTS
            read(15,*) i,ikvec(ik,:),wtk(ik)
          end do
          close(15)
        end if


        ! Normalize wtk
        r1=0.d0
        do ik=1,NKPTS
          r1=r1+wtk(ik)
        end do
        NKTOT=int(r1)
        wtk=wtk/r1


        ham_k=zero
        dham_k=zero
        do ik = 1,NKPTS 
          do ir = 1, NRPTS 
            rdotk = 2.d0*pi*dot_product(ikvec(ik,:),dfloat(irvec(ir,:)))
            fac = zexp(ii*rdotk)
            ctmp(:,:)=dfloat(ndegen(ir))*fac*ham_r(ir,:,:)
            ctmp(:,:)=fac*ham_r(ir,:,:)
            ham_k(ik,:,:) = ham_k(ik,:,:) + ctmp(:,:)
            do iv=1,3   
              dham_k(ik,iv,:,:)=dham_k(ik,iv,:,:)+2.d0*pi*ii*dfloat(irvec(ir,iv))*ctmp(:,:) !Derivative Matrix
            end do
          enddo
        enddo
        if (rank == 0)then
         write(6,*) 'COMPUTED H0_K AND Grad(H0_K)'
        endif

        write_hk=.true.

        if(write_hk) then
          open(unit=17,file='ham_k.dat',status='unknown')
          write(17,*) NKPTS
          do ik=1,NKPTS
            write(17,*) ikvec(ik,:)
            do i=1,Norbs
              do j=1,Norbs
                write(17,*) i,j,dreal(ham_k(ik,i,j)),dimag(ham_k(ik,i,j))
              end do
            end do
          end do
          close(17)

        end if ! write_hk = true



        ! USE THIS HAM_K TO GET THE DOS AND COMPARE WITH THE ORIGINAL
        !open(unit=17,file='Ham_k.dat',status='unknown')
        !read(17,*) NKPTS
        !do ik=1,NKPTS
        !  read(17,*) ikvec(ik,:)
        !  do io=1,Norbs**2
        !        read(17,*) i,j,r1,r2
        !        !ham_k(ik,i,j)=dcmplx(r1,r2)
        !  end do
        !end do
        !close(17)
        

          
        allocate(Gf(-Nfreq:Nfreq,Norbs,Norbs))
        allocate(Gf1(-Nfreq:Nfreq,Norbs,Norbs))
        allocate(ctmp2(Norbs-1,Norbs-1))
        allocate(eigv_k(NKPTS,Norbs))
        allocate(eigv(Norbs))
        allocate(deigv_k(NKPTS,3,Norbs))
        allocate(dos(-Nfreq:Nfreq,Norbs))
        allocate(cofacm(Norbs,Norbs))
        allocate(sig2(NKPTS,Norbs))
        allocate(xw(Norbs))
        allocate(t_eta_k(NKPTS,Norbs,Norbs))
        allocate(eta_k1(NKPTS,Norbs,Norbs))
        allocate(eta_k(NKPTS,Norbs,Norbs))
        allocate(eigvec_k(NKPTS,Norbs,Norbs))
        allocate(eigvec(Norbs,Norbs))
        allocate(ceigvec(Norbs,Norbs))

	! Find eigenvalues of H0(k) for each k        
        do ik=1,NKPTS
            ctmp(:,:)=ham_k(ik,:,:)
            call diagonalize(Norbs,ctmp,eigv)
            eigv_k(ik,:)=eigv(:)              !Eigenvalues
            eigvec_k(ik,:,:)=ctmp(:,:)        !Eigenvectors
        end do

	! Find the k-gradient of the eigenvalues of H0(k) for each k        
        do ik=1,NKPTS
          do io=1,Norbs
            ctmp(:,:)=ham_k(ik,:,:)
            do j=1,Norbs
                ctmp(j,j)=ctmp(j,j)-eigv_k(ik,io)  ! |H0_k - Ek| = 0
            end do
            do i=1,Norbs
              do j=1,Norbs
                call cofactor(Norbs,ctmp,ctmp2,i,j)
                cofacm(i,j)=((-1.d0)**(i+j))*determinant(Norbs-1,ctmp2)
              end do
            end do
            do iv=1,3
              num=zero
	      den=zero
              do i=1,Norbs
                do j=1,Norbs
                  num=num+dham_k(ik,iv,i,j)*cofacm(i,j)
                end do
                den=den+cofacm(i,i)
              end do
              deigv_k(ik,iv,io)=dreal(num/den)
            end do ! iv
          end do   ! io
          !write(6,*) 'K-GRADIENT',ik
        end do     ! ik
        
        if (rank == 0)then
         write(6,*) 'COMPUTED DERIVATIVES OF EIGENVALUES'
        endif
   
        ! Adaptive smearing -PHYSICAL REVIEW B 75, 195121 (2007)   
        t_eta_k=zero
        deltak=(2.d0*pi)**3/dfloat(NKTOT)
        do ik=1,NKPTS
          do io=1,Norbs
            sig2(ik,io)=0.d0
            do iv=1,3
              sig2(ik,io)=sig2(ik,io)+(deigv_k(ik,iv,io))**2
            end do
            sig2(ik,io)=bfac*dsqrt(sig2(ik,io))*deltak 
            t_eta_k(ik,io,io)=sig2(ik,io) !t_eta_k, k-dependent broadening parameter is now a diagonal matrix in the orbital basis
          end do
        end do
        if (rank == 0)then
         write(6,*) 'COMPUTED ADAPTIVE SMEARING'
        endif

!!! Task division for MPI
!-------------------------------------------
        N=Nfreq
        jstart=-N+rank*int((2*N)/nprocs)+rank
        jend=jstart+int((2*N)/nprocs)
        if (rank == nprocs-1)then
          jend=N
        endif
!-------------------------------------------

        dos=0.d0
        if(sm_flag.eq.1) then  ! Gaussian smearing PRB 75,195121(2007).    
          do iw=-Nfreq,Nfreq
            w=dw*dfloat(iw)+w0
            do ik=1,NKPTS
              xw(:)=(w-eigv_k(ik,:))/(dsqrt(2.d0)*sig2(ik,:))
              dos(iw,:)=dos(iw,:) + wtk(ik)*dexp(-(xw(:))**2)/sig2(ik,:)
            end do
          end do
          dos=dos/dsqrt(2.d0*pi)
!-------------------------------------------
        else if(sm_flag.eq.2) then   ! first order Hermite smearing
                                     ! PRB 40,3616(1989)
          A0=1.d0/dsqrt(pi)
          A1=-A0/4.d0
          do iw=-Nfreq, Nfreq     ! Frequency sum starts
            w=dw*dfloat(iw)+w0
            do ik=1,NKPTS
              xw(:)=(w-eigv_k(ik,:))/(dsqrt(2.d0)*sig2(ik,:))
              dos(iw,:)=dos(iw,:) + wtk(ik)*dexp(-(xw(:))**2)*&
                      (A0+A1*((xw(:))**2-1.d0))/sig2(ik,:)
            end do
            end do        ! Frequency sum ends
          dos=dos/dsqrt(2.d0*pi)
!--------------------------------------------
        else	! Lorentzian smearing
          do iw=-Nfreq,Nfreq
            w=dw*dfloat(iw)+w0
            do ik=1,NKPTS
              xw(:)=(w-eigv_k(ik,:))/(sig2(ik,:))
              dos(iw,:)=dos(iw,:) + wtk(ik)*1.d0/&
                        (sig2(ik,:)*(xw(:)**2+1.d0))
            end do
          end do
          dos=dos/pi
        end if
!--------------------------------------------
        !Create the eta_k matrix
        eta_k=zero
        do ik=1,NKPTS
           eigvec(:,:)=eigvec_k(ik,:,:) !Unitary matrix U
           ceigvec=dconjg(transpose(eigvec)) ! U_dagger
           ctmp(:,:)=t_eta_k(ik,:,:)
           !Now the product of U_dagger * t_eta_k * U is eta_k(:,:)
           ctmp=MATMUL(ctmp,eigvec)
           ctmp=MATMUL(ceigvec,ctmp)
           eta_k(ik,:,:)=ctmp(:,:)
        end do

       !############ Writing eta_k #################

        open(unit=18,file='eta_k.dat',status='unknown')
        write(18,*) NKPTS
        do ik=1,NKPTS
          write(18,*) ikvec(ik,:)
          do i=1,Norbs
            do j=1,Norbs
              write(18,*) i,j,dreal(eta_k(ik,i,j)),dimag(eta_k(ik,i,j))
            end do
          end do
        end do
        close(18)

!        call cpu_time(t1)

! Non-Interacting Green's function as (w + ii*eta_k - H_0(k))^{-1}
        Gf1=zero
        do iw = jstart, jend  !Frequency sum 
          w=dw*dfloat(iw)+w0
          z1=w+ii*eta
          do ik=1,NKPTS !k-summation
            ctmp(:,:)=ii*eta_k(ik,:,:)-ham_k(ik,:,:)
            do io=1,Norbs
              ctmp(io,io)=z1+ctmp(io,io)
            end do
            call cmatinv(Norbs,ctmp)
            Gf1(iw,:,:)=Gf1(iw,:,:) + wtk(ik)*ctmp(:,:)
          end do  ! k-summ end here.
        end do   ! Frequency sum ends here.
         
        Gf=Gf1
        nm=(2*N+1)*Norbs*Norbs

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(Gf1,Gf, nm, MPI_DOUBLE_COMPLEX,&
                                 MPI_SUM, MPI_COMM_WORLD, ierr)

      
!        if(task_id==0)   write(6,*) '_____________'
!        call cpu_time(t2)
!        if(task_id==0)  write(6,*)'time',t2-t1

        do io=1,Norbs
          r1=0.d0
          do iw=-Nfreq,Nfreq
            r1=r1+dos(iw,io)*dw
          end do
          do iw=-Nfreq,Nfreq
            dos(iw,io)=dos(iw,io)/r1
          end do
          if(rank == 0) write(6,*) 'Spectral Sum = ',io,sum(dos(:,io))*dw
        end do

        open(unit=16,file='dos_Bloch.dat',status='unknown')  ! Total dos in Bloch basis
        open(unit=17,file='dos_Gf.dat',status='unknown')     ! Total dos in Wannier basis
!        open(18, file='dos_test_Gf.dat', status='unknown')

        do iw=-Nfreq,Nfreq
          w=dw*dfloat(iw)+w0
          z1=zero
          r1=0.d0
          do io=1,Norbs
            z1=z1+Gf(iw,io,io)/dfloat(Norbs)
            r1=r1+dos(iw,io)/dfloat(Norbs)
          end do

          write(16,'(2f15.6)') w,r1
          write(17,'(2f15.6)') w,-dimag(z1)/pi
!          write(18,*) w, dreal(z1), -dimag(z1)/pi

        end do
        close(16)
        close(17)
!        close(18)
        
        deallocate(ham_r)
        deallocate(irvec)
        deallocate(ndegen)

        deallocate(ham_k)
        deallocate(dham_k)
        deallocate(ikvec)
        deallocate(wtk)

        deallocate(Gf)
        deallocate(Gf1)
        deallocate(ctmp)
        deallocate(ctmp2)
        deallocate(eigv_k)
        deallocate(deigv_k)
        deallocate(dos)
        deallocate(sig2)
        deallocate(xw)
        deallocate(eigvec_k)
        deallocate(eigvec)
        deallocate(ceigvec)
        deallocate(t_eta_k)
        deallocate(eta_k1)
        deallocate(eta_k)

        call MPI_FINALIZE(ierr)   

end program read_ham


!#####################################################################

subroutine cmatinv(n,a)
  implicit none
  integer :: n,INFO,NRHS,LDA,LDB,i,j
  complex*16 :: a(n,n),x(n,n),b(n,n)
  integer :: IPIV(n)
  complex*16,parameter :: zero=dcmplx(0.d0,0.d0),one=cmplx(1.d0,0.d0)

  NRHS=N
  
  LDA=N

  B=A

  LDB=N

  do i=1,n
    do j=1,n
       A(i,j)=zero
    end do
    A(i,i)=one
  end do
  
  call ZGESV( N, NRHS, B, LDA, IPIV, A, LDB, INFO )

return
end subroutine cmatinv

!****************************************************
subroutine diagonalize(N, A, evals)

character*1 :: JOBZ, UPLO
integer :: N, LDA, LWORK, INFO
real*8 :: evals(N), RWORK(3*N-2)
complex*16 :: A(N,N), WORK(2*N-1)

JOBZ='V'
UPLO='U'
LDA=N
LWORK=2*N-1


call ZHEEV(JOBZ, UPLO, N, A, LDA, evals, WORK, LWORK, RWORK, INFO)

!write(6,*) 'INFO=',INFO

return
end subroutine diagonalize
!****************************************************
complex*16 function determinant(N, A)

integer :: N
integer :: LDA, IPIV(N), INFO, i
complex*16 :: A(N,N),z1,one

one=dcmplx(1.d0,0.d0)

LDA=N
call zgetrf(N,N,A,LDA, IPIV, INFO)

z1=one
do i=1,N
  if(IPIV(i).ne.i) then
    z1=-z1*A(i,i)
  else
    z1=z1*A(i,i)
  end if
end do

determinant = z1

return
end function determinant


!****************************************************

subroutine cofactor(N, A, C_A, MI, MJ)

integer :: N,MI,MJ
complex*16 :: A(N,N), C_A(N-1,N-1)
integer :: i,j,k,l

    k = 1
    do i=1, n
      if (i .ne. mI) then
         l = 1
         do j=1, n
           if (j .ne. mJ) then
             C_A(k,l) = A(i,j)
             l = l+ 1
           end if
         end do
         k = k+ 1
      end if
    end do
return
end subroutine cofactor
