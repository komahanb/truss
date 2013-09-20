! Copyright (C) 2002, 2010 Carnegie Mellon University and others.
! All Rights Reserved.
! This code is published under the Eclipse Public License.
!
!    $Id: hs071_f.f.in 1861 2010-12-21 21:34:47Z andreasw $
!
! =============================================================================
!
!
!  This file contains routines to define a small Rosenbrock test problem.
!
!  f=0 at the optimal solution for x_i=1 for all i!
!
! =============================================================================
!
!
! =============================================================================
!
!                            Main driver program
!
! =============================================================================
!
      program example
!
      implicit none
!
!     include the Ipopt return codes
!
      include 'IpReturnCodes.inc'
!
!     Size of the problem (number of variables and equality constraints)
!
      integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  (N = 3, M = 3, NELE_JAC = 6, NELE_HESS = 6)
      parameter  (IDX_STY = 1 )

!
!
!     Space for multipliers and constraints
!
      double precision LAM(M)
      double precision G(M)
!
!     Vector of variables
!
      double precision X(N)
!
!     Vector of lower and upper bounds
!
      double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
      double precision G_L(M), G_U(M)
!
!     Private data for evaluation routines
!     This could be used to pass double precision and integer arrays untouched
!     to the evaluation subroutines EVAL_*
!
      double precision DAT(2)
      integer IDAT(1)
!
!     Place for storing the Ipopt Problem Handle
!
      integer*8 IPROBLEM
      integer*8 IPCREATE
!
      integer IERR
      integer IPSOLVE, IPADDSTROPTION
      integer IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
!
      double precision F
      integer i

      double precision  infbound
      parameter        (infbound = 1.d+20)

!
!     The following are the Fortran routines for computing the model
!     functions and their derivatives - their code can be found further
!     down in this file.
!
      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS
!!
!!     The next is an optional callback method.  It is called once per
!!     iteration.
!!
      external ITER_CB
!
!     Set initial point and bounds:
!

      do i=1,N
          X(i) = 2.0
          X_L(i) = 0.0
          X_U(i) = infbound
      end do
!!$      data X   / -1.2d0, 1.d0 /
!!$      data X_L / -1d1, -1d1 /
!!$      data X_U / 1d1, 1d1 /
!
!     Set bounds for the constraints
!

      do i=1,M
         G_L(i)=-infbound
         G_U(i)=0.d0
      end do

!!$      data G_L / 25d0, 40d0 /
!!$      data G_U / 1d40, 40d0 /

!
!     First create a handle for the Ipopt problem (and read the options
!     file)
!
      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif
!
!     Open an output file
!
      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 50)
      if (IERR.ne.0 ) then
         write(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif
!

!!
!!     Set a callback function to give you control once per iteration.
!!     You can use it if you want to generate some output, or to stop
!!     the optimization early.
!!
      call IPSETCALLBACK(IPROBLEM, ITER_CB)

!
!     As a simple example, we pass the constants in the constraints to
!     the EVAL_C routine via the "private" DAT array.
!
      DAT(1) = 1.5 !FS
!      DAT(2) = 0.d0
!
!     Call optimization routine
!

      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
!
!     Output:
!
      if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
         write(*,*)
         write(*,*) 'The solution was found.'
         write(*,*)
         write(*,*) 'The final value of the objective function is ',F
         write(*,*)
         write(*,*) 'The optimal values of X are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'X  (',i,') = ',X(i)
         enddo
!!$         write(*,*)
!!$         write(*,*) 'The multipliers for the lower bounds are:'
!!$         write(*,*)
!!$         do i = 1, N
!!$            write(*,*) 'Z_L(',i,') = ',Z_L(i)
!!$         enddo
!!$         write(*,*)
!!$         write(*,*) 'The multipliers for the upper bounds are:'
!!$         write(*,*)
!!$         do i = 1, N
!!$            write(*,*) 'Z_U(',i,') = ',Z_U(i)
!!$         enddo
         write(*,*)
         write(*,*) 'The multipliers for the equality constraints are:'
         write(*,*)
         do i = 1, M
            write(*,*) 'LAM(',i,') = ',LAM(i)
         enddo
         write(*,*)
      else
         write(*,*)
         write(*,*) 'An error occoured.'
         write(*,*) 'The error code is ',IERR
         write(*,*)
      endif
!

 9000 continue

!
!     Clean up
!
      call IPFREE(IPROBLEM)
      stop
!
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000
    end program example
!
! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X,I
      double precision F, X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR


      double precision  rho, L, sigmay, pi, p, E, Fs

!!$      write(*,*)
!!$      print *,'F',X,NEW_X
!!$      write(*,*)

      rho=0.2836
      sigmay=36260.0
      p=25000.0
      L=5.0
      E=30e6
      pi=4.0*atan(1.0)

      fs=DAT(1)

!!$
!!$
!!$      F=0.d0
!!$      do i=1,N-1
!!$         F=F+100.d0*(X(i+1)-X(i)**2)**2+(1-X(i))**2
!!$      end do
!!$

      f = rho*x(1)*L+rho*x(2)*sqrt(L**2+x(3)**2)

      IERR = 0

      return
      end
!
! =============================================================================
!
!                Computation of gradient of objective function
!
! =============================================================================
!
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X,i
      double precision GRAD(N), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      double precision  rho, L, sigmay, pi, p, E, Fs

      rho=0.2836
      sigmay=36260.0
      p=25000.0
      L=5.0
      E=30e6
      pi=4.0*atan(1.0)
      fs=dat(1)

!!$
!!$      GRAD(1)=-200.d0*(x(2)-x(1)**2)*2.d0*x(1)-2.d0*(1-x(1))
!!$      do i=2,N-1
!!$         GRAD(i)=200.d0*(x(i)-x(i-1)**2)-200.d0*(x(i+1)-x(i)**2)*2.d0*x(i)-2.d0*(1-x(i))
!!$      end do
!!$      GRAD(n)=200.d0*(x(n)-x(n-1)**2)
!!$

!      print*,x
      grad(1) = rho*L
      grad(2) = rho*sqrt(L**2+x(3)**2)
      grad(3) = rho*x(2)*x(3) / sqrt(L**2+x(3)**2)
 !     print*,GRAD
      IERR = 0
      return
      end
!
! =============================================================================
!
!                     Computation of equality constraints
!
! =============================================================================
!
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X, M
      double precision G(M), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      double precision  rho, L, sigmay, pi, p, E, Fs 

      rho=0.2836
      sigmay=36260.0
      p=25000.0
      L=5.0
      E=30e6
      pi=4.0*atan(1.0)
      FS=dat(1)

!      G(1) = X(1)*X(2)*X(3)*X(4) - DAT(1)
!      G(2) = X(1)**2 + X(2)**2 + X(3)**2 + X(4)**2 - DAT(2)

!      print*,x
      G(1) = p*Fs*sqrt(L**2+x(3)**2) / (x(2)*x(3)*sigmay) - 1.0
      G(2) = p*Fs*L / (x(1)*x(3)*sigmay) - 1.0
      G(3) = 4.0*p*Fs*L**3 / (x(1)**2*x(3)*E*pi) - 1.0
!      print*,G
      IERR = 0
      return
      end
!
! =============================================================================
!
!                Computation of Jacobian of equality constraints
!
! =============================================================================
!
      subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,IDAT, DAT, IERR)
      integer TASK, N, NEW_X, M, NZ
      double precision X(N), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*),dc(M,N)
      integer IDAT(*)
      integer IERR
!
!     structure of Jacobian:
!
      integer AVAR1(8), ACON1(8)

      double precision  rho, L, sigmay, pi, p, E, Fs 
      
      FS=dat(1)
      rho=0.2836
      sigmay=36260.0
      p=25000.0
      L=5.0
      E=30e6
      pi=4.0*atan(1.0)

      if( TASK.eq.0 ) then 
         !
         !     structure of Jacobian:
         !
         ACON(1) = 1
         AVAR(1) = 2


         ACON(2) = 1
         AVAR(2) = 3

         ACON(3) = 2
         AVAR(3) = 1

         ACON(4) = 2
         AVAR(4) = 3

         ACON(5) = 3
         AVAR(5) = 1

         ACON(6) = 3
         AVAR(6) = 3
         
      else

         !---- GRADIENT OF CONSTRAINTS

         dc(:,:)=0.0

         dc(1,1) =  0.0
         dc(1,2) = -p*Fs*sqrt(L**2+x(3)**2) / (x(2)**2*x(3)*sigmay)
         dc(1,3) = -p*Fs*L**2/sqrt(L**2/x(3)**2+1.0)/ (x(2)*x(3)**3*sigmay)

         dc(2,1) = -p*Fs*L / (x(1)**2*x(3)*sigmay)
         dc(2,2) =  0.0
         dc(2,3) = -p*Fs*L / (x(1)*x(3)**2*sigmay)

         dc(3,1) = -8.0*p*Fs*L**3 / (pi*E*x(1)**3*x(3))
         dc(3,2) =  0.0
         dc(3,3) = -4.0*p*Fs*L**3 / (pi*E*x(1)**2*x(3)**2)

         A(1)=dc(1,2)
         A(2)=dc(1,3)

         A(3)=dc(2,1)
         A(4)=dc(2,3)

         A(5)=dc(3,1)
         A(6)=dc(3,3)    

      end if

      IERR = 0
      return
    end subroutine EV_JAC_G

!
! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
      subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NEW_LAM, NNZH, i, ir,j
      double precision X(N), OBJFACT, LAM(M), HESS(NNZH),OBJHESS(NNZH),CONHESS(M,NNZH)

      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      double precision  rho, L, sigmay, pi, p, E, Fs 
      double precision :: hesstmp

      rho=0.2836
      sigmay=36260.0
      p=25000.0
      L=5.0
      E=30e6
      pi=4.0*atan(1.0)

      FS=dat(1)

      if( TASK.eq.0 ) then
!
!     structure of sparse Hessian (lower triangle):
!

         IRNH(1) = 1
         ICNH(1) = 1

         IRNH(2) = 2
         ICNH(2) = 2

         IRNH(3) = 3
         ICNH(3) = 3

         IRNH(4) = 2
         ICNH(4) = 1

         IRNH(5) = 3
         ICNH(5) = 2
         
         IRNH(6) = 3
         ICNH(6) = 1
        
      else
!
!     calculate Hessian:
!
         !Objective function
        objhess(1)=0.0
        objhess(2)=0.0
        objhess(3)=(rho*x(2)/sqrt((L**2) + (x(3)**2))) - (rho*x(2)*(x(3)**2)/(((L**2) + (x(3)**2))**(3.0/2.0)))
        objhess(4)=0.0
        objhess(5)=(rho*x(3))/sqrt((L**2) + (x(3)**2))
        objhess(6)=0.0


        ! first constraint

        conhess(1,1)=0.0
        conhess(1,2)=(2.0*Fs*p*((L**2)+(x(3)**2))**(1.0/2.0))/(sigmay*x(2)**3*x(3))
        conhess(1,3)=(2.0*Fs*p*((L**2) + (x(3)**2))**(1.0/2.0))/(sigmay*x(2)*x(3)**3) - (Fs*p)/(sigmay*x(2)*x(3)*((L**2) + (x(3)**2))**(1.0/2.0)) - (Fs*p*x(3))/(sigmay*x(2)*((L**2) + (x(3)**2))**(3.0/2.0))
        conhess(1,4)=0.0
        conhess(1,5)=(Fs*p*((L**2) + (x(3)**2))**(1.0/2.0))/(sigmay*(x(2)**2)*(x(3)**2)) - (Fs*p)/(sigmay*(x(2)**2)*(L**2 + (x(3)**2))**(1.0/2.0))
        conhess(1,6)=0.0
        
        ! Second constraint

        conhess(2,1)=(2.0*Fs*L*p)/(sigmay*(x(1)**3)*x(3))
        conhess(2,2)=0.0
        conhess(2,3)=(2.0*Fs*L*p)/(sigmay*x(1)*(x(3)**3))
        conhess(2,4)=0.0
        conhess(2,5)=0.0
        conhess(2,6)=(Fs*L*p)/(sigmay*(x(1)**2)*(x(3)**2))

        ! Third constraint

        conhess(3,1)=(24.0*Fs*(L**3)*p)/(E*pi*(x(1)**4)*x(3))
        conhess(3,2)=0.0
        conhess(3,3)=(8.0*Fs*(L**3)*p)/(E*pi*(x(1)**2)*(x(3)**3))
        conhess(3,4)=0.0
        conhess(3,5)=0.0
        conhess(3,6)=(8.0*Fs*(L**3)*p)/(E*pi*(x(1)**3)*(x(3)**2))
        
        ! Assemble

        HESS(:)=0.0
        do i=1,NNZH
           hesstmp=0.0
           do j=1,m
              hesstmp=hesstmp+lam(j)*conhess(j,i)
           end do
           hess(i)=hesstmp+objhess(i)
        end do

      IERR = 0

      endif
      return
    end subroutine EV_HESS
!
! =============================================================================
!
!                   Callback method called once per iteration
!
! =============================================================================
!
      subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,DAT, ISTOP)
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*)
      integer IDAT(*)
      integer ISTOP

      if (ITER_COUNT .eq.0) then
         write(*,*) 
         write(*,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
      end if

      write(*,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU
      if (ITER_COUNT .gt. 1 .and. DNORM.le.1D-8) ISTOP = 1

      return
    end subroutine ITER_CB
