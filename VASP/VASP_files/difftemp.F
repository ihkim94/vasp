!    ferroelectriparaelectric phase boundary
!    temperature change incorporated

      program FE_PE
#include "petsc/finclude/petsc.h"
      USE petsc
      implicit none

!-------------------- PETSc -------------------------------------------
!#include "petsc/finclude/petsc.h"
!----------------------------------------------------------------------

      integer maxnd,maxele,maxfc,maxsol,maxdbc
      integer maxcon,ndim,nnds,nints,ndpn,maxibcs,minibcs
      parameter (maxnd=167000,maxele=160000,maxfc=3,maxcon=2,
     1           maxdbc=170000,maxsol=100,
     2           ndim=3,nnds=8,nints=8,ndpn=8)
*     ndim = number of spatial dimensions
*     nnds = number of nodes per element
*     nints = number of integration stations per element
*     ndpn = number of degrees of freedom per node
      real*8 node(maxnd,ndim),dbcval(maxdbc),fbcval(maxfc)
      real*8 u(maxnd,ndim),phi(maxnd),P(3,maxnd,ndim),theta(maxnd)
      real*8 NN(nnds,nints),B(maxele,nnds,nints,ndim),det(maxele,nints)
      real*8 strn(maxele,nints,ndim,ndim),E(maxele,nints,ndim)
      real*8 Pint(maxele,nints,ndim),PintPrev(maxele,nints,ndim)
      real*8 gradP(maxele,nints,ndim,ndim)
      real*8 thint(maxele,nints),gradTheta(maxele,nints,ndim)
      real*8 strs(maxele,nints,ndim+1,ndim+1),svm(maxele,nints)
      real*8 D(maxele,nints,ndim)
      real*8 psiP(maxele,nints,ndim),ksi(maxele,nints,ndim,ndim)
      real*8 s(maxele,nints),qt(maxele,nints,ndim)
      real*8 nE(maxnd,ndim),nD(maxnd,ndim),nS(maxnd)
      real*8 nstrn(maxnd,ndim,ndim),nstrs(maxnd,ndim,ndim),nsvm(maxnd)
      real*8 ngradT(maxnd,ndim),nflux(maxnd,ndim),nmat(maxnd)
      real*8 conA(maxdbc,maxcon),conB(maxdbc,maxcon)
      real*8 unitconB(maxdbc,maxcon)
!      real*8 A(maxwk)
      real*8 drhs(maxnd*ndpn),f(maxnd*ndpn),finte(maxnd*ndpn)
      real*8 unitdbc(maxdbc),unitfbc(maxfc),sumsqu,sumsqf,ssu1
      real*8 t,dt,beta_P,T0
      real*8 drhs0(maxnd*ndpn),a0
      real*8 uapp,phiapp,Qapp,avgstrain,Tapp,Fapp
      real*8 avg_stress(ndim,ndim),avg_D(ndim),avg_E(ndim),avg_S,avg_T
      real*8 Px,Py,Pz,P0,theta_init,Psi(maxele,nints),avgTtot
      real*8 exx,eyy,exy,es,ezz,PsiTot,h,w,depth,area!,elthick
      integer nodes,neles,ele(maxele,8),ndofs,mpc(maxdbc,maxcon,2)
      integer mpcnd(maxdbc,2),ncon,gcon(maxnd,ndpn),mpcnum(maxdbc)
      integer fbcnd(maxfc),fbcdof(maxfc),ndbcs,nfbcs
!      integer maxa(maxnd*ndpn+1),mht(maxnd*ndpn),nwk,
      integer EC0cop1(maxele)!,fillet
      integer isol,control
      integer dbcnd(maxdbc),dbcdof(maxdbc),ibcs,ic
      integer nx,ny,gdof1,flagT,flagNR,controlstart
      integer i,j,k,ii,jj,opt,PsiPrint,SPrint,resultfilecount!,n,m
      !,rhsfile
      integer modwrite
      character(7) fnum
      character(4) ou

      integer cpu(nnds,ndpn),dof,loc(maxnd*ndpn),its
      logical eleproc(maxele)

      !----- PETSC variable---------
      PC pc
      Mat PetA
      Vec Petx,Petb
      Vec,POINTER :: xSEQ
!      Vec DD
      KSP ksp
      PetscScalar zero
      PetscScalar Petdrhs(maxnd*ndpn)
      PetscInt Istart,Iend
      PetscErrorCode ierr
      PetscMPIInt rank,sizeof
      PetscScalar x_array(1)
      PetscOffset i_x
      VecScatter sca
!      PetscViewer rhsview
      !-----------------------------

      zero=0.d0

      ! output written to file results & paraview with following suffix
      write(ou,910)'-AR1'
 910  format(A)

      !----- intialize PETSC -------------------------------
      CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      CHKERRA(ierr)
      CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      !rank is the processor number for each processor, zero based
      CALL MPI_Comm_size(PETSC_COMM_WORLD,sizeof,ierr)
      !sizeof=#of processors
      print*, "rank is ", rank, "size is", sizeof
      !------------------------------------------------------


      ! CURIE TEMPERATURE
      T0=273.00000000d0+124.99725d0

      a0=1.00d0
      P0=0.181144132267376282d0

      Px=1.d0
      Py=1.d0
      Pz=1.d0

      es=0.00d0
      exy=0.d0

      exx=-0.d0
      eyy=-0.d0
      ezz=0.0d0

      theta_init=2.32d0
      Qapp=0.d0

      dt=1.d0

      ! read input from nodes, elements, fbcs, dbcs
      call readinp(ndim,nnds,ndpn,maxnd,maxele,maxfc,maxdbc,
     1             maxcon,nodes,neles,ndbcs,nfbcs,ndofs,gcon,ele,
     2             fbcnd,fbcdof,dbcnd,dbcdof,node,fbcval,dbcval,
     3             mpc,conA,conB,mpcnum,ncon,mpcnd,nx,ny,rank)
      if (rank .EQ. 0) then
        print*,'********************************'
        print*,'****** INPUT READ SUCCESS ******'
        print*,'********************************'
        print*
      endif

 922  format(A10,I6,A10,I6,A10,I7)
 923  format(A10,F6.2,A10,F6.2,A10,F6.2)
      ! width, height, area of model, given that it is rect or parallelogram
      h=0.d0
      depth=0.d0
      w=0.d0
      !we assume min is 0,0
      do i=1,nodes
        if (node(i,3).gt.h) h=node(i,3)
        if (node(i,2).gt.depth) depth=node(i,2)
        if (node(i,1).gt.w) w=node(i,1)
      enddo
      !elthick=1.d0
      !fillet=1
      !w=w-2.d0*elthick
      !(w gives the width of the ec material, doesn't include electrode)

      if (rank .EQ. 0) then
        print*,'************************************************'
        print 923,'w=',w,'d=',depth,'h=',h
        print*
        print 922,'nodes',nodes,'neles',neles
        print*
        print 922,'ncon',ncon,'ndofs',ndofs
        print*,'************************************************'
        print*
      endif

      do i=1,neles
        EC0cop1(i)=0
      enddo

!      if (fillet.eq.0) then
!        call choosemat(ndim,nnds,maxnd,maxele,neles,ele,
!     1             node,elthick,h,w,EC0cop1)
!      else
!        call choosematFillet(ndim,nnds,maxnd,maxele,neles,ele,
!     1             node,elthick,h,w,EC0cop1)
!        if (rank.eq.0) then
!          print*,"using fillet choosemat regime"
!          print*
!        endif
!      endif

!      n=0
!      m=0
!      do i=1,neles
!        if (EC0cop1(i).eq.1) n=n+1
!        if (EC0cop1(i).eq.0) m=m+1
!      enddo

!      print*,"number of elements in electrode=",n
!      print*,"number of elements in EC material=",m
!      print*,m+n,"=",neles



      !-----------------------------------------------
      !------ Declare PETsc components ---------------
      !-----------------------------------------------
      if (rank.eq.0) print*,'before matmpiaij'
      CALL MatCreate(PETSC_COMM_WORLD,PetA,ierr)
      CALL MatSetSizes(PetA,PETSC_DECIDE,PETSC_DECIDE,
     1                 ndofs,ndofs,ierr)
      if (sizeof.eq.1) then
        CALL MatSetType(PetA,MATSEQAIJ,ierr)
        CALL MatSeqAIJSetPreallocation(PetA,
!     1        int(ndofs/1),PETSC_NULL_INTEGER,ierr)
     1       int(ndofs/10),PETSC_NULL_INTEGER,ierr)
        !2nd entry is number of nonzeros per row
      else
          CALL MatSetType(PetA,MATMPIAIJ,ierr)
!          CALL MatMPIAIJSetPreallocation(PetA,
!     1        int(ndofs/(sizeof*2)),PETSC_NULL_INTEGER,
!     2       int(ndofs*(sizeof-1)/(sizeof*100)),PETSC_NULL_INTEGER,ierr)
          CALL MatMPIAIJSetPreallocation(PetA,
     1        int(ndofs/(sizeof)),PETSC_NULL_INTEGER,
     2      int(ndofs*(sizeof-1)/(5*sizeof)),PETSC_NULL_INTEGER,ierr)
      endif
      if (rank.eq.0) print*,'after matmpiaij'
      CALL MatSetFromOptions(PetA,ierr)
      CALL MatSetUp(PetA,ierr)

      !--------------------------------------------------------------
      !  Create vectors.  Note that we form 1 vector from scratch and
      !  then duplicate as needed.
!      CALL MatCreateVecs(PetA,PETSC_NULL_OBJECT,Petb,ierr)
      CALL MatCreateVecs(PetA,PETSC_NULL_VEC,Petb,ierr)
      CALL VecSetFromOptions(Petb,ierr)
      CALL VecDuplicate(Petb,Petx,ierr)
!      CALL VecDuplicate(Petb,DD,ierr)

      CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
!      CALL KSPGetPC(ksp,pc,ierr)
!      CALL PCSetType(pc,PCLU,ierr)
!      CALL KSPSetFromOptions(ksp,ierr)
!      CALL KSPSetOptionsPrefix(ksp,"kspu_",ierr)
!      CALL PCSetOptionsPrefix(pcu,"pcu_",ierr)

      CALL MatGetOwnershipRange(PetA,IStart,IEnd,ierr)
      print*,'Processor',rank,'starts at',IStart,'and ends at',IEnd-1

      !-----------------------------------------------
      ! define which elements are "owned" by which CPU
      !-----------------------------------------------
      do 157 j=1,ndpn
        do 157 i=1,nnds
          cpu(i,j)=ndofs+1 !initialize all to maxdof+1
 157  continue

      if (sizeof.eq.1) then
        do 158 i=1,neles
          eleproc(i)=.TRUE. !if 1 processor, it owns all elements
 158    continue
      else
        do 159 i=1,neles
          eleproc(i)=.FALSE. !assume processor doesn't own this element
          do 160 j=1,nnds
            do 160 k=1,ndpn
              !find gdof for all the dofs in this element
              gdof1=gcon(ele(i,j),k)
              if (gdof1.le.ndofs) cpu(j,k)=gdof1
 160      continue
          dof=MINVAL(cpu)-1 !single value that indicates MINIMUM gdof in element
          !if the minimum gdof is in the range of Istart->IEnd-1 then
          !this processor owns this element
          if ((dof.ge.IStart).and.(dof.lt.IEnd)) then
            eleproc(i)=.TRUE.
          endif
 159    continue
      endif

      !check single processor case has all true
      do 180 i=1,neles
          if ((eleproc(i).eqv..FALSE.).and.(sizeof.eq.1)) then
            print*,eleproc(i)
            print*,"ERROR"
            stop
          endif
 180  continue

      if (rank.eq.0) print*,"ABOUT TO CALL CALCB"
      call calcB(ndim,nnds,nints,maxnd,maxele,neles,ele,
     1           node,det,B,NN,rank)
 
      CALL MPI_Barrier(PETSC_COMM_WORLD,ierr)
      !CALL PetscBarrier(PETSC_NULL_INTEGER,ierr)
      ! options for initial condition
      ! 0: use initial conditions
      ! 1: pure paraelectric (use es and theta_init)
      ! 2: Purely z-Poled ferroelectric
      opt=0


      if (opt.eq.0) then

        if (rank.eq.0) then
          print*,"READING IN initial state"
          print*
        endif

        open(unit=24,file='initial')
        do 25 i=1,nodes
            read(24,*) u(i,1),u(i,2),u(i,3),phi(i),P(3,i,1),P(3,i,2),
     1             P(3,i,3),theta(i),P(2,i,1),P(2,i,2),P(2,i,3)
!            P(2,i,1)=P(3,i,1)
!            P(2,i,2)=P(3,i,2)
            P(1,i,1)=P(2,i,1)
            P(1,i,2)=P(2,i,2)
            P(1,i,3)=P(2,i,3)
 25     continue
        close(unit=24)
      else
        if (rank.eq.0) then
          print*
          print*,"GENERATING initial state"
        endif

        call initialState(opt,w,Pz,es,exx,eyy,ezz,theta_init,0.d0,
     1           node,u,phi,P,theta,maxnd,ndim,nodes,rank)

        if (rank.eq.0) then
          print*,"done with initial state"
          print*
        endif

      endif

      ! initialize f() and drhs()
      do 30 i=1,nodes*ndpn
        f(i)=0.d0
        drhs(i)=0.d0
 30   continue

      ! initialize unit boundary conditions
      do 556 i=1,ndbcs
        unitdbc(i)=dbcval(i)
 556  continue
      do 557 i=1,nfbcs
        unitfbc(i)=fbcval(i)
 557  continue
      do 558 i=1,ncon
        do 559 j=1,mpcnum(i)
          unitconB(i,j)=conB(i,j)
 559    continue
 558  continue

***** START GLOBAL IBCS LOOP
      minibcs=1
      maxibcs=3
      resultfilecount=1
      do 777 ibcs=minibcs,maxibcs
        flagT=0

        if (rank.eq.0) print*,'ibcs=',ibcs

        !phiapp=dble(ibcs-1)*1.d0
        !if (ibcs.lt.4) phiapp=dble(ibcs-1)
        !if (ibcs.ge.4) phiapp=2.d0+dble(ibcs-3)*2.d0
        phiapp=1.d0
        Qapp=1.d0
        uapp=1.d0
        Fapp=1.d0
        Tapp=1.d0
!        Tapp=0.97d0-0.01*(dble(ibcs-1))

        if (rank.eq.0) then
          print*
          print*,'*****************************'
          print*,'uapp   =',uapp
          print*,'Qapp   =',Qapp
          print*,'Tapp   =',Tapp
          print*,'phiapp =',phiapp
          print*,'*****************************'
        endif

        do 778 i=1,nfbcs
          if (fbcdof(i).eq.ndim+1) then
            fbcval(i)=unitfbc(i)*Qapp
          elseif (fbcdof(i).le.ndim) then
            fbcval(i)=unitfbc(i)*Fapp
          endif
 778    continue

        do 779 i=1,ndbcs
          if (dbcdof(i).le.ndim) then
            dbcval(i)=unitdbc(i)*uapp
          else
            if (dbcdof(i).eq.ndim+1) then
              dbcval(i)=unitdbc(i)*phiapp
            elseif (dbcdof(i).eq.2*ndim+2) then
              dbcval(i)=unitdbc(i)*Tapp
            else
              dbcval(i)=unitdbc(i)
            endif
          endif
 779    continue

******* CONTROL LOOP (TIME)
        ssu1=1.d0
        t=0.d0 ! start time at zero
        ic=301 ! max time steps for 1 ibc
        controlstart=0
        !if (ibcs.eq.18) controlstart=11

        do 900 control=controlstart,ic

          if (rank.eq.0) print*,'control (timestep)=',control,
     1                 ', ibcs=',ibcs,', resultfile =',resultfilecount
          flagNR=0

          !for isothermal sections
          if (1.eq.1) then
            if (control.le.0.and.ibcs.eq.1) then
              beta_P=1.d2
            elseif (control.gt.0.and.ibcs.eq.1) then
              beta_P=1.d1
            elseif (ibcs.eq.2) then
              beta_P=1.d-3
            elseif (ibcs.eq.3) then
              beta_P=0.d-1
            elseif (control.le.4.and.ibcs.gt.5) then
              beta_P=1.d-1
            elseif (control.le.12.and.ibcs.gt.5) then
              beta_P=5.d-1
            elseif (control.le.200.and.ibcs.gt.5) then
              beta_P=1.d-1
            endif
          endif
          !ar40: 14,1.d-1  87,1.d1  113,5,d0  250,1.d0
          !ar20: 12,1.d-1 rest 1.d1

          !coolin AR1
          if (1.eq.0) then
          if (control.lt.20) beta_P=1.d-1
          if (control.ge.20) beta_P=1.d1
          if (control.ge.630) beta_P=1.d2
          if (control.ge.634) beta_P=1.d1
          if (control.ge.665) beta_P=1.d0
          if (ibcs.ge.2) beta_P=1.d-1
          if (ibcs.eq.18) beta_P=1.d0
          if (ibcs.eq.18.and.control.ge.10) beta_P=1.d-1
          if (ibcs.ge.19) beta_P=1.d0
          if (ibcs.eq.22) beta_P=1.d1
          if (ibcs.eq.22.and.control.ge.2) beta_P=1.d0
          if (ibcs.eq.22.and.control.ge.27) beta_P=1.d-1
          endif

          !coolin AR10
          if (1.eq.0) then
          if (control.lt.20) beta_P=1.d-1
          if (control.ge.20) beta_P=1.d1
          if (control.ge.630) beta_P=1.d2
          if (control.ge.634) beta_P=1.d1
          if (control.ge.665) beta_P=1.d0
          if (ibcs.ge.2) beta_P=1.d-1
          if (ibcs.eq.18) beta_P=1.d0
          if (ibcs.eq.18.and.control.ge.10) beta_P=1.d-1
          if (ibcs.ge.19) beta_P=1.d0
          if (ibcs.eq.22) beta_P=1.d1
          if (ibcs.eq.22.and.control.ge.2) beta_P=1.d0
          if (ibcs.eq.22.and.control.ge.27) beta_P=1.d-1
          endif

          ! electrode AR10
          if (1.eq.0) then
          if (ibcs.eq.1.and.control.le.1) beta_P=1.d1
          if (ibcs.eq.1.and.control.ge.2) beta_P=1.d-1
          if (ibcs.eq.2.and.control.le.1) beta_P=1.d2
          if (ibcs.eq.2.and.control.ge.2) beta_P=1.d-1
          if (ibcs.eq.2.and.control.ge.58) beta_P=1.d-2
          if (ibcs.eq.3.and.control.le.1) beta_P=1.d1
          if (ibcs.eq.3.and.control.ge.2) beta_P=1.d-1
          if (ibcs.eq.3.and.control.ge.8) beta_P=1.d-2
          if (ibcs.eq.4.and.control.le.0) beta_P=1.d1
          if (ibcs.eq.4.and.control.ge.1) beta_P=1.d-1
          if (ibcs.eq.5.and.control.le.80) beta_P=1.d1
          if (ibcs.eq.5.and.control.ge.81) beta_P=1.d-1
          if (ibcs.ge.6.and.control.le.1) beta_P=1.d1
          if (ibcs.ge.6.and.control.ge.2) beta_P=1.d-1
          if (ibcs.eq.18.and.control.le.1) beta_P=1.d4
          if (ibcs.eq.18.and.control.ge.2) beta_P=1.d2
          if (ibcs.eq.18.and.control.ge.11) beta_P=1.d0
          if (ibcs.eq.19.and.control.le.1) beta_P=1.d2
          if (ibcs.eq.19.and.control.ge.2) beta_P=1.d0
          if (ibcs.ge.20.and.control.le.0) beta_P=1.d4
          if (ibcs.ge.20.and.control.ge.1) beta_P=1.d3
          if (ibcs.ge.20.and.control.ge.3) beta_P=1.d2
          if (ibcs.ge.20.and.control.ge.4) beta_P=1.d0
          if (ibcs.ge.21.and.control.le.1) beta_P=1.d3
          if (ibcs.ge.21.and.control.ge.2) beta_P=1.d1
          if (ibcs.ge.21.and.control.ge.3) beta_P=1.d0
          endif

          !cooling AR5
          if (1.eq.0) then
            if (control.le.22) beta_P=1.d-1
            if (control.ge.23) beta_P=1.d1
            if (control.ge.36) beta_P=1.d2
            if (control.ge.52) beta_P=1.d0
            if (ibcs.ge.9) beta_P=1.d-1
            if (ibcs.ge.9.and.control.ge.32) beta_P=1.d0
          endif

          if (rank.eq.0) print 742,'beta =',beta_P

 742      format(A9,E10.4E1,E10.4E1)
         ! reset rhs for next set of NR calcs
          do 81 ii=1,nodes*ndpn
            drhs(ii)=0.d0
 81      continue

          ! P(1,_,_) holds values from 2 time steps back
          ! P(2,_,_) holds values from 1 time step back
          if (control.ne.controlstart) then
            do 50 i=1,nodes
              P(1,i,1)=P(2,i,1)
              P(1,i,2)=P(2,i,2)
              P(1,i,3)=P(2,i,3)
              P(2,i,1)=P(3,i,1)
              P(2,i,2)=P(3,i,2)
              P(2,i,3)=P(3,i,3)
  50        continue
          endif

          ! compute grads & stresses for external force calcs
          if ((control.eq.controlstart).and.(ibcs.eq.minibcs)) then
            call calcdofs(ndim,ndpn,maxnd,maxdbc,maxcon,nodes,ndofs,
     1                 ncon,gcon,dbcval,drhs,u,phi,P,theta,
     2                 mpc,ndbcs,dbcnd,dbcdof,conA,conB,mpcnd,mpcnum)
            call gradients(ndim,nnds,nints,maxele,maxnd,neles,ele,B,
     1                    NN,u,strn,phi,E,P,Pint,PintPrev,gradP,
     2                    theta,thint,gradTheta)
            call stresses(ndim,nints,maxele,neles,strn,E,Pint,gradP,
     1                 thint,gradTheta,EC0cop1,strs,D,psiP,ksi,s,qt,a0,
     2                 beta_P,Psi,det,area,avg_stress,avg_E,
     3                 avg_D,avg_S,avg_T,avgTtot,dt,svm)

            ! write initial state as step "0"
            if (1.eq.1) then
              if ((control.eq.0).and.(ibcs.eq.1)) then
              call nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     1                    ele,E,D,s,strn,strs,svm,gradTheta,qt,nE,nD,nS,
     2                    nstrn,nstrs,ngradT,nflux,nsvm,nmat,EC0cop1)
                ! only have 1 cpu print initial state to file
                if (rank.eq.0) then
                call writeresults(fnum,ou,0,ndim,maxnd,maxele,nints,
     1                nnds,nodes,neles,ele,node,u,phi,P,nE,nD,nS,nstrn,
     2                nstrs,ngradT,nflux,beta_P,exx,eyy,
     3                avg_stress,avg_E,avg_D,avg_S,avg_T,PsiTot,theta,
     4                nsvm,EC0cop1,nmat)
                print*,'finished writing initial state'
                endif
                CALL MPI_Barrier(PETSC_COMM_WORLD,ierr)
                !CALL PetscBarrier(PETSC_NULL_VEC,ierr)
              endif
            endif
          endif

          call applybc(control,ndim,nnds,nints,ndpn,maxfc,fbcnd,
     1                 fbcdof,maxnd,nodes,nfbcs,gcon,fbcval,f,dt,
     2                 beta_P,maxele,neles,ele,det,NN,B,T0,strs,
     3                 D,psiP,ksi,s,Pint,PintPrev,thint,node,EC0cop1)
          ! modifies f using fbcs & prev stresses


********* BEGIN NR LOOP
          ! if maxsol reached, exit without finding solution
          do 100 isol=1,maxsol

            do 120 i=1,nodes*ndpn
              drhs0(i)=drhs(i)
 120        continue

            call calcdofs(ndim,ndpn,maxnd,maxdbc,maxcon,nodes,ndofs,
     1                ncon,gcon,dbcval,drhs0,u,phi,P,theta,
     2                mpc,ndbcs,dbcnd,dbcdof,conA,conB,mpcnd,mpcnum)
            !modifies u,phi,P(3,..,...),theta

!            do i=1,nodes
!              print*,i,u(i,1),u(i,2),u(i,3),phi(i),P(3,i,1),P(3,i,2),
!     1        P(3,i,3),theta(i)
!            enddo
            call gradients(ndim,nnds,nints,maxele,maxnd,neles,ele,B,
     1                     NN,u,strn,phi,E,P,Pint,PintPrev,gradP,
     2                     theta,thint,gradTheta)
            !modifies E,Pint,gradP,strn,thint,gradTheta

            call stresses(ndim,nints,maxele,neles,strn,E,Pint,gradP,
     1                 thint,gradTheta,EC0cop1,strs,D,psiP,ksi,s,qt,a0,
     2                 beta_P,Psi,det,area,avg_stress,avg_E,
     3                 avg_D,avg_S,avg_T,avgTtot,dt,svm)
            !modifies strs,s,D,phiP,ksi
            !computes Psi,avg_stress,avg_E,avg_D

            PsiPrint=0
            call intforce(ndim,nints,nnds,ndpn,maxnd,nodes,ny,maxele,
     1                    neles,ele,gcon,det,B,NN,strs,D,psiP,ksi,s,qt,
     2                    finte,dt,beta_P,Pint,PintPrev,Psi,
     3                    PsiTot,PsiPrint,Sprint,area,h,node,control,
     4                    EC0cop1)
            !finte now populated with values from updated dofs, finte(t+dt)

            do 90 i=1,nodes*ndpn
              drhs(i)=f(i)-finte(i)
              if (finte(i).gt.1.d20) print*,i
 90         continue

            do 96 i=1,ncon
              j=ndofs+ndbcs+i
              do 97 jj=1,mpcnum(i)
                k=gcon(mpc(i,jj,1),mpc(i,jj,2))
                drhs(k)=drhs(k)+drhs(j)*conA(i,jj)
97            continue
96          continue

            !Debugging loop
            if (1.eq.0) then
              do 111 i=1,nodes
                do 111 j=1,ndpn
                  gdof1=gcon(i,j)
                  if (((abs(f(gdof1)-finte(gdof1)).gt.1.d-4).and.
     1             (gdof1.le.ndofs)).and.(rank.eq.0.and.isol.gt.3)) then
                  print*,i,j,f(gdof1),finte(gdof1)
                  endif
 111          continue
            endif!

            sumsqf=0.d0
            do 110 i=1,ndofs
              sumsqf=sumsqf+drhs(i)*drhs(i)
 110        continue
            !-----------------------------------------------------------
            !----------location reordering for petsc--------------------

            ! write vector loc=[0,1,2,....]
            do 122 i=1,ndofs
              loc(i)=i-1
              Petdrhs(i)=drhs(i)
 122        continue

            !----------set vectors -------------------------------------
            ! initialize x to zero (from Ax=b)
            CALL VecSet(Petx,zero,ierr)

            ! set b to drhs, done by single processor
            if (rank.eq.0) then
              CALL VecSetValues(Petb,ndofs,loc,Petdrhs,
     1                               INSERT_VALUES,ierr)
            endif
            CALL VecAssemblyBegin(Petb,ierr)
            CALL VecAssemblyEnd(Petb,ierr)
            !CALL PetscBarrier(PETSC_NULL_INTEGER,ierr)
            CALL MPI_Barrier(PETSC_COMM_WORLD,ierr)
            call VecNorm(Petb,NORM_2,sumsqf,ierr)
            sumsqf=sumsqf**2
            !write drhs before solving
!            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"rhsb.m",
!     &                  rhsview,ierr)
!            call PetscViewerPushFormat(rhsview,
!     1                     PETSC_VIEWER_ASCII_MATLAB,ierr)
!            call VecView(Petb,rhsview,ierr)
!            call PetscViewerDestroy(rhsview)
            !end write

            !-----------------------------------------------------------
            !------- ASSEMBLY IN STF ROUTINE ---------------------------
!            if (rank.eq.0) print*,'entering stf'
            call stf(ndim,nnds,nints,ndpn,B,NN,det,
     1               maxele,maxnd,maxdbc,maxcon,ndofs,neles,
     2               ele,gcon,strn,mpc,conA,ndbcs,ncon,mpcnum,
     3               dt,beta_P,Pint,PintPrev,thint,a0,control,
     4               nodes,PetA,eleproc,Petx,Petb,rank,EC0cop1)
            !CALL PetscBarrier(PETSC_NULL_VEC,ierr)
            CALL MPI_Barrier(PETSC_COMM_WORLD,ierr)
!            if (rank.eq.0) print*,'done w/ stf for all processors'
           !-------------------------------------------------------------
           !-------------- SOLVER ---------------------------------------
           CALL KSPSetOperators(ksp,PetA,PetA,ierr)
           CALL KSPGetPC(ksp,pc,ierr)
           CALL PCSetType(pc,PCLU,ierr)
!          CALL KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
           CALL KSPSetFromOptions(ksp,ierr)

           CALL KSPsolve(ksp,Petb,Petx,ierr)
           CHKERRA(ierr)
           CALL KSPGetIterationNumber(ksp,its,ierr)
!           if (its.eq.9900) then
!             open(unit=447,file='KSPwarning_u.txt')
!             write(447,*) "reach 9900!!!!!"
!           close(unit=447)
!           endif
!           call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)

!           CALL KSPGetIterationNumber(ksp,its,ierr)
!           print*,'its=',its,'on processor',rank

           if (sizeof.eq.1) then
             CALL VecGetArrayRead(Petx,x_array,i_x,ierr)
             do i=1,ndofs
               drhs(i)=x_array(i_x+i)
             enddo
!            print*
!            print*,i_x,x_array(i_x+1),x_array(i_x+ndofs)
             CALL VecRestoreArrayRead(Petx,x_array,i_x,ierr)

           else
             allocate (xSEQ)
             CALL VecScatterCreateToAll(Petx,sca,xSEQ,ierr)
             CALL VecScatterBegin(sca,Petx,xSEQ
     1                     ,INSERT_VALUES,SCATTER_FORWARD,ierr)
             CALL VecScatterEnd(sca,Petx,xSEQ
     1                     ,INSERT_VALUES,SCATTER_FORWARD,ierr)
             CALL VecGetArrayRead(xSEQ,x_array,i_x,ierr)
             do i=1,ndofs
               drhs(i)=x_array(i_x+i)
             enddo
             CALL VecRestoreArrayRead(xSEQ,x_array,i_x,ierr)
             !CALL PetscBarrier(PETSC_NULL_VEC,ierr)
             CALL MPI_Barrier(PETSC_COMM_WORLD,ierr)
             CALL VecDestroy(xSEQ,ierr)
             CALL VecScatterDestroy(sca,ierr)
           endif

          !-------------------------------------------------------------
          !------- PETSC SOLVER IS DONE! -------------------------------
          !-------------------------------------------------------------

           !call sols(ndpn,ndofs,maxwk,maxnd,maxa,A,drhs)
!      rhsfile=2
!      open(unit=rhsfile,file='rhs2.txt')

!      do 222 i=1,nodes*ndpn
!        write(rhsfile,*) drhs(i)
! 222  continue
!      close(unit=rhsfile)
!      stop
           ! find sum of squares of displacement residual
           sumsqu=0.d0
           do 130 i=1,ndofs
             sumsqu=sumsqu+drhs(i)*drhs(i)
             !print*,i,drhs(i)
             if (isnan(sumsqu).and.rank.eq.0) then
              print*,'stopping: problem with dof',i,'  drhs(i)=',drhs(i)
              stop
             endif
             !if (i.eq.1.and.rank.eq.0) print*,i,sumsqu
 130       continue
           if (isol.eq.1) ssu1=sumsqu

           if (rank.eq.0) then
             print 789, 'isol=',isol,'sumsqf=',sumsqf,
     1           'sumsqu=',sumsqu,'Savg=',avg_S,'Tavg=',avg_T,
     2           'py',P(3,85,2),phi(1052)
           endif
 789       format (A7,I4,'  ',A7,E10.2E2,'  ',A7,E10.2E2,'  ',
     1          A7,F10.3,'  ',A7,F10.8,'  ',A7,F10.6,'  ',F10.2)

           if (isnan(sumsqf)) stop
           if (isnan(sumsqu)) stop

           ! if converged, exit nr loop, if diverged, stop
           if (sumsqf.lt.1.d-13.and.sumsqu.lt.1.d-13) then
             if ((isol.eq.1)) then
               if (rank.eq.0) then
                 print*,"converged in time, continue to next ibc"
               endif
               flagT=1 !flagT indicated convergence in control loop
               !reset w/ every new ibc
             endif
             flagNR=1 !flag NR indicates convergence in NR loop
             !rest w/ every new time step
             goto 200
           endif
           if (sumsqf.gt.1.d40.or.sumsqu.gt.1.d40) then
             if (rank.eq.0) print*,'Solution diverging'
             stop
           endif

 100     continue !end NR loop
         if (rank.eq.0) print*,'Convergence not achieved'
         stop

 200     continue !go here if converged in NR loop

         !!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if ((mod(ibcs,1).eq.0).and.(rank.eq.0)) then
           if (flagT.eq.1) then
             call calcdofs(ndim,ndpn,maxnd,maxdbc,maxcon,nodes,ndofs,
     1                ncon,gcon,dbcval,drhs0,u,phi,P,theta,
     2                mpc,ndbcs,dbcnd,dbcdof,conA,conB,mpcnd,mpcnum)

             call gradients(ndim,nnds,nints,maxele,maxnd,neles,ele,B,
     1                     NN,u,strn,phi,E,P,Pint,PintPrev,gradP,
     2                     theta,thint,gradTheta)

             call stresses(ndim,nints,maxele,neles,strn,E,Pint,gradP,
     1                 thint,gradTheta,EC0cop1,strs,D,psiP,ksi,s,qt,a0,
     2                 beta_P,Psi,det,area,avg_stress,avg_E,
     3                 avg_D,avg_S,avg_T,avgTtot,dt,svm)

             ! average y strain in whole model
             avgstrain=0.d0
             do 280 j=1,neles
               do 280 k=1,nints
                 avgstrain=avgstrain+det(j,k)*strn(j,k,2,2)
 280         continue
             avgstrain=avgstrain/(node(nodes,1)-node(1,1))
     1                       /(node(nodes,2)-node(1,2))

             Qapp=0.d0
             do 1003 i=1,ndbcs
               if (dbcdof(i).eq.3) then
                 if (dabs(unitdbc(i)).gt.1.d-4) then
                   Qapp=Qapp+finte(gcon(dbcnd(i),3))
                 endif
               endif
 1003        continue
           endif

           print*,Qapp

           if (1.eq.1) then

           if (1.eq.1) then
             if (control.le.20) then
               modwrite=1
             elseif (control.le.100) then
               modwrite=1
             elseif (control.le.1000) then
               modwrite=15
             elseif (control.le.10000) then
               modwrite=500
             endif
           else
             modwrite=1
           endif

!                     (mod(control,modwrite).eq.0)) then
           if (((flagNR.eq.1)).and.(mod(control,modwrite).eq.0)) then
!     1              .and.((ssu1.lt.1.d-15).or.(control.eq.ic))) then

             call nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     1                    ele,E,D,s,strn,strs,svm,gradTheta,qt,nE,nD,nS,
     2                    nstrn,nstrs,ngradT,nflux,nsvm,nmat,EC0cop1)

             !!!!!!!!!!! WRITING RESULTS TO FILE !!!!!!!!!!!!!!!!!!!!!!!
             !nint(dble(control/1.d2))+1,(ibcs-1)*6+control+1
             !resultfilecount
             print*,"writing results"
             call writeresults(fnum,ou,resultfilecount,ndim,
     1       maxnd,maxele,nints,nnds,nodes,neles,ele,node,u,phi,P,nE,nD,
     2                nS,nstrn,nstrs,ngradT,nflux,beta_P,exx,eyy,
     3                avg_stress,avg_E,avg_D,avg_S,avg_T,PsiTot,theta,
     4                nsvm,EC0cop1,nmat)
             resultfilecount=resultfilecount+1

           endif
           endif
         endif

         !if we are taken enough iterations in time loop
         !move on to next ibc
         if (ssu1.lt.1.d-15) then
           if (rank.eq.0) print*,"going on"
           flagT=1
         endif
         if (flagT.eq.1) goto 777

 900   continue ! end control/time loop
 777  continue ! end ibc loop

      CALL VecDestroy(Petb,ierr)
      CALL VecDestroy(Petx,ierr)
      call MatDestroy(PetA,ierr)
      call KSPDestroy(ksp,ierr)
      !call PCDestroy(pc,ierr)

      call PetscFinalize(ierr)
      CHKERRA(ierr)

      end program

************************************************************************
****** END OF PROGRAM **************************************************
****** BEGIN SUBROUTINES ***********************************************
************************************************************************

      subroutine readinp(ndim,nnds,ndpn,maxnd,maxele,maxfc,maxdbc,
     1             maxcon,nodes,neles,ndbcs,nfbcs,ndofs,gcon,ele,
     2             fbcnd,fbcdof,dbcnd,dbcdof,node,fbcval,dbcval,
     3             mpc,conA,conB,mpcnum,ncon,mpcnd,nx,ny,rank)

      real*8 node(maxnd,ndim),dv,dbcval(maxdbc),fbcval(maxfc)
      real*8 conA(maxdbc,maxcon),conB(maxdbc,maxcon)
      integer fbcnd(maxfc),fbcdof(maxfc),ndofs,mpc(maxdbc,maxcon,2)
      integer dbcnd(maxdbc),dbcdof(maxdbc)
      integer nx,ny,mpcnum(maxdbc)
      integer i,j,k,nodes,neles,ndbcs,nfbcs,ele(maxele,8)
      integer maxnd,maxele,maxdbc,maxfc,n,d,inds,iele,idbc,ifbc
      integer ncon,num,mpcnd(maxdbc,2),gcon(maxnd,ndpn),rank

      inds=15
      iele=16
      idbc=20
      ifbc=19

      open(unit=inds,file='nodes')
      read(inds,*)nodes
      if (nodes.gt.maxnd) then
         print*,'too many nodes'
         print*,'nodes=',nodes,'maxnd=',maxnd
         stop
      endif

      do 250 i=1,nodes
        read(inds,*)k,(node(i,j),j=1,ndim)
 250  continue
      close(unit=inds)
      ndofs=nodes*ndpn

!      print*,'Nodes:'
!      do 222 i=1,nodes
!        print*,node(i,1),node(i,2)
! 222  continue

      open(unit=iele,file='elements')
      read(iele,*)neles
      if (neles.gt.maxele) then
         print*,'too many elements'
         print*,'neles',neles,'maxele',maxele
         stop
      endif
      do 300 i=1,neles
        read(iele,*)k,(ele(i,j),j=1,nnds)
 300  continue
      close(unit=iele)

!      do i=1,neles
!        print*,ele(i,1),ele(i,2),ele(i,3),ele(i,4),ele(i,5),ele(i,6),
!     1            ele(i,7),ele(i,8)
!      enddo

      do 390 i=1,nodes
        do 390 j=1,ndpn
          gcon(i,j)=(i-1)*ndpn+j
 390  continue

      open(unit=idbc,file='dbcs')
      read(idbc,*)ndbcs,ncon
      if (ndbcs.gt.maxdbc) then
        print*,'too many dbcs'
        stop
      endif
      if (ncon.gt.maxdbc) then
        print*,'too many multi-point constraints'
        stop
      endif

      do 400 i=1,ndbcs
        read(idbc,*)n,d,dv
        dbcnd(i)=n
        dbcdof(i)=d
        dbcval(i)=dv
        num=gcon(n,d)

        do 410 j=1,nodes
          do 420 k=1,ndpn
            if((gcon(j,k)).gt.num) gcon(j,k)=gcon(j,k)-1
 420      continue
 410    continue
        gcon(n,d)=nodes*ndpn
        ndofs=ndofs-1
 400  continue
!      print*,'ndofs = ',ndofs

      do 425 i=1,ncon
        read(idbc,*)mpcnum(i),mpcnd(i,1),mpcnd(i,2)
        do 426 j=1,mpcnum(i)
          read(idbc,*)mpc(i,j,1),mpc(i,j,2),conA(i,j),conB(i,j)
 426    continue
        num=gcon(mpcnd(i,1),mpcnd(i,2))
        do 427 j=1,nodes
          do 428 k=1,ndpn
            if((gcon(j,k)).gt.num) gcon(j,k)=gcon(j,k)-1
 428      continue
 427    continue
        gcon(mpcnd(i,1),mpcnd(i,2))=nodes*ndpn
        ndofs=ndofs-1
 425  continue

      if (rank.eq.0) then
        print*
        print*,'*******************************************************'
        print*,'dimension of stiffness matrix:  ',nodes*ndpn
        print*
        print*,'ndofs = ',ndofs
        print*,'ndbcs = ',ndbcs
        print*,'ncon  = ',ncon
        print*
        print*,'nodes*ndpn-ndbcs      = ',nodes*ndpn-ndbcs
        print*,'nodes*ndpn-ncon       = ',nodes*ndpn-ncon
        print*,'nodes*ndpn-ndbcs-ncon = ',nodes*ndpn-ndbcs-ncon
        print*,'*******************************************************'
        print*
      endif

      close(unit=idbc)

!      do 501 i=1,nodes
!      print*,gcon(i,1),gcon(i,2),gcon(i,3),gcon(i,4),
!     1         gcon(i,5),gcon(i,6),gcon(i,7),gcon(i,8)
! 501  continue

      open(unit=ifbc,file='fbcs')
      read(ifbc,*)nfbcs
      if (nfbcs.gt.maxfc) then
        print*,'too many force bcs'
        stop
      endif
      do 500 i=1,nfbcs
        read(ifbc,*)n,d,dv
        fbcnd(i)=n
        fbcdof(i)=d
        fbcval(i)=dv
 500  continue
      close(unit=ifbc)

      return
      end

***********************************************************************

      subroutine choosemat(ndim,nnds,maxnd,maxele,neles,ele,
     1             node,elthick,h,w,EC0cop1)

      real*8 node(maxnd,ndim),elthick,xi,yi,h,w
      integer i,n,neles,ele(maxele,8)
      integer maxnd,maxele,nnds
      integer EC0cop1(maxele)

      ! in array EC0cop1, 0= electrocaloric material, 1=electrode
      do 100 n=1,neles
        xi=0.d0
        yi=0.d0
        ! find an interior point
        do 200 i=1,nnds
          xi=xi+node(ele(n,i),1)
          yi=yi+node(ele(n,i),2)
 200    continue
        xi=xi/dble(nnds)
        yi=yi/dble(nnds)

        if (xi.lt.elthick.or.xi.gt.w+elthick) then
          EC0cop1(n)=1
        elseif (yi.lt.0.75d0*h+elthick/2.d0.and.
     1                   yi.gt.0.75d0*h-elthick/2.d0) then
          if ((xi.lt.0.95d0*w+elthick).or.
     1             ((yi-0.75d0*h)**2.d0+(xi-0.95d0*w-elthick)**2.d0.le.
     2                          (elthick/2.d0)**2.d0)) then
            EC0cop1(n)=1
          endif
        elseif (yi.lt.0.25d0*h+elthick/2.d0.and.
     1                   yi.gt.0.25d0*h-elthick/2.d0) then
          if ((xi.gt.0.05d0*w+elthick).or.
     1             ((yi-0.25d0*h)**2.d0+(xi-0.05d0*w-elthick)**2.d0.le.
     2                          (elthick/2.d0)**2.d0)) then
            EC0cop1(n)=1
          endif
        else
          EC0cop1(n)=0
        endif

 100  continue

      return
      end

***********************************************************************

      subroutine choosematFillet(ndim,nnds,maxnd,maxele,neles,ele,
     1             node,elthick,h,w,EC0cop1)

      real*8 node(maxnd,ndim),elthick,xi,yi,h,w
      integer i,n,neles,ele(maxele,8)
      integer maxnd,maxele,nnds
      integer EC0cop1(maxele)

      ! in array EC0cop1, 0= electrocaloric material, 1=electrode
      do 100 n=1,neles

        xi=0.d0
        yi=0.d0
        ! find an interior point
        do 200 i=1,nnds
          xi=xi+node(ele(n,i),1)
          yi=yi+node(ele(n,i),2)
 200    continue
        xi=xi/dble(nnds)
        yi=yi/dble(nnds)

        !on either rhs or lhs edge
        if (xi.lt.elthick.or.xi.gt.w+elthick) then
          EC0cop1(n)=1
        !between 3h/4 +/- elthick/2
        elseif (yi.lt.0.75d0*h+elthick/2.d0.and.
     1                   yi.gt.0.75d0*h-elthick/2.d0) then
          if ((xi.lt.0.95d0*w+elthick).or.
     1             ((yi-0.75d0*h)**2.d0+(xi-0.95d0*w-elthick)**2.d0.le.
     2                          (elthick/2.d0)**2.d0)) then
            EC0cop1(n)=1
          endif
        !between h/4 +/- elthick/2
        elseif (yi.lt.0.25d0*h+elthick/2.d0.and.
     1                   yi.gt.0.25d0*h-elthick/2.d0) then
          if ((xi.gt.0.05d0*w+elthick).or.
     1             ((yi-0.25d0*h)**2.d0+(xi-0.05d0*w-elthick)**2.d0.le.
     2                          (elthick/2.d0)**2.d0)) then
            EC0cop1(n)=1
          endif
        !in main EC portion
        else
          EC0cop1(n)=0
        endif

        if (xi.lt.2.d0*elthick.and.xi.gt.elthick) then
          if (((yi.lt.0.75d0*h+1.5d0*elthick
     1        .and.yi.gt.0.75d0*h+0.5d0*elthick)
     2      .and.((yi-(0.75*h+1.5d0*elthick))**2.d0+
     3         (xi-(2.d0*elthick))**2.d0.gt.elthick**2.d0))
     4    .or.((yi.lt.0.75d0*h-0.5d0*elthick
     5        .and.yi.gt.0.75d0*h-1.5d0*elthick)
     6      .and.((yi-(0.75*h-1.5d0*elthick))**2.d0+
     7        (xi-(2.d0*elthick))**2.d0.gt.elthick**2.d0))) then

            EC0cop1(n)=1
          endif

        elseif (xi.lt.w+elthick.and.xi.gt.w) then
          if (((yi.lt.0.25d0*h+1.5d0*elthick
     1          .and.yi.gt.0.25d0*h+0.5d0*elthick)
     2        .and.((yi-(0.25*h+1.5d0*elthick))**2.d0+
     3           (xi-w)**2.d0.gt.elthick**2.d0))
     4      .or.((yi.lt.0.25d0*h-0.5d0*elthick
     5         .and.yi.gt.0.25d0*h-1.5d0*elthick)
     6        .and.((yi-(0.25*h-1.5d0*elthick))**2.d0+
     7           (xi-w)**2.d0.gt.elthick**2.d0))) then

            EC0cop1(n)=1
          endif
        endif

 100  continue

      return
      end

*************************************************************


*************************************************************

*************************************************************
*     This subroutine will construct the bele(4,4,2) for    *
*     each element and store them in B(maxele,4,4,2).       *
*************************************************************
      subroutine calcB(ndim,nnds,nints,maxnd,maxele,neles,ele,
     1                 node,det,B,NN,rank)
      real*8 node(maxnd,ndim),B(maxele,nnds,nints,ndim)
      real*8 det(maxele,nints),detp,ted
      real*8 rst(2),r,s,t,p(nnds,ndim),aj(ndim,ndim),aji(ndim,ndim)
      real*8 NN(nnds,nints),A(ndim,ndim)
      integer ndim,nnds,nints,nc,rank
      integer i,j,k,m,ele(maxele,8),neles,maxnd,maxele
      integer k1,k2,k3,k123

      !cube root of nints gives # of integration stations / dim
      nc=idint(dble(nints)**(1/3.d0))

      if (nints.eq.8) then
        rst(1)=-1.d0/dsqrt(3.d0)
        rst(2)=-rst(1)
      else
        print*,'Number of integration stations is not supported'
      endif

      if (nnds.ne.8) then
        print*,'Number of nodes is not supported'
        stop
      endif

      do 10 k1=1,nc
        r=rst(k1)
        do 20 k2=1,nc
          s=rst(k2)
          do 30 k3=1,nc
            t=rst(k3)

            k123=(k1-1)*nc*nc+(k2-1)*nc+k3

            if (nnds.eq.8) then
              NN(1,k123)=(1-r)*(1-s)*(1-t)/8.d0
              NN(2,k123)=(1+r)*(1-s)*(1-t)/8.d0
              NN(3,k123)=(1+r)*(1+s)*(1-t)/8.d0
              NN(4,k123)=(1-r)*(1+s)*(1-t)/8.d0
              NN(5,k123)=(1-r)*(1-s)*(1+t)/8.d0
              NN(6,k123)=(1+r)*(1-s)*(1+t)/8.d0
              NN(7,k123)=(1+r)*(1+s)*(1+t)/8.d0
              NN(8,k123)=(1-r)*(1+s)*(1+t)/8.d0
            endif

 30       continue
 20     continue
 10   continue

      do 100 m=1,neles

        do 310 k1=1,nc
          r=rst(k1)
          do 320 k2=1,nc
            s=rst(k2)
            do 330 k3=1,nc
              t=rst(k3)

              k123=(k1-1)*nc*nc+(k2-1)*nc+k3

              if (nnds.eq.8) then
                p(1,1)=-(1-s)*(1-t)/8.d0
                p(2,1)=(1-s)*(1-t)/8.d0
                p(3,1)=(1+s)*(1-t)/8.d0
                p(4,1)=-(1+s)*(1-t)/8.d0
                p(5,1)=-(1-s)*(1+t)/8.d0
                p(6,1)=(1-s)*(1+t)/8.d0
                p(7,1)=(1+s)*(1+t)/8.d0
                p(8,1)=-(1+s)*(1+t)/8.d0

                p(1,2)=-(1-r)*(1-t)/8.d0
                p(2,2)=-(1+r)*(1-t)/8.d0
                p(3,2)=(1+r)*(1-t)/8.d0
                p(4,2)=(1-r)*(1-t)/8.d0
                p(5,2)=-(1-r)*(1+t)/8.d0
                p(6,2)=-(1+r)*(1+t)/8.d0
                p(7,2)=(1+r)*(1+t)/8.d0
                p(8,2)=(1-r)*(1+t)/8.d0

                p(1,3)=-(1-r)*(1-s)/8.d0
                p(2,3)=-(1+r)*(1-s)/8.d0
                p(3,3)=-(1+r)*(1+s)/8.d0
                p(4,3)=-(1-r)*(1+s)/8.d0
                p(5,3)=(1-r)*(1-s)/8.d0
                p(6,3)=(1+r)*(1-s)/8.d0
                p(7,3)=(1+r)*(1+s)/8.d0
                p(8,3)=(1-r)*(1+s)/8.d0
              endif

              do 402 i=1,ndim
                do 410 j=1,ndim
                  aj(i,j)=0.d0
                  do 420 k=1,nnds
                    aj(i,j)=aj(i,j)+p(k,i)*node(ele(m,k),j)
 420              continue
 410            continue
 402          continue

              if (ndim.ne.3) then
                print*,'Number of dimensions is not supported'
                stop
              endif

              A(1,1)=aj(2,2)*aj(3,3)-aj(2,3)*aj(3,2)
              A(2,1)=aj(2,3)*aj(3,1)-aj(2,1)*aj(3,3)
              A(1,2)=aj(3,2)*aj(1,3)-aj(1,2)*aj(3,3)
              A(3,1)=aj(2,1)*aj(3,2)-aj(3,1)*aj(2,2)
              A(1,3)=aj(1,2)*aj(2,3)-aj(1,3)*aj(2,2)
              A(2,2)=aj(3,3)*aj(1,1)-aj(3,1)*aj(1,3)
              A(3,2)=aj(3,1)*aj(1,2)-aj(3,2)*aj(1,1)
              A(2,3)=aj(1,3)*aj(2,1)-aj(2,3)*aj(1,1)
              A(3,3)=aj(1,1)*aj(2,2)-aj(1,2)*aj(2,1)
              detp=aj(1,1)*(aj(2,2)*aj(3,3)-aj(3,2)*aj(2,3))
     1                 - aj(1,2)*(aj(2,1)*aj(3,3)-aj(3,1)*aj(2,3))
     2                 + aj(1,3)*(aj(2,1)*aj(3,2)-aj(3,1)*aj(2,2))

              ted=1.d0/detp

              do 678 i=1,ndim
                do 679 j=1,ndim
                  aji(i,j)=ted*A(i,j)
 679            continue
 678          continue

              do 500 i=1,ndim
                do 510 j=1,nnds
                  B(m,j,k123,i)=0.d0
                  do 520 k=1,ndim
                    B(m,j,k123,i)=B(m,j,k123,i)+aji(i,k)*p(j,k)
 520              continue
 510            continue
 500          continue

              det(m,k123)=detp

 330        continue
 320        continue
 310      continue

 100  continue
      if (rank.eq.0) print*,"end, calcB"

      return
      end

***********************************************************************
      subroutine applybc(control,ndim,nnds,nints,ndpn,maxfc,fbcnd,
     1                  fbcdof,maxnd,nodes,nfbcs,gcon,fbcval,f,dt,
     2                  beta_P,maxele,neles,ele,det,NN,B,T0,strs,
     3                  D,psiP,ksi,s,Pint,PintPrev,thint,node,EC0cop1)
      implicit none
      real*8 fbcval(maxfc),f(maxnd*ndpn)
      real*8 fele(nnds*ndpn),det(maxele,nints),ww(nints),dt
      real*8 beta_P,T0,Tcn
      real*8 NN(nnds,nints),B(maxele,nnds,nints,ndim)
      real*8 PintPrev(maxele,nints,ndim),Pint(maxele,nints,ndim)
      real*8 thint(maxele,nints)
      real*8 strs(maxele,nints,ndim+1,ndim+1),D(maxele,nints,ndim)
      real*8 ksi(maxele,nints,ndim,ndim),psiP(maxele,nints,ndim)
      real*8 s(maxele,nints),node(maxnd,ndim)
      integer maxfc,fbcnd(maxfc),control
      integer nfbcs,nodes,gcon(maxnd,ndpn),maxnd,fbcdof(maxfc)
      integer i,j,k,m,gdof,n,gp,nodeP
      integer ndim,nnds,nints,ndpn
      integer maxele,ele(maxele,8),neles
      integer EC0cop1(maxele)

      if (nints.eq.8) then
        do i=1,nints
          ww(i)=1.d0
        enddo
      endif

      Tcn=T0/T0

      do 50 i=1,nodes*ndpn
        f(i)=0.d0
 50   continue

      ! first step, time deriv of force is force/dt, otherwise change=0
      if (control.eq.0) then
        do 100 i=1,nfbcs
          gdof=gcon(fbcnd(i),fbcdof(i))
          f(gdof)=fbcval(i)/dt
 100    continue
      endif

      do 99 n=1,neles
        do 101 m=1,ndpn*nnds
          fele(m)=0.d0
 101    continue

        if (control.ne.0) then
          ! F(ux,uy)
          do 102 nodeP=1,nnds
            do 202 i=1,ndim
              m=(nodeP-1)*ndpn+i
              do 302 gp=1,nints
                do 402 j=1,ndim
                  fele(m)=fele(m)+strs(n,gp,j,i)*B(n,nodeP,gp,j)*
     1                                    det(n,gp)*ww(gp)/dt
 402            continue
 302          continue
 202        continue
 102      continue

        ! F(phi)
          do 103 nodeP=1,nnds
            m=(nodeP-1)*ndpn+ndim+1
            do 203 gp=1,nints
              do 303 i=1,ndim
                fele(m)=fele(m)+D(n,gp,i)*B(n,nodeP,gp,i)*
     1                          det(n,gp)*ww(gp)/dt
 303          continue
 203        continue
 103      continue

          ! F(Px,Py)
!          if (EC0cop1(n).eq.0) then
          do 104 nodeP=1,nnds
            do 204 i=1,ndim
              m=(nodeP-1)*ndpn+ndim+1+i

              do 304 gp=1,nints
                fele(m)=fele(m)+(psiP(n,gp,i)/dt+beta_P/(dt**2)*
     1                           (2.d0*Pint(n,gp,i)-PintPrev(n,gp,i)))*
     2                        NN(nodeP,gp)*det(n,gp)*ww(gp)
                do 404 j=1,ndim
                  fele(m)=fele(m)+ksi(n,gp,j,i)*B(n,nodeP,gp,j)*
     1                            det(n,gp)*ww(gp)/dt
 404            continue
 304          continue
 204        continue
 104      continue
!          endif
        else
!          if (EC0cop1(n).eq.0) then
          do 107 nodeP=1,nnds
            do 207 i=1,ndim
              m=(nodeP-1)*ndpn+ndim+1+i
              do 307 gp=1,nints
                   fele(m)=fele(m)+beta_P/(dt**2)*PintPrev(n,gp,i)*
     1                        NN(nodeP,gp)*det(n,gp)*ww(gp)
 307          continue
 207        continue
 107      continue
!          endif
        endif

        ! F(T)
        do 105 nodeP=1,nnds
          m=nodeP*ndpn

          do 205 gp=1,nints
            fele(m)=fele(m)-s(n,gp)*NN(nodeP,gp)*det(n,gp)*ww(gp)/dt

!            if (EC0cop1(n).eq.0) then
              do 305 i=1,ndim
                fele(m)=fele(m)-beta_P*PintPrev(n,gp,i)*PintPrev(n,gp,i)
     1                        *NN(nodeP,gp)*det(n,gp)*ww(gp)/(Tcn*dt**2)
 305          continue
!            endif
 205      continue
 105    continue

        do 106 nodeP=1,nnds
          do 206 k=1,ndpn
            m=(nodeP-1)*ndpn+k
            gdof=gcon(ele(n,nodeP),k)
            f(gdof)=f(gdof)+fele(m)
!            if (k.eq.3) print*,f(gdof),fele(m)
 206      continue
 106    continue

 99   continue

      return
      end
***********************************************************************

      subroutine calcdofs(ndim,ndpn,maxnd,maxdbc,maxcon,nodes,ndofs,
     1                 ncon,gcon,dbcval,drhs,u,phi,P,theta,
     2                 mpc,ndbcs,dbcnd,dbcdof,conA,conB,mpcnd,mpcnum)
      implicit none
      real*8 drhs(maxnd*ndpn)
      real*8 u(maxnd,ndim),phi(maxnd),P(3,maxnd,ndim),theta(maxnd)
      real*8 dbcval(maxdbc),conA(maxdbc,maxcon),conB(maxdbc,maxcon)
      integer maxnd,maxdbc,nodes,mpc(maxdbc,maxcon,2),dbcnd(maxdbc)
      integer i,j,gcon(maxnd,ndpn),gdof,ndof,ncon,ndofs,ndbcs,ndim,ndpn
      integer mpcnd(maxdbc,2),maxcon,mpcnum(maxdbc),dbcdof(maxdbc)

      do 100 i=1,nodes
        do 200 j=1,ndpn
          gdof=gcon(i,j) !global degree of freedom
          if (gdof.le.ndofs) then
            if (j.le.ndim) then
              u(i,j)=u(i,j)+drhs(gdof) !update u=u_prev+du
            elseif (j.eq.ndim+1) then
              phi(i)=phi(i)+drhs(gdof) !update phi...
            elseif (j.le.2*ndim+1) then
              P(3,i,j-ndim-1)=P(3,i,j-ndim-1)+drhs(gdof) !update P...
            else
              theta(i)=theta(i)+drhs(gdof) !update theta...

              if (theta(i).lt.0.d0) then
                print*,'theta<0',drhs(gdof),theta(i)
                stop
              endif

            endif
          endif
 200    continue
 100  continue

      !enforce essential boundary conditions
      do 250 i=1,ndbcs !number of displacement bcs
        ndof=dbcdof(i) !nodal degree of freedom (1-ndof) of each disp bc
        if (ndof.le.ndim) then
          u(dbcnd(i),ndof)=dbcval(i)
        else
          if (ndof.eq.ndim+1) then
            phi(dbcnd(i))=dbcval(i)
          elseif (ndof.le.ndpn-1) then
            P(3,dbcnd(i),ndof-ndim-1)=dbcval(i)
          else
            theta(dbcnd(i))=dbcval(i)
          endif
        endif
 250  continue

      !enforce multipoint constraints (loop 300)
      do 300 i=1,ncon

        !dof 1
        if (mpcnd(i,2).eq.1) then
          u(mpcnd(i,1),1)=0.d0
          do 310 j=1,mpcnum(i)
            if(mpc(i,j,2).eq.1) u(mpcnd(i,1),1)
     1        =u(mpcnd(i,1),1)+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.2) u(mpcnd(i,1),1)
!     1        =u(mpcnd(i,1),1)+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.3) u(mpcnd(i,1),1)
!     1        =u(mpcnd(i,1),1)+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.4) u(mpcnd(i,1),1)
!     1        =u(mpcnd(i,1),1)+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
!            if(mpc(i,j,2).eq.5) u(mpcnd(i,1),1)
!     1        =u(mpcnd(i,1),1)+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.6) u(mpcnd(i,1),1)
!     1        =u(mpcnd(i,1),1)+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.7) u(mpcnd(i,1),1)
!     1        =u(mpcnd(i,1),1)+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.8) u(mpcnd(i,1),1)
!     1        =u(mpcnd(i,1),1)+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 310      continue
        endif

        if (mpcnd(i,2).eq.2) then
          u(mpcnd(i,1),2)=0.d0
          do 320 j=1,mpcnum(i)
!            if(mpc(i,j,2).eq.1) u(mpcnd(i,1),2)
!     1        =u(mpcnd(i,1),2)+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
            if(mpc(i,j,2).eq.2) u(mpcnd(i,1),2)
     1        =u(mpcnd(i,1),2)+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.3) u(mpcnd(i,1),2)
!     1        =u(mpcnd(i,1),2)+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.4) u(mpcnd(i,1),2)
!     1        =u(mpcnd(i,1),2)+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
!            if(mpc(i,j,2).eq.5) u(mpcnd(i,1),2)
!     1        =u(mpcnd(i,1),2)+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.6) u(mpcnd(i,1),2)
!     1        =u(mpcnd(i,1),2)+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.7) u(mpcnd(i,1),2)
!     1        =u(mpcnd(i,1),2)+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.8) u(mpcnd(i,1),2)
!     1        =u(mpcnd(i,1),2)+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 320      continue
        endif

        if (mpcnd(i,2).eq.3) then
          u(mpcnd(i,1),3)=0.d0
          do 330 j=1,mpcnum(i)
!            if(mpc(i,j,2).eq.1) u(mpcnd(i,1),3)
!     1        =u(mpcnd(i,1),3)+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.2) u(mpcnd(i,1),3)
!     1        =u(mpcnd(i,1),3)+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
            if(mpc(i,j,2).eq.3) u(mpcnd(i,1),3)
     1        =u(mpcnd(i,1),3)+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.4) u(mpcnd(i,1),3)
!     1        =u(mpcnd(i,1),3)+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
!            if(mpc(i,j,2).eq.5) u(mpcnd(i,1),3)
!     1        =u(mpcnd(i,1),3)+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.6) u(mpcnd(i,1),3)
!     1        =u(mpcnd(i,1),3)+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.7) u(mpcnd(i,1),3)
!     1        =u(mpcnd(i,1),3)+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.8) u(mpcnd(i,1),3)
!     1        =u(mpcnd(i,1),3)+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 330      continue
        endif

        if (mpcnd(i,2).eq.4) then
          phi(mpcnd(i,1))=0.d0
          do 340 j=1,mpcnum(i)
!            if(mpc(i,j,2).eq.1) phi(mpcnd(i,1))
!     1        =phi(mpcnd(i,1))+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.2) phi(mpcnd(i,1))
!     1        =phi(mpcnd(i,1))+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.3) phi(mpcnd(i,1))
!     1        =phi(mpcnd(i,1))+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
            if(mpc(i,j,2).eq.4) phi(mpcnd(i,1))
     1        =phi(mpcnd(i,1))+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
!            if(mpc(i,j,2).eq.5) phi(mpcnd(i,1))
!     1        =phi(mpcnd(i,1))+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.6) phi(mpcnd(i,1))
!     1        =phi(mpcnd(i,1))+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.7) phi(mpcnd(i,1))
!     1        =phi(mpcnd(i,1))+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.8) phi(mpcnd(i,1))
!     1        =phi(mpcnd(i,1))+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 340      continue
        endif

        if (mpcnd(i,2).eq.5) then
          P(3,mpcnd(i,1),1)=0.d0
          do 350 j=1,mpcnum(i)
!            if(mpc(i,j,2).eq.1) P(3,mpcnd(i,1),1)
!     1        =P(3,mpcnd(i,1),1)+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.2) P(3,mpcnd(i,1),1)
!     1        =P(3,mpcnd(i,1),1)+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.3) P(3,mpcnd(i,1),1)
!     1        =P(3,mpcnd(i,1),1)+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.4) P(3,mpcnd(i,1),1)
!     1        =P(3,mpcnd(i,1),1)+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
            if(mpc(i,j,2).eq.5) P(3,mpcnd(i,1),1)
     1        =P(3,mpcnd(i,1),1)+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.6) P(3,mpcnd(i,1),1)
!     1        =P(3,mpcnd(i,1),1)+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.7) P(3,mpcnd(i,1),1)
!     1        =P(3,mpcnd(i,1),1)+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.8) P(3,mpcnd(i,1),1)
!     1        =P(3,mpcnd(i,1),1)+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 350      continue
        endif

        if (mpcnd(i,2).eq.6) then
          P(3,mpcnd(i,1),2)=0.d0
          do 360 j=1,mpcnum(i)
!            if(mpc(i,j,2).eq.1) P(3,mpcnd(i,1),2)
!     1        =P(3,mpcnd(i,1),2)+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.2) P(3,mpcnd(i,1),2)
!     1        =P(3,mpcnd(i,1),2)+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.3) P(3,mpcnd(i,1),2)
!     1        =P(3,mpcnd(i,1),2)+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.4) P(3,mpcnd(i,1),2)
!     1        =P(3,mpcnd(i,1),2)+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
!            if(mpc(i,j,2).eq.5) P(3,mpcnd(i,1),2)
!     1        =P(3,mpcnd(i,1),2)+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
            if(mpc(i,j,2).eq.6) P(3,mpcnd(i,1),2)
     1        =P(3,mpcnd(i,1),2)+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.7) P(3,mpcnd(i,1),2)
!     1        =P(3,mpcnd(i,1),2)+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.8) P(3,mpcnd(i,1),2)
!     1        =P(3,mpcnd(i,1),2)+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 360      continue
        endif

        if (mpcnd(i,2).eq.7) then
          P(3,mpcnd(i,1),3)=0.d0
          do 370 j=1,mpcnum(i)
!            if(mpc(i,j,2).eq.1) P(3,mpcnd(i,1),3)
!     1        =P(3,mpcnd(i,1),3)+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.2) P(3,mpcnd(i,1),3)
!     1        =P(3,mpcnd(i,1),3)+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.3) P(3,mpcnd(i,1),3)
!     1        =P(3,mpcnd(i,1),3)+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.4) P(3,mpcnd(i,1),3)
!     1        =P(3,mpcnd(i,1),3)+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
!            if(mpc(i,j,2).eq.5) P(3,mpcnd(i,1),3)
!     1        =P(3,mpcnd(i,1),3)+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.6) P(3,mpcnd(i,1),3)
!     1        =P(3,mpcnd(i,1),3)+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
            if(mpc(i,j,2).eq.7) P(3,mpcnd(i,1),3)
     1        =P(3,mpcnd(i,1),3)+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.8) P(3,mpcnd(i,1),3)
!     1        =P(3,mpcnd(i,1),3)+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 370      continue
        endif

        if (mpcnd(i,2).eq.8) then
          theta(mpcnd(i,1))=0.d0
          do 380 j=1,mpcnum(i)
!            if(mpc(i,j,2).eq.1) theta(mpcnd(i,1))
!     1        =theta(mpcnd(i,1))+conA(i,j)*u(mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.2) theta(mpcnd(i,1))
!     1        =theta(mpcnd(i,1))+conA(i,j)*u(mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.3) theta(mpcnd(i,1))
!     1        =theta(mpcnd(i,1))+conA(i,j)*u(mpc(i,j,1),3)+conB(i,j)
!            if(mpc(i,j,2).eq.4) theta(mpcnd(i,1))
!     1        =theta(mpcnd(i,1))+conA(i,j)*phi(mpc(i,j,1))+conB(i,j)
!            if(mpc(i,j,2).eq.5) theta(mpcnd(i,1))
!     1        =theta(mpcnd(i,1))+conA(i,j)*P(3,mpc(i,j,1),1)+conB(i,j)
!            if(mpc(i,j,2).eq.6) theta(mpcnd(i,1))
!     1        =theta(mpcnd(i,1))+conA(i,j)*P(3,mpc(i,j,1),2)+conB(i,j)
!            if(mpc(i,j,2).eq.7) theta(mpcnd(i,1))
!     1        =theta(mpcnd(i,1))+conA(i,j)*P(3,mpc(i,j,1),3)+conB(i,j)
            if(mpc(i,j,2).eq.8) theta(mpcnd(i,1))
     1        =theta(mpcnd(i,1))+conA(i,j)*theta(mpc(i,j,1))+conB(i,j)
 380      continue
        endif

 300  continue

!      print*,'displacements at nodes'
!      do i=1,nodes
!        print 966, i,u(i,1),u(i,2),u(i,3)
!      enddo
! 966  format(I2,F8.3,F8.3,F8.3)
      return
      end

*************************************************************************

      subroutine gradients(ndim,nnds,nints,maxele,maxnd,neles,ele,B,
     1                     NN,u,strn,phi,E,P,Pint,PintPrev,gradP,
     2                     theta,thint,gradTheta)
      implicit none
      real*8 B(maxele,nnds,nints,ndim),NN(nnds,nints)
      real*8 u(maxnd,ndim),strn(maxele,nints,ndim,ndim)
      real*8 phi(maxnd),E(maxele,nints,ndim)
      real*8 P(3,maxnd,ndim),Pint(maxele,nints,ndim)
      real*8 PintPrev(maxele,nints,ndim),gradP(maxele,nints,ndim,ndim)
      real*8 theta(maxnd),thint(maxele,nints)
      real*8 gradTheta(maxele,nints,ndim)
      integer maxele,maxnd,neles,ele(maxele,8)
      integer i,j,k,m,jj,q,l,ndim,nnds,nints

      do 701 m=1,neles
        do 702 j=1,nints
          ! initialize all derived values to zero
          thint(m,j)=0.d0
          do 703 i=1,ndim
            Pint(m,j,i)=0.d0
            PintPrev(m,j,i)=0.d0
            E(m,j,i)=0.d0
            gradTheta(m,j,i)=0.d0
            do 704 jj=1,ndim
              if(i.eq.1.AND.jj.eq.1) then
                strn(m,j,i,jj)=0.d0
              elseif(i.eq.2.AND.jj.eq.2) then
                strn(m,j,i,jj)=0.d0
              else
                strn(m,j,i,jj)=0.d0
              endif
              gradP(m,j,i,jj)=0.d0
 704        continue
 703      continue

          do 705 i=1,nnds

            k=ele(m,i) !gives node number

            ! format: B(maxele/m,nnds/i,nints/j,ndim/l)
            do l=1,ndim

              Pint(m,j,l)=Pint(m,j,l)+NN(i,j)*P(3,k,l)
              PintPrev(m,j,l)=PintPrev(m,j,l)+NN(i,j)*P(2,k,l)

              E(m,j,l)=E(m,j,l)-B(m,i,j,l)*phi(k)

              gradTheta(m,j,l)=gradTheta(m,j,l)+B(m,i,j,l)*theta(k)
           !if (l.eq.1) print*,m,l,gradTheta(m,j,l),B(m,i,j,l),theta(k)
              do q=1,ndim
                strn(m,j,q,l)=strn(m,j,q,l)+B(m,i,j,l)*u(k,q)
                gradP(m,j,q,l)=gradP(m,j,q,l)+B(m,i,j,l)*P(3,k,q)
              enddo

            enddo

            thint(m,j)=thint(m,j)+NN(i,j)*theta(k)

            if (isnan(Pint(m,j,l))) then
            print*,'Pint',m,j,l,'gives NAN'
            stop
            endif

            if (isnan(strn(m,j,2,2))) then
            print*,'y-strain',m,j,'gives NAN'
            stop
            endif

 705      continue

!          print*,gradTheta(m,j,1),gradTheta(m,j,2),gradTheta(m,j,3)

!            if (abs(thint(m,j)-1.d0).gt.1.d-10) then
!            print*,'HEY',m,j,thint(m,j)
!            stop
!            endif

 702    continue
 701  continue

!      print*
!      print 967,'strains:',strn(1,1,3,3),strn(1,1,1,3),
!     1        strn(1,1,3,1),strn(1,1,2,3)
! 967  format(A10,F8.3,F8.3,F8.3,F8.3)
      return
      end

*************************************************************************

      subroutine stresses(ndim,nints,maxele,neles,strn,E,Pint,gradP,
     1                 thint,gradTheta,EC0cop1,strs,D,psiP,ksi,s,qt,a0,
     2                 beta_P,Psi,det,area,avg_stress,avg_E,
     3                 avg_D,avg_S,avg_T,avgTtot,dt,svm)
      implicit none

      integer maxele,neles,nints,ndim
      real*8 strn(maxele,nints,ndim,ndim),E(maxele,nints,ndim)
      real*8 Pint(maxele,nints,ndim),gradP(maxele,nints,ndim,ndim)
      real*8 thint(maxele,nints),gradTheta(maxele,nints,ndim)

      real*8 strs(maxele,nints,ndim+1,ndim+1),D(maxele,nints,ndim)
      real*8 ksi(maxele,nints,ndim,ndim),psiP(maxele,nints,ndim)
      real*8 s(maxele,nints),qt(maxele,nints,ndim)

      real*8 c11,c12,c44
      real*8 a0,a1c,a11,a12,a111,a112,a123,a1111,a1112,a1122,a1123
      real*8 Q11,Q12,Q44,Q11n
      real*8 cep,Ta,Tb,beta1,ktherm
      real*8 T0,E0,s0,P0,ee0,k0,dt,l0,beta_P
      real*8 Px,Py,Pz,exx,eyy,ezz,exy,eyz,ezx,Ex,Ey,Ez,T
      integer m,l,i,j,k

      real*8 lambda,mu,cpcop,Tacop,Tbcop,betacop,kthermcop
      integer EC0cop1(maxele)

      real*8 svm(maxele,nints)

      real*8 Psi(maxele,nints)
      real*8 avg_stress(ndim,ndim),avg_E(ndim),avg_D(ndim),avg_S,avg_T
      real*8 det(maxele,nints),area,ww(nints),avgTtot

! There are two material, so we must account for the material
! properties of each one

*** electrocaloric material
! CURIE TEMPERATURE
      T0=273.d0+124.99725d0

! VALUE OF P AT FREE ENERGY MIN IN TETRAGONAL STATE
      P0=0.181144132267376282d0

! FROM CHEN APPENDIX
      Q11n=0.10d0
! Spontaneous strain = Q11*spontaneous polarization^2
      ee0=Q11n*P0**2.d0

! SPONTANEOUS ELECTRIC FIELD FROM d_PHI/d_P VS E PLOT
      E0=3.205504976500d5

! sigma_o, a derived quantity used for normalization
      s0=E0*P0/ee0

! characteristic length scale (for normalization)
      l0=sqrt(a0*P0/E0)

! Permittivity of free space (known value)
      k0=-8.854d-12*E0/P0

! TEMPERATURE DEPENDENT TERMS
! from li et al. 2005
      a1c=4.124d5*T0*P0/E0 !constant of a1 term (see above comment)
      a11=(-2.097d8)*P0**3.d0/E0
      a12=(7.974d8)*P0**3.d0/E0
      a111=(1.294d9)*P0**5.d0/E0
      a112=-1.950d9*P0**5.d0/E0
      a123=-2.500d9*P0**5.d0/E0
      a1111=3.863d10*P0**7.d0/E0
      a1112=2.529d10*P0**7.d0/E0
      a1122=1.637d10*P0**7.d0/E0
      a1123=1.367d10*P0**7.d0/E0

! Chen Appendix numbers
      c11=1.78d11*ee0/s0
      c12=0.964d11*ee0/s0
      c44=1.22d11*ee0/s0

      Q11=0.10d0*P0**2.d0/ee0
      Q12=-0.034d0*P0**2.d0/ee0
      Q44=0.029d0*P0**2.d0/ee0

! Heat capacity
      cep=2.5324d6*T0/(E0*P0) !2.5324×106 J/(m3 K)

! Thermal expansion (=0 for now)
      beta1=3.d-6*T0/s0

! Some reference temperature for beta1 and cep
      Ta=1.d0*T0/T0
      Tb=1.d0*T0/T0

! Thermal conductivity
      ktherm=2.61d0*dt*T0/(P0*E0*l0**2)

*** electrode material
       lambda=71.d9*ee0/s0
       mu=33.4d9*ee0/s0
       betacop=18.d-6*T0/s0
       Tacop=1.d0*T0/T0
       Tbcop=1.d0*T0/T0
       cpcop=3.40252d6*T0/(E0*P0)
       kthermcop=401.d0*dt*T0/(P0*E0*l0**2)

! find derived values at integration stations
      do 100 m=1,neles
        if (EC0cop1(m).eq.0) then
        do 200 l=1,nints

          Px=Pint(m,l,1)
          Py=Pint(m,l,2)
          Pz=Pint(m,l,3)
          exx=strn(m,l,1,1)
          eyy=strn(m,l,2,2)
          ezz=strn(m,l,3,3)
          exy=(strn(m,l,1,2)+strn(m,l,2,1))/2.d0
          eyz=(strn(m,l,2,3)+strn(m,l,3,2))/2.d0
          ezx=(strn(m,l,3,1)+strn(m,l,1,3))/2.d0
          Ex=E(m,l,1)
          Ey=E(m,l,2)
          Ez=E(m,l,3)
          T=thint(m,l)

          strs(m,l,1,1)= (T - Ta)*beta1 -
     1     c11*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) -
     2     c12*(2.d0*Px**2*Q12 + Py**2*(Q11 + Q12) +
     3       Pz**2*(Q11 + Q12) - eyy - ezz)

          strs(m,l,2,2)= (T - Ta)*beta1 -
     1     c11*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) -
     2     c12*(2.d0*Py**2*Q12 + Px**2*(Q11 + Q12) +
     3       Pz**2*(Q11 + Q12) - exx - ezz)

          strs(m,l,3,3)= (T - Ta)*beta1 -
     1     c12*(2.d0*Pz**2*Q12 + Px**2*(Q11 + Q12) +
     2       Py**2*(Q11 + Q12) - exx - eyy) -
     3     c11*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)

          strs(m,l,1,2)=c44*(2.d0*exy - Px*Py*Q44)
          strs(m,l,2,3)=c44*(2.d0*eyz - Py*Pz*Q44)
          strs(m,l,3,1)=c44*(2.d0*ezx - Px*Pz*Q44)
          strs(m,l,2,1)=strs(m,l,1,2) !due to symmetry
          strs(m,l,3,2)=strs(m,l,2,3) !due to symmetry
          strs(m,l,1,3)=strs(m,l,3,1) !due to symmetry

          D(m,l,1)=Px+Ex*k0
          D(m,l,2)=Py+Ey*k0
          D(m,l,3)=Pz+Ez*k0

          psiP(m,l,1)=-Ex + 4.d0*Px**3*a11 + 6.d0*Px**5*a111 +
     1  8.d0*Px**7*a1111 + (2.d0*Px*Py**6 + 2.d0*Px*Pz**6 +
     2     6.d0*Px**5*(Py**2 + Pz**2))*a1112 +
     3  (4.d0*Px**3*Py**2 + 4.d0*Px**3*Pz**2 +
     4     2.d0*Px*(Py**4 + Pz**4))*a112 +
     5  4.d0*Px**3*(Py**4 + Pz**4)*a1122 +2.d0*Px**3*Py**2*Pz**2*a1123 +
     6  2.d0*Px*Py**2*Pz**2*(Px**2 + Py**2 + Pz**2)*a1123 +
     7  2.d0*Px*(Py**2 + Pz**2)*a12 + 2.d0*Px*Py**2*Pz**2*a123 +
     8  2.d0*Px*(T - 388/T0)*a1c +(c44*(2.d0*Px*(Py**2 + Pz**2)*Q44**2 -
     9       2.d0*Q44*(Py*(exy + exy) + Pz*(ezx + ezx))))/2.d0 +
     a  c12*(4*Px*Q12*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) +
     b     2.d0*Px*Q11*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) +
     c     2.d0*Px*Q12*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) +
     d     2.d0*Px*Q11*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz) +
     e     2.d0*Px*Q12*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)) +
     f  (c11*(4.d0*Px*Q11*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) +
     g     4.d0*Px*Q12*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy)+
     h     4.d0*Px*Q12*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)))/2.d0

          psiP(m,l,2)=-Ey + 4.d0*Py**3*a11 + 6.d0*Py**5*a111 +
     1    8.d0*Py**7*a1111 + (2.d0*Px**6*Py + 2.d0*Py*Pz**6 +
     2     6.d0*Py**5*(Px**2 + Pz**2))*a1112 +
     3  (4.d0*Px**2*Py**3 + 4.d0*Py**3*Pz**2 +
     4     2.d0*Py*(Px**4 + Pz**4))*a112 +
     5  (4.d0*Px**4*Py**3 + 4.d0*Py**3*Pz**4)*a1122 +
     6  2.d0*Px**2*Py**3*Pz**2*a1123 +
     7  2.d0*Px**2*Py*Pz**2*(Px**2 + Py**2 + Pz**2)*a1123 +
     8   (2.d0*Px**2*Py + 2.d0*Py*Pz**2)*a12 +
     9  2.d0*Px**2*Py*Pz**2*a123 + 2.d0*Py*(T - 388/T0)*a1c +
     a  (c44*(2.d0*Px**2*Py*Q44**2 + 2.d0*Py*Pz**2*Q44**2 -
     b       2.d0*Px*Q44*(exy + exy) - 2.d0*Pz*Q44*(eyz + eyz)))/2.d0 +
     c   c12*(2.d0*Py*Q11*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) +
     d     2.d0*Py*Q12*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) +
     e     4.d0*Py*Q12*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) +
     f     2.d0*Py*Q11*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz) +
     g     2.d0*Py*Q12*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)) +
     h  (c11*(4.d0*Py*Q12*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx)+
     i      4.d0*Py*Q11*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy)+
     j      4.d0*Py*Q12*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)))/2.d0

          psiP(m,l,3)= -Ez + 4.d0*Pz**3*a11 + 6.d0*Pz**5*a111 +
     1  8.d0*Pz**7*a1111 + (2.d0*Px**6*Pz + 2.d0*Py**6*Pz +
     2     6.d0*(Px**2 + Py**2)*Pz**5)*a1112 +
     3  (2.d0*(Px**4 + Py**4)*Pz + 4.d0*Px**2*Pz**3 +
     4     4.d0*Py**2*Pz**3)*a112 +
     5  (4.d0*Px**4*Pz**3 + 4.d0*Py**4*Pz**3)*a1122 +
     6  2.d0*Px**2*Py**2*Pz**3*a1123 +
     7  2.d0*Px**2*Py**2*Pz*(Px**2 + Py**2 + Pz**2)*a1123 +
     8    (2.d0*Px**2*Pz + 2.d0*Py**2*Pz)*a12 +
     9  2.d0*Px**2*Py**2*Pz*a123 + 2.d0*Pz*(T - 388/T0)*a1c +
     a  (c44*(2.d0*Px**2*Pz*Q44**2 + 2.d0*Py**2*Pz*Q44**2 -
     b       2.d0*Px*Q44*(ezx + ezx) - 2.d0*Py*Q44*(eyz + eyz)))/2.d0 +
     c  (c11*(4.d0*Pz*Q12*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) +
     d     4.d0*Pz*Q12*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) +
     e     4.d0*Pz*Q11*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)))/2.d0+
     f  c12*(2.d0*Pz*Q11*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) +
     g     2.d0*Pz*Q12*(Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx) +
     h     2.d0*Pz*Q11*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) +
     i     2.d0*Pz*Q12*(Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) +
     j     4.d0*Pz*Q12*(Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz))

          do i=1,ndim
            qt(m,l,i)=-ktherm*gradTheta(m,l,i)
            do j=1,ndim
              ksi(m,l,i,j)=a0*gradP(m,l,j,i)
            enddo
          enddo

          ! ENTROPY
          s(m,l)= -(Px**2 + Py**2 + Pz**2)*a1c
     1      - beta1*(exx+eyy+ezz) + cep*log(T/Tb)

          if (isnan(s(m,l))) then
            print*,Px,Py,Pz,strn(m,l,2,2),log(T/Tb),T,Tb
            stop
          endif

 200    continue
        else
        print*,'SOMETHING WRONG, IN COPPER STRESS ROUTINE'
        do 300 l=1,nints

          exx=strn(m,l,1,1)
          eyy=strn(m,l,2,2)
          ezz=strn(m,l,3,3)
          exy=(strn(m,l,1,2)+strn(m,l,2,1))/2.d0
          eyz=(strn(m,l,2,3)+strn(m,l,3,2))/2.d0
          ezx=(strn(m,l,3,1)+strn(m,l,1,3))/2.d0
          Ex=E(m,l,1)
          Ey=E(m,l,2)
          Ez=E(m,l,3)
          T=thint(m,l)

          strs(m,l,1,1)= lambda*(exx+eyy+ezz)
     1     +2.d0*mu*exx+ betacop*(T-Tacop)

          strs(m,l,2,2)= lambda*(exx+eyy+ezz)
     1     +2.d0*mu*eyy+ betacop*(T-Tacop)

          strs(m,l,3,3)= lambda*(exx+eyy+ezz)
     1     +2.d0*mu*ezz+ betacop*(T-Tacop)

          strs(m,l,1,2)= 2.d0*mu*exy
          strs(m,l,2,3)= 2.d0*mu*eyz
          strs(m,l,3,1)= 2.d0*mu*ezx
          strs(m,l,2,1)=strs(m,l,1,2) !due to symmetry
          strs(m,l,3,2)=strs(m,l,1,2) !due to symmetry
          strs(m,l,1,3)=strs(m,l,1,2) !due to symmetry

          D(m,l,1)=Ex*k0
          D(m,l,2)=Ey*k0
          D(m,l,3)=Ez*k0

          do i=1,ndim
            psiP(m,l,i)= 0.d0
            qt(m,l,i)=-kthermcop*gradTheta(m,l,i)
            do j=1,ndim
              ksi(m,l,i,j)=0.d0
            enddo
          enddo

          ! ENTROPY
          s(m,l)=
     1      -betacop*(exx+eyy+ezz)
     2      +cpcop*log(T/Tbcop)

          if (isnan(s(m,l))) then
          print*,'ERROR in copper entropy',log(T/Tbcop),T,Tbcop
          stop
          endif

 300    continue
        endif

 100  continue

      !Find von mises stress
      do m=1,neles
        do l=1,nints
          svm(m,l)=1.d0/sqrt(2.d0)*
     1                sqrt((strs(m,l,1,1)-strs(m,l,2,2))**2.d0
     2                    + (strs(m,l,2,2)-strs(m,l,3,3))**2.d0
     3                    + (strs(m,l,3,3)-strs(m,l,1,1))**2.d0
     4                    + 6.d0*(strs(m,l,1,2)**2.d0
     5                        + strs(m,l,2,3)**2.d0
     6                        + strs(m,l,3,1)**2.d0))
        enddo
      enddo

 ! Find Helmholtz Free Energy density (Psi) at int stats.
      do 101 m=1,neles
        if (EC0cop1(i).eq.0) then
        do 102 l=1,nints

          exx=strn(m,l,1,1)
          eyy=strn(m,l,2,2)
          ezz=strn(m,l,3,3)
          exy=(strn(m,l,1,2)+strn(m,l,2,1))/2.d0
          eyz=(strn(m,l,2,3)+strn(m,l,3,2))/2.d0
          ezx=(strn(m,l,3,1)+strn(m,l,1,3))/2.d0
          Ex=E(m,l,1)
          Ey=E(m,l,2)
          Ez=E(m,l,3)
          T=thint(m,l)

           Psi(i,j)= (Px**4 + Py**4 + Pz**4)*a11 +
     1      (Px**6 + Py**6 + Pz**6)*a111 +
     2      (Px**8 + Py**8 + Pz**8)*a1111 +
     3      ((Px**2 + Py**2)*Pz**6 + Py**6*(Px**2 + Pz**2) +
     4         Px**6*(Py**2 + Pz**2))*a1112 +
     5      ((Px**4 + Py**4)*Pz**2 + Py**2*(Px**4 + Pz**4) +
     6         Px**2*(Py**4 + Pz**4))*a112 +
     7      (Py**4*Pz**4 + Px**4*(Py**4 + Pz**4))*a1122 +
     8      Px**2*Py**2*Pz**2*(Px**2 + Py**2 + Pz**2)*a1123 +
     9      (Py**2*Pz**2 + Px**2*(Py**2 + Pz**2))*a12 +
     a      Px**2*Py**2*Pz**2*a123 +
     b      (Px**2 + Py**2 + Pz**2)*(T - 388/T0)*a1c +
     c      (c44*(Py**2*Pz**2*Q44**2 + Px**2*(Py**2 + Pz**2)*Q44**2 +
     d           exy**2 + ezx**2 + 2.d0*exy*exy + exy**2 + eyz**2 +
     e           2.d0*ezx*ezx + ezx**2 -
     f           2.d0*Px*Q44*(Py*(exy + exy) + Pz*(ezx + ezx)) +
     g           2.d0*eyz*eyz + eyz**2 -2.d0*Py*Pz*Q44*(eyz + eyz)))/2.+
     h      c12*((Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx)*
     i          (Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy) +
     j         (Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx)*
     k          (Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz) +
     l         (Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy)*
     m          (Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)) +
     n      (c11*((Px**2*Q11 + Py**2*Q12 + Pz**2*Q12 - exx)**2 +
     o           (Py**2*Q11 + Px**2*Q12 + Pz**2*Q12 - eyy)**2 +
     p           (Pz**2*Q11 + Px**2*Q12 + Py**2*Q12 - ezz)**2))/2. +
     q      (T - Ta)*beta1*(exx + eyy + ezz) +
     r      ((Ex**2 + Ey**2 + Ez**2)*k0)/2. -
     s      cep*(-T + Tb + T*Log(T/Tb))

 102    continue
        endif
 101  continue

      if (nints.eq.8) then
        do i=1,nints
          ww(i)=1.d0
        enddo
      endif

      ! find the average stress,D,E,S,T in the EC material
      avg_S=0.d0
      avg_T=0.d0
      area=0.d0
      avgTtot=0.d0
      do 109 k=1,ndim
        avg_D(k)=0.d0
        avg_E(k)=0.d0
        do 110 l=1,ndim
          avg_stress(k,l)=0.d0
 110    continue
 109  continue

      ! calculate avg stress
      do 103 i=1,neles
        if (EC0cop1(i).eq.0) then
        do 104 j=1,nints
          avg_S=avg_S+s(i,j)*ww(j)*det(i,j)
          avg_T=avg_T+thint(i,j)*ww(j)*det(i,j)
          area=area+ww(j)*det(i,j)
          do 105 k=1,ndim
            do 106 l=1,ndim
            avg_stress(k,l)=avg_stress(k,l)+strs(i,j,k,l)*ww(j)*det(i,j)
 106        continue
            avg_D(k)=avg_D(k)+D(i,j,k)*ww(j)*det(i,j)
            avg_E(k)=avg_E(k)+E(i,j,k)*ww(j)*det(i,j)
 105      continue
 104    continue
        endif
 103  continue

      do 129 i=1,neles
        do 129 j=1,nints
          avgTtot=avgTtot+thint(i,j)*ww(j)*det(i,j)
 129  continue
      avgTtot=avgTtot/area

      do 107 i=1,ndim
        do 108 j=1,ndim
          avg_stress(i,j)=avg_stress(i,j)/area
 108    continue
        avg_D(i)=avg_D(i)/area
        avg_E(i)=avg_E(i)/area
 107  continue
      avg_S=avg_S/area
      avg_T=avg_T/area

!      print*
!      print*,'avg values'
!      print*,avg_S,avg_T
!      print*,avg_E(1),avg_E(2),avg_E(3)
!      print*,avg_stress(1,1),avg_stress(2,2),avg_stress(3,3)
!      print*,avg_stress(1,2),avg_stress(2,3),avg_stress(3,1)

!      print*,Psi(1,1)

      return
      end

************************************************************************

       subroutine intforce(ndim,nints,nnds,ndpn,maxnd,nodes,ny,maxele,
     1                    neles,ele,gcon,det,B,NN,strs,D,psiP,ksi,s,qt,
     2                    finte,dt,beta_P,Pint,PintPrev,Psi,
     3                    PsiTot,PsiPrint,Sprint,area,h,node,control,
     4                     EC0cop1)
      implicit none

      real*8 node(maxnd,ndim),det(maxele,nints),NN(nnds,nints)
      real*8 B(maxele,nnds,nints,ndim),ww(nints)

      real*8 strs(maxele,nints,ndim+1,ndim+1),D(maxele,nints,ndim)
      real*8 psiP(maxele,nints,ndim),ksi(maxele,nints,ndim,ndim)
      real*8 s(maxele,nints),qt(maxele,nints,ndim)
      real*8 Pint(maxele,nints,ndim),PintPrev(maxele,nints,ndim)
      real*8 dt,beta_P,T0,Tcn,finte(maxnd*ndpn)
      real*8 fu(nnds,ndim),fphi(nnds),fP(nnds,ndim),ftheta(nnds)

      integer maxnd,maxele,nodes,neles
      integer ndim,nints,nnds,ndpn
      integer ele(maxele,8),gcon(maxnd,ndpn)
      integer ny,gdof
      integer EC0cop1(maxele)

      real*8 Psi(maxele,nints),PsiTot
      !real*8 PsiFerro,PsiPara
      !real*8 af,ap,
      real*8 h,area
      real*8 STot

      integer i,j,k,l,ii,jj,n
      integer SFile,SPrint,PsiTotFile,PsiPrint,control

      PsiTotFile=11
      SFile=12

! CURIE TEMPERATURE
      T0=273.d0+124.99725d0
      Tcn=T0/T0

      if (nints.eq.8) then
        do i=1,nints
          ww(i)=1.d0
        enddo
      endif

      !initialize internal forces to zero
      do 100 i=1,nodes*ndpn
        finte(i)=0.d0
!        print*,i,finte(i),'challah'
 100  continue

      do 200 i=1,neles

        !initialize internal forces at for given element to zero
        do 210 j=1,nnds
          fphi(j)=0.d0
          ftheta(j)=0.d0
          do 211 k=1,ndim
            fu(j,k)=0.d0
            fP(j,k)=0.d0
 211      continue
 210    continue

        do 220 k=1,nints
          do 230 l=1,nnds

           ftheta(l)=ftheta(l)-ww(k)*det(i,k)*s(i,k)*NN(l,k)/dt

           ! B structure: B(maxele,nnds,nints,ndim)
           do 221 j=1,ndim
             fphi(l)=fphi(l)+det(i,k)*ww(k)*D(i,k,j)*B(i,l,k,j)/dt
!             if (EC0cop1(i).eq.0) then
             fP(l,j)=fP(l,j)+det(i,k)*ww(k)*NN(l,k)*
     1                       (psiP(i,k,j)/dt+beta_P*Pint(i,k,j)/(dt**2))
             ftheta(l)=ftheta(l)+(beta_P*Pint(i,k,j)*(Pint(i,k,j) -
     1                        2.d0*PintPrev(i,k,j))*NN(l,k)/(Tcn*dt**2)+
     2                        qt(i,k,j)*B(i,l,k,j)/Tcn)*det(i,k)*ww(k)
!             else
!               ftheta(l)=ftheta(l)+(qt(i,k,j)*B(i,l,k,j)/Tcn)*
!     1                                               det(i,k)*ww(k)
!             endif

 221       continue

           do 222 ii=1,ndim
             do 223 jj=1,ndim
               fu(l,ii)=fu(l,ii)+det(i,k)*ww(k)*B(i,l,k,jj)
     1                 *strs(i,k,jj,ii)/dt
!               if (EC0cop1(i).eq.0) then
                 fP(l,ii)=fP(l,ii)+det(i,k)*ww(k)*ksi(i,k,jj,ii)*
     1                 B(i,l,k,jj)/dt
!               endif
 223         continue
 222       continue

 230      continue
 220    continue

        do 300 n=1,nnds !local nodes
            gdof=gcon(ele(i,n),ndim+1) !For phi for element i, locnode n
            finte(gdof)=finte(gdof)+fphi(n)

            gdof=gcon(ele(i,n),ndpn)
            finte(gdof)=finte(gdof)+ftheta(n)

            do 310 l=1,ndim
              gdof=gcon(ele(i,n),l)
              finte(gdof)=finte(gdof)+fu(n,l)

              gdof=gcon(ele(i,n),ndim+1+l)
              finte(gdof)=finte(gdof)+fP(n,l)

 310        continue
 300    continue
 200  continue

!      do 123 i=1,nodes*ndpn
!        print*,i,finte(i),'challah'
! 123  continue
      PsiTot=0.d0
      STot=0.d0
      do 101 i=1,neles
        do 102 j=1,nints
          Stot=Stot+s(i,j)*ww(j)*det(i,j)/area
          PsiTot=PsiTot+Psi(i,j)*ww(j)*det(i,j)/area
 102    continue
 101  continue

!      if (PsiPrint.eq.1) then
!        open(unit=PsiTotFile,file='PsiTot',access='append')
!        write(PsiTotFile,*) PsiTot,PsiFerro/af,PsiPara/ap
!        close(unit=PsiTotFile)
!      endif
!      print*,'Psi=',PsiTot

!      do 321 i=6,10
!         print*,finte(gcon(i,5))
! 321  continue
      return
      end

*************************************************************************
******* STIFFNESS ROUTINES
*************************************************************************

      subroutine stf(ndim,nnds,nints,ndpn,B,NN,det,
     1               maxele,maxnd,maxdbc,maxcon,ndofs,neles,
     2               ele,gcon,strn,mpc,conA,ndbcs,ncon,mpcnum,
     3               dt,beta_P,Pint,PintPrev,thint,a0,control,
     4               nodes,PetA,eleproc,Petx,Petb,rank,EC0cop1)

#include "petsc/finclude/petsc.h" 
      USE petsc     
      implicit none

      real*8 B(maxele,nnds,nints,ndim),det(maxele,nints)
      real*8 kele(nnds*ndpn*(nnds*ndpn+1)/2)
      real*8 keTemp1(nnds*ndpn,nnds*ndpn)
      real*8 keTemp2(nnds*ndpn*2,nnds*ndpn*2)
      real*8 bele(nnds,nints,ndim),detele(nints)
      real*8 strn(maxele,nints,ndim,ndim)
      real*8 conA(maxdbc,maxcon),dt,beta_P,NN(nnds,nints),kcon
      real*8 thint(maxele,nints),a0,Pint(maxele,nints,ndim)
      real*8 PintPrev(maxele,nints,ndim)
      integer maxele,maxnd,maxdbc,maxcon,ndbcs,ncon
!      integer removerows
      integer neles,ele(maxele,8),gcon(maxnd,ndpn)
      integer i,j,k,l,n,nloci,nlocj,dofi,dofj,ndofs,coni,conj
      integer mpc(maxdbc,maxcon,2),nci,ncj,ii,jj
      integer mpcnum(maxdbc),ndim,nnds,nints,ndpn,control,nodes
      integer EC0cop1(maxele),EC0cop1ele

      logical eleproc(maxele),flgi,flgj
      integer tmpi,tmpj,checkij,lploci

      Mat PetA
      Vec Petb,Petx
!      Vec DD
      PetscErrorCode ierr
      PetscInt pkn
      PetscInt ploci(nnds*ndpn*2)
      PetscScalar zero
      PetscMPIInt rank
!      PetscViewer Aview

      pkn=nnds*ndpn*2
      zero=0.d0

      !-------------initialize matrix-----------------
      CALL MatZeroEntries(PetA,ierr)

      do 100 l=1,neles

        lploci=nnds*ndpn !start dim of  local stiffness
        !only do computations if we are in the correct processor for the element
        if (eleproc(l)) then
          !the elements are unique to the processor, however, multiple
          !processors could be adding to the same kij in the stiffness
          !because kij loc is based on the two nodes i and j
          !if (mod(l,1000).eq.0) print*,'element',l
        do 110 i=1,nnds
          do 110 k=1,nints
            do 110 j=1,ndim
              bele(i,k,j)=B(l,i,k,j)
 110    continue

        do 120 i=1,nints
 120      detele(i)=det(l,i)

        EC0cop1ele=EC0cop1(l)
        call stfele(ndim,nnds,nints,ndpn,dt,beta_P,maxele,l,kele,
     1           bele,NN,Pint,PintPrev,detele,strn,thint,a0,control,
     2           EC0cop1ele)

         !kele is a linear storage vector for a 24x24 symmetric matrix
         !only stores diagonal and above
         !n gives the storage location in rank 1
         !i and j are locations for rank 2 (note that i>=j)

        !keTemp1 gives the 24x24 symmetric matrix in rank 2 form
        do 111 i=1,nnds*ndpn
          do 111 j=1,nnds*ndpn
            if (i.ge.j) then
              n=i*(i-1)/2+j
            else
              n=j*(j-1)/2+i
            endif
            keTemp1(i,j)=kele(n)
 111    continue

        !negative integers in ploci ignored by matsetvalue
        !keTemp2 gives the matrix w/ mcp calcs included
        do 114 i=1,nnds*ndpn*2
          ploci(i)=-1
          do 214 j=1,nnds*ndpn*2
            keTemp2(i,j)=0.d0
 214      continue
 114    continue

         do 130 i=1,nnds*ndpn
          do 140 j=1,i
            !n=i*(i-1)/2+j !this is for linear storage kele
            nloci=(i-1)/ndpn+1 !local node for i
            dofi=mod(i-1,ndpn)+1 !local dof for i
            nlocj=(j-1)/ndpn+1 !local node for j
            dofj=mod(j-1,ndpn)+1 !local dof for j
            dofi=gcon(ele(l,nloci),dofi) !global dof for idof
            dofj=gcon(ele(l,nlocj),dofj) !global dof for jdof
            coni=dofi-ndofs-ndbcs
            conj=dofj-ndofs-ndbcs
            if(coni.gt.0) then
              nci=mpcnum(coni)
              do 150 ii=1,nci
                dofi=gcon(mpc(coni,ii,1),mpc(coni,ii,2))
                if (dofi.gt.ndofs) goto 150
                if(conj.gt.0) then
                  ncj=mpcnum(conj)
                  do 155 jj=1,ncj
                    dofj=gcon(mpc(conj,jj,1),mpc(conj,jj,2))
                    if (dofj.gt.ndofs) goto 155
                    if (i.eq.j.and.dofi.gt.dofj) goto 155
                    if((coni.ne.conj).and.(dofi.eq.dofj)) then
                      kcon=2.d0*keTemp1(i,j)
                    else
                      kcon=keTemp1(i,j)
                    endif
                    kcon=conA(coni,ii)*conA(conj,jj)*kcon

!                    if (dofj.eq.14984+1.and.dofi.eq.305+1) print*,
!     1                       l,i,j,ploci(i),ploci(j),kcon

                    flgi=.false.
                    flgj=.false.
                    do 11 checkij=1,nnds*ndpn*2
                     if (dofi-1.eq.ploci(checkij)) then
                        flgi=.true.
                        tmpi=checkij
                        goto 12
                      endif
  11                continue
  12                continue
                    if (.NOT.flgi) then
                      if (ploci(i).eq.-1) then
                        ploci(i)=dofi-1
                        tmpi=i
                      else
                        lploci=lploci+1
                        ploci(lploci)=dofi-1
                        tmpi=lploci
                      endif
                    endif

                    do 13 checkij=1,nnds*ndpn*2
                     if (dofj-1.eq.ploci(checkij)) then
                        flgj=.true.
                        tmpj=checkij
                        goto 14
                      endif
  13                continue
  14                continue
                    if (.NOT.flgj) then
                      if (ploci(j).eq.-1) then
                        tmpj=j
                        ploci(j)=dofj-1
                      else
                        lploci=lploci+1
                        ploci(lploci)=dofj-1
                        tmpj=lploci
                      endif
                    endif

                    keTemp2(tmpi,tmpj)=keTemp2(tmpi,tmpj)+kcon
!                    if (dofj.eq.14984+1.and.dofi.eq.305+1) print*,
!     1            l,ploci(tmpi),ploci(tmpj),tmpi,tmpj,keTemp2(tmpi,tmpj)
!              if (dofi.eq.14984+1.and.dofj.eq.305+1) print*,
!     1            l,ploci(tmpi),ploci(tmpj),tmpi,tmpj,keTemp2(tmpi,tmpj)

 155              continue
                else
                  if (dofj.gt.ndofs) goto 150
                  if(dofi.eq.dofj) then
                    kcon=2.d0*keTemp1(i,j)
                  else
                    kcon=keTemp1(i,j)
                  endif
                  kcon=conA(coni,ii)*kcon

                  flgi=.false.
                  flgj=.false.
                  do 15 checkij=1,nnds*ndpn*2
                    if (dofi-1.eq.ploci(checkij)) then
                      flgi=.true.
                      tmpi=checkij
                      goto 16
                    endif
  15              continue
                  if (.NOT.flgi) then
                    if (ploci(i).eq.-1) then
                      ploci(i)=dofi-1
                      tmpi=i
                    else
                      lploci=lploci+1
                      ploci(lploci)=dofi-1
                      tmpi=lploci
                    endif
                  endif
  16              continue

                  do 17 checkij=1,nnds*ndpn*2
                    if (dofj-1.eq.ploci(checkij)) then
                      flgj=.true.
                      tmpj=checkij
                      goto 18
                    endif
  17              continue
                  if (.NOT.flgj) then
                    if (ploci(j).eq.-1) then
                      tmpj=j
                      ploci(j)=dofj-1
                    else
                      lploci=lploci+1
                      ploci(lploci)=dofj-1
                      tmpj=lploci
                    endif
                  endif
  18              continue

                  keTemp2(tmpi,tmpj)=keTemp2(tmpi,tmpj)+kcon
!                  if (dofj.eq.14984+1.and.dofi.eq.305+1) print*,
!     1            l,ploci(tmpi),ploci(tmpj),tmpi,tmpj,keTemp2(tmpi,tmpj)
!              if (dofi.eq.14984+1.and.dofj.eq.305+1) print*,
!     1            l,ploci(tmpi),ploci(tmpj),tmpi,tmpj,keTemp2(tmpi,tmpj)

                endif
 150          continue
            elseif(conj.gt.0) then
              if (dofi.gt.ndofs) goto 140
              ncj=mpcnum(conj)
              do 160 jj=1,ncj
                dofj=gcon(mpc(conj,jj,1),mpc(conj,jj,2))
                if (dofj.gt.ndofs) goto 160
                if(dofi.eq.dofj) then
                  kcon=2.d0*keTemp1(i,j)
                else
                  kcon=keTemp1(i,j)
                endif
                kcon=conA(conj,jj)*kcon

                flgi=.false.
                flgj=.false.
                do 19 checkij=1,nnds*ndpn*2
                 if (dofi-1.eq.ploci(checkij)) then
                    flgi=.true.
                    tmpi=checkij
                    goto 20
                  endif
  19            continue
  20            continue
                if (.NOT.flgi) then
                  if (ploci(i).eq.-1) then
                    ploci(i)=dofi-1
                    tmpi=i
                  else
                    lploci=lploci+1
                    ploci(lploci)=dofi-1
                    tmpi=lploci
                  endif
                endif

                do 21 checkij=1,nnds*ndpn*2
                 if (dofj-1.eq.ploci(checkij)) then
                    flgj=.true.
                    tmpj=checkij
                    goto 22
                  endif
  21            continue
  22            continue
                if (.NOT.flgj) then
                  if (ploci(j).eq.-1) then
                    tmpj=j
                    ploci(j)=dofj-1
                  else
                    lploci=lploci+1
                    ploci(lploci)=dofj-1
                    tmpj=lploci
                  endif
                endif
                keTemp2(tmpi,tmpj)=keTemp2(tmpi,tmpj)+kcon
!                if (dofj.eq.14984+1.and.dofi.eq.305+1) print*,
!     1            l,ploci(tmpi),ploci(tmpj),tmpi,tmpj,keTemp2(tmpi,tmpj)
!              if (dofi.eq.14984+1.and.dofj.eq.305+1) print*,
!     1            l,ploci(tmpi),ploci(tmpj),tmpi,tmpj,keTemp2(tmpi,tmpj)


 160          continue
            else
              if (dofi.gt.ndofs) goto 140
              if (dofj.gt.ndofs) goto 140
              kcon=keTemp1(i,j)

              flgi=.false.
              flgj=.false.
              do 23 checkij=1,nnds*ndpn*2
                if (dofi-1.eq.ploci(checkij)) then
                  flgi=.true.
                  tmpi=checkij
                  goto 24
                endif
  23          continue
  24          continue
              if (.NOT.flgi) then
                if (ploci(i).eq.-1) then
                  ploci(i)=dofi-1
                  tmpi=i
                else
                  lploci=lploci+1
                  ploci(lploci)=dofi-1
                  tmpi=lploci
                endif
              endif

              do 25 checkij=1,nnds*ndpn*2
                if (dofj-1.eq.ploci(checkij)) then
                  flgj=.true.
                  tmpj=checkij
                  goto 26
                endif
  25          continue
  26          continue
              if (.NOT.flgj) then
                if (ploci(j).eq.-1) then
                  tmpj=j
                  ploci(j)=dofj-1
                else
                  lploci=lploci+1
                  ploci(lploci)=dofj-1
                  tmpj=lploci
                endif
              endif
              keTemp2(tmpi,tmpj)=keTemp2(tmpi,tmpj)+kcon
            endif
 140      continue
 130    continue

        !keTemp2 needs the symmetric side
        do 113 i=1,nnds*ndpn*2
          do 113 j=i+1,nnds*ndpn*2
            if (abs(keTemp2(i,j)).gt.1.d-18.or.
     1                  abs(keTemp2(j,i)).gt.1.d-18) then
               keTemp2(i,j)=keTemp2(i,j)+keTemp2(j,i)
               keTemp2(j,i)=keTemp2(i,j)
            endif
            if (isnan(keTemp2(i,j))) print*,"ISNAN IN KETEMP2(",i,j,")"
 113    continue

        !------------- ASSEMBLE IN PETSC -------------------------------
        CALL MatSetValues(PetA,pkn,ploci,pkn,ploci,keTemp2,
     1                          ADD_VALUES,ierr)
        !---------------------------------------------------------------
!      else
!        print*,"HHHHHHHEEEEEEEEYYYYYYYYY, you done messed up"
!        stop
      endif

 100  continue
      !end of element loop
!      if (rank.eq.0) print*,'out of element loop'
      CALL MatSetOption(PetA,MAT_SYMMETRIC,PETSC_FALSE,ierr)
      CALL MatAssemblyBegin(PetA,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(PetA,MAT_FINAL_ASSEMBLY,ierr)
!      call MatGetDiagonal(PetA,DD,ierr)
      !write diagonal before solving
!            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"diag.m",
!     &                  Aview,ierr)
!            call PetscViewerPushFormat(Aview,
!     1                     PETSC_VIEWER_ASCII_MATLAB,ierr)
!            call VecView(DD,Aview,ierr)
!            call PetscViewerDestroy(Aview)
            !end write
!      removerows=ndbcs+ncon
!      do 115 i=1,removerows
!         dbcrow(i)=ndofs-1+i
!         if (dbcrow(i).gt.ndpn*nodes) print*,"SIZE ERROR IN DBCROW"
! 115  continue
!      CALL VecZeroEntries(Petx,ierr)
!      CALL VecZeroEntries(Petb,ierr)
!      CALL MatZeroRows(PetA,removerows,dbcrow,1.d0,Petx,Petb,ierr)
      CALL MatSetOption(PetA,MAT_NEW_NONZERO_LOCATIONS,
     1                  PETSC_FALSE,ierr)
      !write a matrix
!      call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Aa.m",
!     &                  Aview,ierr)
!      call PetscViewerPushFormat(Aview, PETSC_VIEWER_ASCII_MATLAB,ierr)
!      call MatView(PetA,Aview,ierr)
!      call PetscViewerDestroy(Aview)

      return
      end
******************************************************************************

      subroutine stfele(ndim,nnds,nints,ndpn,dt,beta_P,maxele,
     1                  iele,kele,bele,nele,Pint,PintPrev,detele,strn,
     2                  thint,a0,control,EC0cop1ele)
      implicit none
      real*8 kele(nnds*ndpn*(nnds*ndpn+1)/2)
      real*8 detele(nints),bele(nnds,nints,ndim),nele(nnds,nints)
      real*8 Pint(maxele,nints,ndim),PintPrev(maxele,nints,ndim)
      real*8 strn(maxele,nints,ndim,ndim)
      real*8 thint(maxele,nints)
      real*8 a0,dt,beta_P,ww(nints),delta(ndim,ndim)
      real*8 deijdekl(ndim,ndim,ndim,ndim),deijdPk(ndim,ndim,ndim)
      real*8 deijdT(ndim,ndim),dEidEj(ndim,ndim),dEidPj(ndim,ndim)
      real*8 dPidPk(ndim,ndim),dPidT(ndim)
      real*8 dxiijdxikl(ndim,ndim,ndim,ndim),dTdT
      real*8 ktherm,kthermcop,T0,Tcn
      integer maxele,iele,ia,i,j,k,l,nodeP,nodeQ,kdofP,kdofQ
      integer n,ndim,nnds,nints,ndpn,control
      integer EC0cop1ele

      if (nints.eq.8) then
        do i=1,nints
          ww(i)=1.d0
        enddo
      endif

      do 999 i=1,ndim
        do 998 j=1,ndim
          delta(i,j)=0.d0
          if (i.eq.j) delta(i,j)=1.d0
 998    continue
 999  continue

      !stiffness terms on or above diagonal
      do 100 i=1,nnds*ndpn*(nnds*ndpn+1)/2
 100    kele(i)=0.d0

      do 200 ia=1,nints

        ! Second derivatives of h
        if (EC0cop1ele.eq.0) then
          call materialEC(strn(iele,ia,1,1),strn(iele,ia,2,2),
     1        strn(iele,ia,3,3),(strn(iele,ia,1,2)+strn(iele,ia,2,1))/2,
     2        (strn(iele,ia,2,3)+strn(iele,ia,3,2))/2,
     3        (strn(iele,ia,3,1)+strn(iele,ia,1,3))/2,
     4        Pint(iele,ia,1),Pint(iele,ia,2),Pint(iele,ia,3),
     5        thint(iele,ia),a0,ndim,deijdekl,deijdPk,deijdT,dEidEj,
     6        dEidPj,dPidPk,dPidT,dxiijdxikl,dTdT,ktherm,T0,beta_P,dt)
        else
          call materialcop(thint(iele,ia),ndim,deijdekl,deijdPk,
     1               deijdT,dEidEj,dEidPj,dPidPk,dPidT,dxiijdxikl,
     2               dTdT,kthermcop,T0,a0,dt)
        endif
        Tcn=T0/T0

        ! del u(k)
        do 300 nodeP=1,nnds
          do 400 i=1,ndim
            kdofP=(nodeP-1)*ndpn+i !u
            do 500 nodeQ=1,nnds
              do 600 k=1,ndim
                kdofQ=(nodeQ-1)*ndpn+k !u
                if (kdofQ.gt.kdofP) goto 500
                n=kdofP*(kdofP-1)/2+kdofQ
                do 700 j=1,ndim
                  do 800 l=1,ndim
                    kele(n)=kele(n)+deijdekl(i,j,k,l)*
     1                             bele(nodeP,ia,j)*bele(nodeQ,ia,l)*
     2                             detele(ia)*ww(ia)/dt
 800                continue
 700            continue
 600          continue
 500        continue
 400      continue
 300    continue

        do 306 nodeP=1,nnds
          do 406 i=1,ndim
            kdofP=(nodeP-1)*ndpn+i !u
            do 506 nodeQ=1,nnds
              do 606 k=1,ndim
                kdofQ=(nodeQ-1)*ndpn+ndim+1+k !P
                if (kdofQ.gt.kdofP) goto 506
                n=kdofP*(kdofP-1)/2+kdofQ
                do 706 j=1,ndim
                  kele(n)=kele(n)+deijdPk(i,j,k)*bele(nodeP,ia,j)*
     1                             nele(nodeQ,ia)*detele(ia)*ww(ia)/dt
 706            continue
 606          continue
 506        continue
 406      continue
 306    continue

        do 307 nodeP=1,nnds
          do 407 i=1,ndim
            kdofP=(nodeP-1)*ndpn+i !u
            do 507 nodeQ=1,nnds
              kdofQ=nodeQ*ndpn !T
              if (kdofQ.gt.kdofP) goto 507
              n=kdofP*(kdofP-1)/2+kdofQ
              do 607 j=1,ndim
                kele(n)=kele(n)+deijdT(i,j)*bele(nodeP,ia,j)*
     1                             nele(nodeQ,ia)*detele(ia)*ww(ia)/dt
 607          continue
 507        continue
 407      continue
 307    continue

        !del phi
        do 301 nodeP=1,nnds
          kdofP=(nodeP-1)*ndpn+ndim+1
          do 401 nodeQ=1,nnds
            kdofQ=(nodeQ-1)*ndpn+ndim+1
            if (kdofQ.gt.kdofP) goto 401
            n=kdofP*(kdofP-1)/2+kdofQ
            do 501 i=1,ndim
              do 601 j=1,ndim
                kele(n)=kele(n)+dEidEj(i,j)*bele(nodeP,ia,i)*
     1                           bele(nodeQ,ia,j)*detele(ia)*ww(ia)/dt
 601          continue
 501        continue
 401      continue
 301    continue

        do 302 nodeP=1,nnds
          kdofP=(nodeP-1)*ndpn+ndim+1
          do 402 nodeQ=1,nnds
            do 502 j=1,ndim
              kdofQ=(nodeQ-1)*ndpn+ndim+1+j
              if (kdofQ.gt.kdofP) goto 402
              n=kdofP*(kdofP-1)/2+kdofQ
              do 602 i=1,ndim
                kele(n)=kele(n)-dEidPj(i,j)*bele(nodeP,ia,i)*
     1                             nele(nodeQ,ia)*detele(ia)*ww(ia)/dt
 602          continue
 502        continue
 402      continue
 302    continue

        !del P(k)
        do 308 nodeP=1,nnds
          do 408 i=1,ndim
            kdofP=(nodeP-1)*ndpn+ndim+1+i !P
            do 508 nodeQ=1,nnds
              do 608 k=1,ndim
                kdofQ=(nodeQ-1)*ndpn+k !u
                if (kdofQ.gt.kdofP) goto 608
                n=kdofP*(kdofP-1)/2+kdofQ
                do 708 l=1,ndim
                  kele(n)=kele(n)+deijdPk(k,l,i)*bele(nodeQ,ia,l)*
     1                             nele(nodeP,ia)*detele(ia)*ww(ia)/dt
 708            continue
 608          continue
 508        continue
 408      continue
 308    continue

        do 309 nodeP=1,nnds
          do 409 i=1,ndim
            kdofP=(nodeP-1)*ndpn+ndim+1+i !P
            do 509 nodeQ=1,nnds
              kdofQ=(nodeQ-1)*ndpn+ndim+1 !phi
              if (kdofQ.gt.kdofP) goto 509
              n=kdofP*(kdofP-1)/2+kdofQ
              do 609 k=1,ndim
                kele(n)=kele(n)-dEidPj(k,i)*bele(nodeQ,ia,k)*
     1                             nele(nodeP,ia)*detele(ia)*ww(ia)/dt
 609          continue
 509        continue
 409      continue
 309    continue

        do 303 nodeP=1,nnds
          do 403 i=1,ndim
          kdofP=(nodeP-1)*ndpn+ndim+1+i !P
            do 503 nodeQ=1,nnds
              do 603 k=1,ndim
                kdofQ=(nodeQ-1)*ndpn+ndim+1+k !P
                if (kdofQ.gt.kdofP) goto 603
                n=kdofP*(kdofP-1)/2+kdofQ
                kele(n)=kele(n)+(dPidPk(i,k)/dt+beta_P
     1                        *delta(i,k)/(dt**2))*nele(nodeP,ia)
     2                        *nele(nodeQ,ia)*detele(ia)*ww(ia)
                do 703 j=1,ndim
                  do 803 l=1,ndim
                   kele(n)=kele(n)+dxiijdxikl(i,j,k,l)*bele(nodeP,ia,j)*
     1                             bele(nodeQ,ia,l)*detele(ia)*ww(ia)/dt
 803              continue
 703            continue
 603          continue
 503        continue
 403      continue
 303    continue

        do 304 nodeP=1,nnds
          do 404 i=1,ndim
          kdofP=(nodeP-1)*ndpn+ndim+1+i !P
            do 504 nodeQ=1,nnds
              kdofQ=nodeQ*ndpn !T
              if (kdofQ.gt.kdofP) goto 504
              n=kdofP*(kdofP-1)/2+kdofQ
              kele(n)=kele(n)+dPidT(i)*nele(nodeP,ia)*nele(nodeQ,ia)*
     1                           detele(ia)*ww(ia)/dt
 504        continue
 404      continue
 304    continue

        !del T
        do 310 nodeP=1,nnds
          ! temp dof
          kdofP=nodeP*ndpn
          do 410 nodeQ=1,nnds
            do 510 j=1,ndim
              ! Polarization dof
              kdofQ=(nodeQ-1)*ndpn+ndim+1+j
              if (kdofQ.gt.kdofP) goto 510
              n=kdofP*(kdofP-1)/2+kdofQ
              kele(n)=kele(n)+(dPidT(j)/dt)*
     1                 nele(nodeP,ia)*nele(nodeQ,ia)*detele(ia)*ww(ia)
 510        continue
 410      continue
 310    continue

        do 305 nodeP=1,nnds
          ! temp dof
          kdofP=nodeP*ndpn
          do 405 nodeQ=1,nnds
            ! temp dof
            kdofQ=nodeQ*ndpn
            if (kdofQ.gt.kdofP) goto 405
            n=kdofP*(kdofP-1)/2+kdofQ
            kele(n)=kele(n)+dTdT*nele(nodeP,ia)*nele(nodeQ,ia)*
     1                           detele(ia)*ww(ia)/dt
            do 505 i=1,ndim
              if (EC0cop1ele.eq.0) then
              kele(n)=kele(n)-ktherm*bele(nodeP,ia,i)*bele(nodeQ,ia,i)*
     1                           detele(ia)*ww(ia)/(T0/T0)
              else
              kele(n)=kele(n)-kthermcop*bele(nodeP,ia,i)*
     1                        bele(nodeQ,ia,i)*detele(ia)*ww(ia)/(T0/T0)
              endif
 505        continue
 405      continue
 305    continue

        do 311 nodeP=1,nnds
          kdofP=nodeP*ndpn !T
          do 411 nodeQ=1,nnds
            do 511 j=1,ndim
              kdofQ=(nodeQ-1)*ndpn+j !u
              if (kdofQ.gt.kdofP) goto 511
              n=kdofP*(kdofP-1)/2+kdofQ
              do 611 k=1,ndim
                kele(n)=kele(n)+deijdT(j,k)*bele(nodeQ,ia,k)*
     1                             nele(nodeP,ia)*detele(ia)*ww(ia)/dt
 611          continue
 511        continue
 411      continue
 311    continue

 200  continue

      do 251 i=1,nnds*ndpn*(nnds*ndpn+1)/2
        if (isnan(kele(i))) print*,iele,i,kele(i)
!        if (abs(kele(i)).gt.1.d-10) print*,iele,i,kele(i)
 251  continue

      return
      end

******************************************************************************

      subroutine materialEC(exx,eyy,ezz,exy,eyz,ezx,Px,Py,Pz,T,a0,ndim,
     1               deijdekl,deijdPk,deijdT,dEidEj,dEidPj,dPidPk,dPidT,
     2               dxiijdxikl,dTdT,ktherm,T0,beta_P,dt)
      implicit none
      real*8 exx,eyy,ezz,exy,eyz,ezx,Px,Py,Pz,T
      real*8 c11,c12,c44
      real*8 a0,a1c,a12,a11,a111,a112,a123,a1111,a1112,a1122,a1123
      real*8 Q11,Q12,Q44,Q11n
      real*8 T0,P0,ee0,E0,s0,dt,l0,beta_P
      real*8 k0,beta1,cep,Ta,ktherm
      real*8 deijdekl(ndim,ndim,ndim,ndim),deijdPk(ndim,ndim,ndim)
      real*8 deijdT(ndim,ndim),dEidEj(ndim,ndim),dEidPj(ndim,ndim)
      real*8 dPidPk(ndim,ndim),dPidT(ndim)
      real*8 dxiijdxikl(ndim,ndim,ndim,ndim),dTdT
      integer i,j,k,l,ndim

*** electrocaloric material
! CURIE TEMPERATURE
      T0=273.d0+124.99725d0

! VALUE OF P AT FREE ENERGY MIN IN TETRAGONAL STATE
      P0=0.181144132267376282d0

! FROM CHEN APPENDIX
      Q11n=0.10d0
! Spontaneous strain = Q11*spontaneous polarization^2
      ee0=Q11n*P0**2.d0

! SPONTANEOUS ELECTRIC FIELD FROM d_PHI/d_P VS E PLOT
      E0=3.205504976500d5

! sigma_o, a derived quantity used for normalization
      s0=E0*P0/ee0

! characteristic length scale (for normalization)
      l0=sqrt(a0*P0/E0)

! Permittivity of free space (known value)
      k0=-8.854d-12*E0/P0

! TEMPERATURE DEPENDENT TERMS
! from li et al. 2005
      a1c=4.124d5*T0*P0/E0 !constant of a1 term (see above comment)
      a11=(-2.097d8)*P0**3.d0/E0
      a12=(7.974d8)*P0**3.d0/E0
      a111=(1.294d9)*P0**5.d0/E0
      a112=-1.950d9*P0**5.d0/E0
      a123=-2.500d9*P0**5.d0/E0
      a1111=3.863d10*P0**7.d0/E0
      a1112=2.529d10*P0**7.d0/E0
      a1122=1.637d10*P0**7.d0/E0
      a1123=1.367d10*P0**7.d0/E0

! Chen Appendix numbers
      c11=1.78d11*ee0/s0
      c12=0.964d11*ee0/s0
      c44=1.22d11*ee0/s0

      Q11=0.10d0*P0**2.d0/ee0
      Q12=-0.034d0*P0**2.d0/ee0
      Q44=0.029d0*P0**2.d0/ee0

! Heat capacity
      cep=2.5324d6*T0/(E0*P0) !2.5324×106 J/(m3 K)

! Thermal expansion (=0 for now)
      beta1=3.d-6*T0/s0

! Some reference temperature for beta1 and cep
      Ta=1.d0*T0/T0

! Thermal conductivity
      ktherm=2.61d0*dt*T0/(P0*E0*l0**2)

! Second derivatives of h
      dTdT=0.d0
      do 100 i=1,ndim
        dPidT(i)=0.d0
        do 200 j=1,ndim
          deijdT(i,j)=0.d0
          dEidEj(i,j)=0.d0
          dEidPj(i,j)=0.d0
          dPidPk(i,j)=0.d0
          do 300 k=1,ndim
            deijdPk(i,j,k)=0.d0
            do 400 l=1,ndim
              deijdekl(i,j,k,l)=0.d0
              dxiijdxikl(i,j,k,l)=0.d0
 400        continue
 300      continue
 200    continue
 100  continue

      deijdekl(1,1,1,1)=c11
      deijdekl(2,2,2,2)=deijdekl(1,1,1,1)
      deijdekl(3,3,3,3)=deijdekl(1,1,1,1)

      deijdekl(1,1,2,2)=c12
      deijdekl(1,1,3,3)=deijdekl(1,1,2,2)
      deijdekl(2,2,1,1)=deijdekl(1,1,2,2)
      deijdekl(2,2,3,3)=deijdekl(1,1,2,2)
      deijdekl(3,3,1,1)=deijdekl(1,1,2,2)
      deijdekl(3,3,2,2)=deijdekl(1,1,2,2)

      deijdekl(1,2,1,2)=c44
      deijdekl(1,2,2,1)=deijdekl(1,2,1,2)
      deijdekl(1,3,1,3)=deijdekl(1,2,1,2)
      deijdekl(1,3,3,1)=deijdekl(1,2,1,2)
      deijdekl(2,3,2,3)=deijdekl(1,2,1,2)
      deijdekl(2,3,3,2)=deijdekl(1,2,1,2)
      deijdekl(2,1,2,1)=deijdekl(1,2,1,2)
      deijdekl(2,1,1,2)=deijdekl(1,2,1,2)
      deijdekl(3,1,3,1)=deijdekl(1,2,1,2)
      deijdekl(3,1,1,3)=deijdekl(1,2,1,2)
      deijdekl(3,2,3,2)=deijdekl(1,2,1,2)
      deijdekl(3,2,2,3)=deijdekl(1,2,1,2)

      deijdPk(1,1,1)=-2.d0*Px*(c11*Q11+2.d0*c12*Q12)
      deijdPk(1,1,2)=-2.d0*Py*(c11*Q12+c12*(Q11+Q12))
      deijdPk(1,1,3)=-2.d0*Pz*(c11*Q12+c12*(Q11+Q12))
      deijdPk(1,2,1)=-c44*Py*Q44
      deijdPk(1,2,2)=-c44*Px*Q44
      deijdPk(1,3,1)=-c44*Pz*Q44
      deijdPk(1,3,3)=-c44*Px*Q44
      deijdPk(2,1,1)=deijdPk(1,2,1)
      deijdPk(2,1,2)=deijdPk(1,2,2)
      deijdPk(2,2,1)=-2.d0*Px*(c11*Q12+c12*(Q11+Q12))
      deijdPk(2,2,2)=-2.d0*Py*(c11*Q11+2.d0*c12*Q12)
      deijdPk(2,2,3)=deijdPk(1,1,3)
      deijdPk(2,3,2)=-c44*Pz*Q44
      deijdPk(2,3,3)=-c44*Py*Q44
      deijdPk(3,1,1)=deijdPk(1,3,1)
      deijdPk(3,1,3)=deijdPk(1,3,3)
      deijdPk(3,2,2)=deijdPk(2,3,2)
      deijdPk(3,2,3)=deijdPk(2,3,3)
      deijdPk(3,3,1)=deijdPk(2,2,1)
      deijdPk(3,3,2)=deijdPk(1,1,2)
      deijdPk(3,3,3)=-2.d0*Pz*(c11*Q11+2.d0*c12*Q12)

      deijdT(1,1)=beta1
      deijdT(2,2)=beta1
      deijdT(3,3)=beta1

      dEidEj(1,1)=-k0
      dEidEj(2,2)=-k0
      dEidEj(3,3)=-k0

      dEidPj(1,1)=-1.d0
      dEidPj(2,2)=-1.d0
      dEidPj(3,3)=-1.d0

      dPidPk(1,1)=(c44*Py**2*Q44**2*T0 + c44*Pz**2*Q44**2*T0 +
     1    12.d0*Px**2*T0*a11 + 30.d0*Px**4*T0*a111 +
     2    56.d0*Px**6*T0*a1111 + 30.d0*Px**4*Py**2*T0*a1112 +
     3    2.d0*Py**6*T0*a1112 + 30.d0*Px**4*Pz**2*T0*a1112 +
     4    2.d0*Pz**6*T0*a1112 + 12.d0*Px**2*Py**2*T0*a112 +
     5    2.d0*Py**4*T0*a112 + 12.d0*Px**2*Pz**2*T0*a112 +
     6    2.d0*Pz**4*T0*a112 + 12.d0*Px**2*Py**4*T0*a1122 +
     7    12.d0*Px**2*Pz**4*T0*a1122 +
     8    12.d0*Px**2*Py**2*Pz**2*T0*a1123 + 2.d0*Py**4*Pz**2*T0*a1123 +
     9    2.d0*Py**2*Pz**4*T0*a1123 + 2.d0*Py**2*T0*a12 +
     a    2.d0*Pz**2*T0*a12 + 2.d0*Py**2*Pz**2*T0*a123 -
     b    776.d0*a1c + 2.d0*T*T0*a1c +
     c    2.d0*c11*T0*(2.d0*Pz**2*Q11*Q12 + Pz**2*Q12**2 +
     d       Py**2*Q12*(2.d0*Q11 + Q12) +
     e       3.d0*Px**2*(Q11**2 + 2.d0*Q12**2) - Q11*exx -
     f       Q12*eyy - Q12*ezz) +
     g    2.d0*c12*T0*(12.d0*Px**2*Q11*Q12 + 6.d0*Px**2*Q12**2 +
     h       Py**2*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) +
     i       Pz**2*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) -
     j       2.d0*Q12*exx - Q11*eyy - Q12*eyy - Q11*ezz - Q12*ezz))/T0
      dPidPk(1,2)= 2.d0*(2.d0*c11*Px*Py*Q12*(2.d0*Q11 + Q12) +
     1    2.d0*c12*Px*Py*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) +
     2    c44*Px*Py*Q44**2 + 6.d0*Px**5*Py*a1112 +
     3    6.d0*Px*Py**5*a1112 + 4.d0*Px**3*Py*a112 +
     4    4.d0*Px*Py**3*a112 + 8.d0*Px**3*Py**3*a1122 +
     5    4.d0*Px**3*Py*Pz**2*a1123 + 4.d0*Px*Py**3*Pz**2*a1123 +
     6    2.d0*Px*Py*Pz**4*a1123 + 2.d0*Px*Py*a12 +
     7    2.d0*Px*Py*Pz**2*a123 - c44*Q44*exy)
      dPidPk(2,1)=dPidPk(1,2)
      dPidPk(1,3)=2.d0*(2.d0*c11*Px*Pz*Q12*(2.d0*Q11 + Q12) +
     1    2.d0*c12*Px*Pz*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) +
     2    c44*Px*Pz*Q44**2 + 6.d0*Px**5*Pz*a1112 +
     3    6.d0*Px*Pz**5*a1112 + 4.d0*Px**3*Pz*a112 +
     4    4.d0*Px*Pz**3*a112 + 8.d0*Px**3*Pz**3*a1122 +
     5    4.d0*Px**3*Py**2*Pz*a1123 + 2.d0*Px*Py**4*Pz*a1123 +
     6    4.d0*Px*Py**2*Pz**3*a1123 + 2.d0*Px*Pz*a12 +
     7    2.d0*Px*Py**2*Pz*a123 - c44*Q44*ezx)
      dPidPk(3,1)=dPidPk(1,3)
      dPidPk(2,2)=(c44*Px**2*Q44**2*T0 + c44*Pz**2*Q44**2*T0 +
     1    12.d0*Py**2*T0*a11 + 30.d0*Py**4*T0*a111 +
     2    56.d0*Py**6*T0*a1111 + 2.d0*Px**6*T0*a1112 +
     3    30.d0*Px**2*Py**4*T0*a1112 +
     4    30.d0*Py**4*Pz**2*T0*a1112 + 2.d0*Pz**6*T0*a1112 +
     5    2.d0*Px**4*T0*a112 + 12.d0*Px**2*Py**2*T0*a112 +
     6    12.d0*Py**2*Pz**2*T0*a112 + 2.d0*Pz**4*T0*a112 +
     7    12.d0*Px**4*Py**2*T0*a1122 + 12.d0*Py**2*Pz**4*T0*a1122 +
     8    2.d0*Px**4*Pz**2*T0*a1123 + 12.d0*Px**2*Py**2*Pz**2*T0*a1123 +
     9    2.d0*Px**2*Pz**4*T0*a1123 + 2.d0*Px**2*T0*a12 +
     a    2.d0*Pz**2*T0*a12 + 2.d0*Px**2*Pz**2*T0*a123 -
     b    776.d0*a1c + 2.d0*T*T0*a1c +
     c    2.d0*c11*T0*(2.d0*Pz**2*Q11*Q12 + Pz**2*Q12**2 +
     d       Px**2*Q12*(2.d0*Q11 + Q12) +
     e       3.d0*Py**2*(Q11**2 + 2.d0*Q12**2) - Q12*exx -
     f       Q11*eyy - Q12*ezz) +
     g    2.d0*c12*T0*(12.d0*Py**2*Q11*Q12 + 6.d0*Py**2*Q12**2 +
     h       Px**2*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) +
     i       Pz**2*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) -
     j       Q11*exx - Q12*exx - 2.d0*Q12*eyy - Q11*ezz - Q12*ezz))/T0
        dPidPk(2,3)=2.d0*(2.d0*c11*Py*Pz*Q12*(2.d0*Q11 + Q12) +
     1    2.d0*c12*Py*Pz*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) +
     2    c44*Py*Pz*Q44**2 + 6.d0*Py**5*Pz*a1112 +
     3    6.d0*Py*Pz**5*a1112 + 4.d0*Py**3*Pz*a112 +
     4    4.d0*Py*Pz**3*a112 + 8.d0*Py**3*Pz**3*a1122 +
     5    2.d0*Px**4*Py*Pz*a1123 + 4.d0*Px**2*Py**3*Pz*a1123 +
     6    4.d0*Px**2*Py*Pz**3*a1123 + 2.d0*Py*Pz*a12 +
     7    2.d0*Px**2*Py*Pz*a123 - c44*Q44*eyz)
      dPidPk(3,2)=dPidPk(2,3)
      dPidPk(3,3)=(c44*Px**2*Q44**2*T0 + c44*Py**2*Q44**2*T0 +
     1    12.d0*Pz**2*T0*a11 + 30.d0*Pz**4*T0*a111 +
     2    56.d0*Pz**6*T0*a1111 + 2.d0*Px**6*T0*a1112 +
     3    2.d0*Py**6*T0*a1112 + 30.d0*Px**2*Pz**4*T0*a1112 +
     4    30.d0*Py**2*Pz**4*T0*a1112 + 2.d0*Px**4*T0*a112 +
     5    2.d0*Py**4*T0*a112 + 12.d0*Px**2*Pz**2*T0*a112 +
     6    12.d0*Py**2*Pz**2*T0*a112 + 12.d0*Px**4*Pz**2*T0*a1122 +
     7    12.d0*Py**4*Pz**2*T0*a1122 + 2.d0*Px**4*Py**2*T0*a1123 +
     8    2.d0*Px**2*Py**4*T0*a1123 + 12.d0*Px**2*Py**2*Pz**2*T0*a1123 +
     9    2.d0*Px**2*T0*a12 + 2.d0*Py**2*T0*a12 +
     a    2.d0*Px**2*Py**2*T0*a123 - 776.d0*a1c + 2.d0*T*T0*a1c +
     b    2.d0*c11*T0*(2.d0*Py**2*Q11*Q12 + Py**2*Q12**2 +
     c       Px**2*Q12*(2.d0*Q11 + Q12) +
     d       3.d0*Pz**2*(Q11**2 + 2.d0*Q12**2) - Q12*exx -
     e       Q12*eyy - Q11*ezz) +
     f    2.d0*c12*T0*(12.d0*Pz**2*Q11*Q12 + 6.d0*Pz**2*Q12**2 +
     g       Px**2*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) +
     h       Py**2*(Q11**2 + 2.d0*Q11*Q12 + 3.d0*Q12**2) -
     i       Q11*exx - Q12*exx - Q11*eyy - Q12*eyy - 2.d0*Q12*ezz))/T0

      dPidT(1)=2.d0*Px*a1c
      dPidT(2)=2.d0*Py*a1c
      dPidT(3)=2.d0*Pz*a1c

      dxiijdxikl(1,1,1,1)=a0
      dxiijdxikl(1,2,1,2)=a0
      dxiijdxikl(2,1,2,1)=a0
      dxiijdxikl(2,2,2,2)=a0
      dxiijdxikl(3,1,3,1)=a0
      dxiijdxikl(3,2,3,2)=a0
      dxiijdxikl(1,3,1,3)=a0
      dxiijdxikl(2,3,2,3)=a0
      dxiijdxikl(3,3,3,3)=a0

      dTdT=-cep/T

      return
      end

******************************************************************************

      subroutine materialcop(T,ndim,deijdekl,deijdPk,deijdT,dEidEj,
     1                    dEidPj,dPidPk,dPidT,dxiijdxikl,
     2                    dTdT,kthermcop,T0,a0,dt)
      implicit none
      real*8 T
      real*8 T0,P0,ee0,E0,s0,dt,l0,ezz,k0,a0,Q11n
      real*8 lambda,mu
      real*8 betacop,cpcop,kthermcop
      real*8 deijdekl(ndim,ndim,ndim,ndim),deijdPk(ndim,ndim,ndim)
      real*8 deijdT(ndim,ndim),dEidEj(ndim,ndim),dEidPj(ndim,ndim)
      real*8 dPidPk(ndim,ndim),dPidT(ndim)
      real*8 dxiijdxikl(ndim,ndim,ndim,ndim),dTdT
      integer i,j,k,l,ndim

! CURIE TEMPERATURE
      T0=273.d0+124.99725d0

! VALUE OF P AT FREE ENERGY MIN IN TETRAGONAL STATE
      P0=0.181144132267376282d0

! FROM CHEN APPENDIX
      Q11n=0.10d0
! Spontaneous strain = Q11*spontaneous polarization^2
      ee0=Q11n*P0**2.d0

! SPONTANEOUS ELECTRIC FIELD FROM d_PHI/d_P VS E PLOT
      E0=3.205504976500d5

! sigma_o, a derived quantity used for normalization
      s0=E0*P0/ee0

! characteristic length scale (for normalization)
      l0=sqrt(a0*P0/E0)

! set value for generalized plane strain
      ezz=0.d0

! Permittivity of free space (known value)
      k0=-8.854d-12*E0/P0

*** electrode material
      lambda=71.d9*ee0/s0
      mu=33.4d9*ee0/s0
      betacop=0.d0!18.d-6*T0/s0
      cpcop=3.40252d6*T0/(E0*P0)
      kthermcop=401.d0*dt*T0/(P0*E0*l0**2)

! Second derivatives of h (set all to zero to start)
      dTdT=0.d0
      do 100 i=1,ndim
        dPidT(i)=0.d0
        do 200 j=1,ndim
          deijdT(i,j)=0.d0
          dEidEj(i,j)=0.d0
          dEidPj(i,j)=0.d0
          dPidPk(i,j)=0.d0
          do 300 k=1,ndim
            deijdPk(i,j,k)=0.d0
            do 400 l=1,ndim
              deijdekl(i,j,k,l)=0.d0
              dxiijdxikl(i,j,k,l)=0.d0
 400        continue
 300      continue
 200    continue
 100  continue

      deijdekl(1,1,1,1)=lambda+2.d0*mu
      deijdekl(1,1,2,2)=lambda
      deijdekl(1,2,1,2)=mu
      deijdekl(1,2,2,1)=mu
      deijdekl(2,1,1,2)=mu
      deijdekl(2,1,2,1)=mu
      deijdekl(2,2,1,1)=lambda
      deijdekl(2,2,2,2)=lambda+2.d0*mu

      deijdT(1,1)=betacop
      deijdT(2,2)=betacop

      dEidEj(1,1)=-k0
      dEidEj(2,2)=-k0

      dTdT=-cpcop/T

      return
      end

************************************************************************

      subroutine nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     1                    ele,E,D,s,strn,strs,svm,gradTheta,qt,nE,nD,nS,
     2                    nstrn,nstrs,ngradT,nflux,nsvm,nmat,EC0cop1)
      implicit none
      integer maxnd,maxele,ndim,nodes,neles,nints,nnds
      integer ele(maxele,8),i,j,k,l,ele_per_node(maxnd)
      real*8 E(maxele,nints,ndim),D(maxele,nints,ndim),s(maxele,nints)
      real*8 strn(maxele,nints,ndim,ndim)
      real*8 strs(maxele,nints,ndim+1,ndim+1),svm(maxele,nints)
      real*8 gradTheta(maxele,nints,ndim),qt(maxele,nints,ndim)
      integer EC0cop1(maxele)
      real*8 nE(maxnd,ndim),nD(maxnd,ndim),nS(maxnd)
      real*8 nstrn(maxnd,ndim,ndim),nstrs(maxnd,ndim,ndim),nsvm(maxnd)
      real*8 ngradT(maxnd,ndim),nflux(maxnd,ndim)
      real*8 nmat(maxnd)

      do 1 i=1,nodes

        ele_per_node(i)=0
        nS(i)=0.d0
        nsvm(i)=0.d0
        nmat(i)=0.d0
        do 3 k=1,ndim
          nE(i,k)=0.d0
          nD(i,k)=0.d0
          ngradT(i,k)=0.d0
          nflux(i,k)=0.d0
          do 4 l=1,ndim
            nstrn(i,k,l)=0.d0
            nstrs(i,k,l)=0.d0
 4        continue
 3      continue

        do 2 j=1,neles
          if (i.eq.ele(j,1)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+s(j,1)
            nsvm(i)=nsvm(i)+svm(j,1)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 11 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,1,k)
              nD(i,k)=nD(i,k)+D(j,1,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,1,k)
              nflux(i,k)=nflux(i,k)+qt(j,1,k)
              do 21 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,1,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,1,k,l)
 21           continue
 11         continue
          elseif (i.eq.ele(j,2)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+s(j,5)
            nsvm(i)=nsvm(i)+svm(j,5)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 12 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,5,k)
              nD(i,k)=nD(i,k)+D(j,5,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,5,k)
              nflux(i,k)=nflux(i,k)+qt(j,5,k)
              do 22 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,5,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,5,k,l)
 22           continue
 12         continue
          elseif (i.eq.ele(j,3)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+s(j,7)
            nsvm(i)=nsvm(i)+svm(j,7)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 13 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,7,k)
              nD(i,k)=nD(i,k)+D(j,7,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,7,k)
              nflux(i,k)=nflux(i,k)+qt(j,7,k)
              do 23 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,7,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,7,k,l)
 23           continue
 13         continue
          elseif (i.eq.ele(j,4)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+S(j,3)
            nsvm(i)=nsvm(i)+svm(j,3)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 14 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,3,k)
              nD(i,k)=nD(i,k)+D(j,3,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,3,k)
              nflux(i,k)=nflux(i,k)+qt(j,3,k)
              do 24 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,3,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,3,k,l)
 24           continue
 14         continue
          elseif (i.eq.ele(j,5)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+s(j,2)
            nsvm(i)=nsvm(i)+svm(j,2)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 15 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,2,k)
              nD(i,k)=nD(i,k)+D(j,2,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,2,k)
              nflux(i,k)=nflux(i,k)+qt(j,2,k)
              do 25 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,2,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,2,k,l)
 25           continue
 15         continue
          elseif (i.eq.ele(j,6)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+s(j,6)
            nsvm(i)=nsvm(i)+svm(j,6)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 16 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,6,k)
              nD(i,k)=nD(i,k)+D(j,6,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,6,k)
              nflux(i,k)=nflux(i,k)+qt(j,6,k)
              do 26 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,6,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,6,k,l)
 26           continue
 16         continue
          elseif (i.eq.ele(j,7)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+s(j,8)
            nsvm(i)=nsvm(i)+svm(j,8)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 17 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,8,k)
              nD(i,k)=nD(i,k)+D(j,8,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,8,k)
              nflux(i,k)=nflux(i,k)+qt(j,8,k)
              do 27 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,8,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,8,k,l)
 27           continue
 17         continue
          elseif (i.eq.ele(j,8)) then
            ele_per_node(i)=ele_per_node(i)+1
            nS(i)=nS(i)+s(j,4)
            nsvm(i)=nsvm(i)+svm(j,4)
            nmat(i)=nmat(i)+dble(EC0cop1(j))
            do 18 k=1,ndim
              nE(i,k)=nE(i,k)+E(j,4,k)
              nD(i,k)=nD(i,k)+D(j,4,k)
              ngradT(i,k)=ngradT(i,k)+gradTheta(j,4,k)
              nflux(i,k)=nflux(i,k)+qt(j,4,k)
              do 28 l=1,ndim
                nstrn(i,k,l)=nstrn(i,k,l)+strn(j,4,k,l)
                nstrs(i,k,l)=nstrs(i,k,l)+strs(j,4,k,l)
 28           continue
 18         continue
          endif

 2      continue

        nS(i)=nS(i)/ele_per_node(i)
        nsvm(i)=nsvm(i)/ele_per_node(i)
        nmat(i)=nmat(i)/ele_per_node(i)
!        if (abs(nS(i)).lt.1.d-16) nS(i)=0.d0
        do 19 k=1,ndim
          nE(i,k)=nE(i,k)/ele_per_node(i)
!          if (abs(nE(i,k)).lt.1.d-16) nE(i,k)=0.d0
          nD(i,k)=nD(i,k)/ele_per_node(i)
!          if (abs(nD(i,k)).lt.1.d-16) nD(i,k)=0.d0
          ngradT(i,k)=ngradT(i,k)/ele_per_node(i)
!          if (abs(ngradT(i,k)).lt.1.d-16) ngradT(i,k)=0.d0
          nflux(i,k)=nflux(i,k)/ele_per_node(i)
!          if (abs(nflux(i,k)).lt.1.d-16) nflux(i,k)=0.d0
          do 29 l=1,ndim
            nstrn(i,k,l)=nstrn(i,k,l)/ele_per_node(i)
!            if (abs(nstrn(i,k,l)).lt.1.d-16) nstrn(i,k,l)=0.d0
            nstrs(i,k,l)=nstrs(i,k,l)/ele_per_node(i)
!            if (abs(nstrs(i,k,l)).lt.1.d-16) nstrs(i,k,l)=0.d0
 29       continue
 19     continue

 1    continue

      return
      end

*****************************************************************

      subroutine writeresults(fnum,ou,ibcs,ndim,maxnd,maxele,nints,nnds,
     1                nodes,neles,ele,node,u,phi,P,nE,nD,nS,nstrn,nstrs,
     2                ngradT,nflux,beta_P,exx,eyy,avg_stress,
     3                avg_E,avg_D,avg_S,avg_T,PsiTot,theta,nsvm,EC0cop1,
     4                nmat)

      implicit none
      integer maxnd,maxele,ndim,nodes,neles,nints,nnds,ibcs
      integer ele(maxele,8),i,resultsfile
      real*8 node(maxnd,ndim),u(maxnd,ndim),phi(maxnd),P(3,maxnd,ndim)
      real*8 nE(maxnd,ndim),nD(maxnd,ndim),nS(maxnd),nsvm(maxnd)
      real*8 nstrn(maxnd,ndim,ndim),nstrs(maxnd,ndim,ndim)
      real*8 ngradT(maxnd,ndim),nflux(maxnd,ndim),nmat(maxnd)
      real*8 avg_stress(ndim,ndim),avg_E(ndim),avg_D(ndim),avg_S,avg_T
      real*8 beta_P,exx,eyy,PsiTot,theta(maxnd)
      real*8 pphi,uux,uuy,uuz,PPx,PPy,PPz
      integer EC0cop1(maxele)
      character(7) fnum
      character(4) ou

******************************************************************
****************   ParaView Post-Processing   ********************

          write(fnum,915)ibcs

          resultsfile=987
          open(unit=resultsfile,file='results'// ou //'/sol'
     1               //trim(adjustl(fnum)))

          write(resultsfile,*)'ux,uy,uz,phi,Px,Py,Pz,theta'
          write(resultsfile,*)'beta T/P ',beta_P," / ",beta_P
          write(resultsfile,*)'exx ',exx,'eyy ',eyy,'(set by ICs)'
          write(resultsfile,*)'Avg. Stress',avg_stress(1,1),
     1              avg_stress(2,2),avg_stress(3,3),
     2              avg_stress(1,2),avg_stress(2,3),avg_stress(3,1)
          write(resultsfile,*)'Avg. E',avg_E(1),avg_E(2),avg_E(3)
          write(resultsfile,*)'Avg. D',avg_D(1),avg_D(2),avg_D(3)
          write(resultsfile,*)'Avg. S',avg_S
          write(resultsfile,*)'Avg. T',avg_T
          write(resultsfile,*)'Psi',PsiTot
          write(resultsfile,*)
          do 851 i=1,nodes
            write(resultsfile,916) u(i,1),u(i,2),u(i,3),phi(i),
     1           P(3,i,1),P(3,i,2),P(3,i,3),theta(i),P(2,i,1),P(2,i,2),
     2           P(2,i,3)
 851      continue
          close(unit=resultsfile)

          open(unit=20,file='paraview'// ou //'/sol'
     1                      //trim(adjustl(fnum))//'.vtk')
          write(20,910) '# vtk DataFile Version 3.0'
          write(20,910) 'electrocaloric device'
          write(20,910) 'ASCII'
          write(20,910) 'DATASET UNSTRUCTURED_GRID'
          write(20,911) 'POINTS ',nodes,' float'
          do 201 i=1,nodes
           write(20,*)node(i,1),node(i,2),node(i,3)
 201      continue
          write(20,913) 'CELLS ',neles,neles*9
          do 202 i=1,neles
            write(20,*)8,ele(i,1)-1,ele(i,2)-1,ele(i,3)-1,
     1        ele(i,4)-1,ele(i,5)-1,ele(i,6)-1,ele(i,7)-1,ele(i,8)-1
 202      continue
          write(20,914) 'CELL_TYPES ',neles
          do 203 i=1,neles
            write(20,*)12
 203      continue
          write(20,914) 'POINT_DATA ',nodes
          write(20,910) 'VECTORS Polarization float'
          do 204 i=1,nodes
            if (abs(P(3,i,1)).lt.1.d-16) then
              PPx=0.d0
            else
              PPx=P(3,i,1)
            endif
            if (abs(P(3,i,2)).lt.1.d-16) then
              PPy=0.d0
            else
              PPy=P(3,i,2)
            endif
            if (abs(P(3,i,3)).lt.1.d-16) then
              PPz=0.d0
            else
              PPz=P(3,i,3)
            endif
            write(20,*)PPx,PPy,PPz
 204      continue
          write(20,910) 'SCALARS Potential float 1'
          write(20,910) 'LOOKUP_TABLE default'
          do 205 i=1,nodes
           if (abs(phi(i)).lt.1.d-16) then
             pphi=0.d0
           else
             pphi=phi(i)
           endif
           write(20,*)pphi
 205      continue
          write(20,910) 'SCALARS Temperature float 1'
          write(20,910) 'LOOKUP_TABLE default'
          do 211 i=1,nodes
           write(20,*)theta(i)
 211      continue
          write(20,910) 'VECTORS Displacement float'
          do 206 i=1,nodes
            if (abs(u(i,1)).lt.1.d-16) then
              uux=0.d0
            else
              uux=u(i,1)
            endif
            if (abs(u(i,2)).lt.1.d-16) then
              uuy=0.d0
            else
              uuy=u(i,2)
            endif
            if (abs(u(i,3)).lt.1.d-16) then
              uuz=0.d0
            else
              uuz=u(i,3)
            endif
            write(20,*)uux,uuy,uuz
 206      continue
          write(20,910) 'VECTORS E float'
          do 207 i=1,nodes
            write(20,*)nE(i,1),nE(i,2),nE(i,3)
 207      continue
          write(20,910) 'VECTORS D float'
          do 208 i=1,nodes
            write(20,*)nD(i,1),nD(i,2),nD(i,3)
 208      continue
          write(20,910) 'SCALARS S float 1'
          write(20,910) 'LOOKUP_TABLE default'
          do 212 i=1,nodes
            write(20,*)nS(i)
 212      continue
          write(20,910) 'TENSORS Strain float'
          do 209 i=1,nodes
            write(20,*)nstrn(i,1,1),(nstrn(i,1,2)+nstrn(i,2,1))/2.d0,
     1            (nstrn(i,1,3)+nstrn(i,3,1))/2.d0
            write(20,*)(nstrn(i,1,2)+nstrn(i,2,1))/2.d0,
     1            nstrn(i,2,2),(nstrn(i,2,3)+nstrn(i,3,2))/2.d0
            write(20,*)(nstrn(i,1,3)+nstrn(i,3,1))/2.d0,
     1            (nstrn(i,2,3)+nstrn(i,3,2))/2.d0,nstrn(i,3,3)
            write(20,*)
 209      continue
          write(20,910) 'TENSORS Stress float'
          do 210 i=1,nodes
           write(20,*)nstrs(i,1,1),nstrs(i,1,2),nstrs(i,1,3)
           write(20,*)nstrs(i,2,1),nstrs(i,2,2),nstrs(i,2,3)
           write(20,*)nstrs(i,3,1),nstrs(i,3,2),nstrs(i,3,3)
           write(20,*)
 210      continue
          write(20,910) 'SCALARS vonMises float 1'
          write(20,910) 'LOOKUP_TABLE default'
          do 215 i=1,nodes
           write(20,*)nsvm(i)
 215      continue
          write(20,910) 'SCALARS material float 1'
          write(20,910) 'LOOKUP_TABLE default'
          do 216 i=1,nodes
            write(20,*)nmat(i)
 216      continue
          write(20,910) 'VECTORS gradT float'
          do 214 i=1,nodes
            write(20,*)ngradT(i,1),ngradT(i,2),ngradT(i,3)
 214      continue
          write(20,910) 'VECTORS HeatFlux float'
          do 213 i=1,nodes
            write(20,*)nflux(i,1),nflux(i,2),ngradT(i,3)
 213      continue


 910      format(A)
 911      format(A,I7,A)
 913      format(A,I7, ' ',I7)
 914      format(A,I7)
 915      format(I7)
 916      format(F16.10,F16.10,F16.10,F20.8,F16.10,F16.10,F16.10,F16.10,
     1                             F16.10,F16.10,F16.10)

          close(unit=20)

      return
      end

*****************************************************************

      subroutine initialState(opt,w,Pz,es,exx,eyy,ezz,theta_init,trange,
     1           node,u,phi,P,theta,maxnd,ndim,nodes,rank)
      implicit none
      real*8 u(maxnd,ndim),phi(maxnd),P(3,maxnd,ndim),theta(maxnd)
      real*8 node(maxnd,ndim)
      real*8 w,Pz,es,exx,eyy,ezz,theta_init
      real*8 trange
      real*8 pola1,pola2,pola3
      integer maxnd,ndim,nodes,i,opt,rank

      do 20 i=1,nodes

! Single domain in paraelectric state
        if (opt.eq.1) then
           u(i,1)=node(i,1)*es
           u(i,2)=node(i,2)*es
           u(i,3)=node(i,3)*es
           phi(i)=0.d0
           P(3,i,1)=1.d-5
           P(3,i,2)=-1.d-5
           P(3,i,3)=1.d-5
           theta(i)=theta_init
        elseif (opt.eq.2) then
           u(i,1)=node(i,1)*exx
           u(i,2)=node(i,2)*eyy
           u(i,3)=node(i,3)*ezz
           phi(i)=0.d0
           CALL random_seed()
           CALL random_number(pola1) 
           CALL random_number(pola2)
           CALL random_number(pola3)          
           P(3,i,1)=pola1*0.00008-0.00004
           P(3,i,2)=pola2*0.00008-0.00004
           P(3,i,3)=pola3*0.00008-0.00004
           theta(i)=theta_init
        endif

!         P(1,i,1)=0.0
!        P(1,i,2)=0.0
!         P(1,i,3)=0.0
!         P(2,i,1)=0.0
!         P(2,i,2)=0.0
!         P(2,i,3)=0.0

        P(1,i,1)=P(3,i,1)
        P(2,i,1)=P(3,i,1)
        P(2,i,2)=P(3,i,2)
        P(1,i,2)=P(2,i,2)
        P(2,i,3)=P(3,i,3)
        P(1,i,3)=P(2,i,3)
 20   continue

      return
      end
*****************************************************************
