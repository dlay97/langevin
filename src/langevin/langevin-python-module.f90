!===================================================================================================
!   Things to feed in from python:
!       -Starting point
!           Include k copies of initial point for k runs
!           Will be serialized here, so probably load balancing is irrelevant
!       -PES/grad/hess
!       -Inverse inertia/grad/hess
!       -Friction and random force (constant)
!       -Random seed
!       -Neck values on grid, for stopping
!       -dt, maxIters, neck criteria
!       -File to output to, and whether to save the paths or not
!
!===================================================================================================

module langevin_2d
    !===============================================================================================
    !   The solver in 2d
    !===============================================================================================
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    implicit none

    real(kind=real64), dimension(:), allocatable                :: uniqueCoords1, uniqueCoords2
    real(kind=real64), dimension(:,:), allocatable              :: pes, neckVals
    real(kind=real64), dimension(:,:,:), allocatable            :: pesGrad, pesHess
    
    !Symmetric matrices are flattened, as
    !   [[M11, M12],
    !    [M12, M22]]
    !
    ! -> [M11, M12, M22]

    !Ordered as (x,y,component)
    real(kind=real64), dimension(:,:,:), allocatable            :: inverseMetric
    !Ordered as (x,y,component,componentGrad/componentHess)
    real(kind=real64), dimension(:,:,:,:), allocatable          :: inverseMetricGrad, inverseMetricHess

    !friction is $\gamma$, randomForce is $g$ with $g.g=\gamma$
    real(kind=real64), dimension(:,:), allocatable              :: friction, randomForce
    
    integer, dimension(2)                                       :: nGridpoints
    ! real(kind=real64), dimension(:,:), allocatable              :: startingPoints

    real(kind=real64), dimension(2)                             :: dx

    !Run parameters
    !TODO: rename massVal -> atomicNumber
    real(kind=real64)                                           :: massVal, startingEneg, dt, maxTime, scissionNeckVal
    real(kind=real64)                                           :: levelDensity
    integer                                                     :: maxIters

    contains
        subroutine setup_arrays(unCoords1In,unCoords2In,&
                                pesIn,pesGradIn,pesHessIn,&
                                invMetIn,invMetGradIn,invMetHessIn,&
                                fricIn,randForceIn,neckIn)

            implicit none
            !maybe a better way of doing this

            !note that, when not supplying the array size, it's important to put arrays of
            !different sizes on different lines. otherwise, it may segfault
            real(kind=real64), dimension(:), intent(in)                         :: unCoords1In
            real(kind=real64), dimension(:), intent(in)                         :: unCoords2In
            real(kind=real64), dimension(:,:), intent(in)                       :: pesIn, neckIn
            real(kind=real64), dimension(:,:,:), intent(in)                     :: pesGradIn
            real(kind=real64), dimension(:,:,:), intent(in)                     :: pesHessIn
            real(kind=real64), dimension(:,:,:), intent(in)                     :: invMetIn
            real(kind=real64), dimension(:,:,:,:), intent(in)                   :: invMetGradIn
            real(kind=real64), dimension(:,:,:,:), intent(in)                   :: invMetHessIn
            real(kind=real64), dimension(:,:), intent(in)                       :: fricIn, randForceIn

            uniqueCoords1 = unCoords1In
            uniqueCoords2 = unCoords2In
            pes = pesIn
            pesGrad = pesGradIn
            pesHess = pesHessIn
            inverseMetric = invMetIn
            inverseMetricGrad = invMetGradIn
            inverseMetricHess = invMetHessIn
            friction = fricIn
            randomForce = randForceIn
            neckVals = neckIn

            ! write(*,*) uniqueCoords1
            ! write(*,*) neckVals

            nGridpoints = shape(pes)
            
            dx(1) = uniqueCoords1(2) - uniqueCoords1(1)
            dx(2) = uniqueCoords2(2) - uniqueCoords2(1)

        end subroutine setup_arrays

        subroutine setup_params(massValIn,startingEnegIn,dtIn,maxTimeIn,scissionNeckValIn,useSeed,seed)
            implicit none

            integer, dimension(33), intent(in)                      :: seed
            logical, intent(in)                                     :: useSeed
            real(kind=real64), intent(in)                           :: massValIn,startingEnegIn,dtIn,maxTimeIn,scissionNeckValIn

            massVal = massValIn !maybe don't need to use this anywhere else in the code
            startingEneg = startingEnegIn
            dt = dtIn
            maxTime = maxTimeIn !same as you
            scissionNeckVal = scissionNeckValIn

            levelDensity = massVal/10.d0
            maxIters = int(maxTime/dt)

            if (useSeed) call random_seed(put=seed)

        end subroutine setup_params

        subroutine get_interpolation_indices(coords,idx,isOutOfBounds)
            implicit none

            integer, dimension(2), intent(out)                      :: idx
            logical, intent(out)                                    :: isOutOfBounds
            real(kind=real64), dimension(2)                         :: coords

            idx(1) = int((coords(1) - uniqueCoords1(1))/dx(1)) + 1
            idx(2) = int((coords(2) - uniqueCoords2(1))/dx(2)) + 1

            isOutOfBounds = .FALSE.

            !Set to '.gt.nGridpoints(1) - 1' because we need one gridpoint above/below/etc this point
            if (idx(1).gt.nGridpoints(1)-1) then
                idx(1) = nGridpoints(1)
                isOutOfBounds = .TRUE.
            else if (idx(1).lt.2) then
                idx(1) = 1
                isOutOfBounds = .TRUE.
            end if

            if (idx(2).gt.nGridpoints(2)-1) then
                idx(2) = nGridpoints(2)
                isOutOfBounds = .TRUE.
            else if (idx(2).lt.2) then
                idx(2) = 1
                isOutOfBounds = .TRUE.
            end if

        end subroutine get_interpolation_indices

        subroutine bilinear_interp(idx,coords,isOutOfBounds,arrIn,val)
            implicit none

            integer, dimension(2)                                   :: idx
            logical                                                 :: isOutOfBounds
            real(kind=real64), dimension(nGridpoints(1),nGridpoints(2))                  :: arrIn
            real(kind=real64), intent(out)                          :: val
            real(kind=real64)                                       :: x1, x2, y1, y2, f00, f01, f10, f11
            real(kind=real64), dimension(2)                         :: coords

            if (isOutOfBounds) then
                val = arrIn(idx(1),idx(2))
                return
            end if

            !Notation and equations from 
            !https://en.wikipedia.org/wiki/Bilinear_interpolation#Repeated_linear_interpolation
            x1 = uniqueCoords1(idx(1))
            x2 = uniqueCoords1(idx(1)+1)

            y1 = uniqueCoords2(idx(2))
            y2 = uniqueCoords2(idx(2)+1)

            f00 = arrIn(idx(1),idx(2))
            f01 = arrIn(idx(1),idx(2)+1)
            f10 = arrIn(idx(1)+1,idx(2))
            f11 = arrIn(idx(1)+1,idx(2)+1)

            val = f00*(x2-coords(1))*(y2-coords(2))
            val = val + f10*(coords(1)-x1)*(y2-coords(2))
            val = val + f01*(x2-coords(1))*(coords(2)-y1)
            val = val + f11*(coords(1)-x1)*(coords(2)-y1)
            val = val/(dx(1)*dx(2))

        end subroutine bilinear_interp

        subroutine interp_wrapper(coords,&
                                  eneg,enegGrad,enegHess,&!6
                                  invMet,invMetGrad,invMetHess,&!3+6+9=18!
                                  neckVal,isOutOfBounds)
                                  !total of 42 evals, matching the previous code, minus (3+6)*2=18 from friction/grad 
                                  !and randForce/grad
            !trust this wrapper against scipy.interpolate.interpn
            implicit none

            integer                                         :: iter1, iter2
            integer, dimension(2)                           :: idx
            logical, intent(out)                            :: isOutOfBounds
            real(kind=real64), dimension(2), intent(in)     :: coords
            real(kind=real64), intent(out)                  :: eneg, neckVal
            real(kind=real64), dimension(2), intent(out)    :: enegGrad
            real(kind=real64), dimension(3), intent(out)    :: invMet, enegHess
            real(kind=real64), dimension(3,2), intent(out)  :: invMetGrad
            real(kind=real64), dimension(3,3), intent(out)  :: invMetHess
            
            call get_interpolation_indices(coords,idx,isOutOfBounds)

            call bilinear_interp(idx,coords,isOutOfBounds,pes,eneg)
            !1 eval

            call bilinear_interp(idx,coords,isOutOfBounds,neckVals,neckVal)

            do iter1=1,2
                call bilinear_interp(idx,coords,isOutOfBounds,pesGrad(:,:,iter1),enegGrad(iter1))
            end do
            !2 evals

            do iter1=1,3
                call bilinear_interp(idx,coords,isOutOfBounds,inverseMetric(:,:,iter1),invMet(iter1))
                call bilinear_interp(idx,coords,isOutOfBounds,pesHess(:,:,iter1),enegHess(iter1))
                !12 evals (6 total removed from friction and randomForce)

                do iter2=1,2
                    call bilinear_interp(idx,coords,isOutOfBounds,inverseMetricGrad(:,:,iter1,iter2),invMetGrad(iter1,iter2))
                end do
                !18 evals (12 total removed from frictionGrad and randomForceGrad)

                do iter2=1,3
                    call bilinear_interp(idx,coords,isOutOfBounds,inverseMetricHess(:,:,iter1,iter2),invMetHess(iter1,iter2))
                end do
                !9 evals

            end do
            !total of 42 evals, as expected
            
        end subroutine interp_wrapper

        subroutine time_evolve_old(coordsIn,allCoords,allMomenta,lastIter,finishedSuccessfully)
            !==================================================================================
            !
            !   Time evolves one trajectory. Should take in:
            !       coords - the starting coords
            !       saveToFile - whether or not to save this run
            !       fileName - if save, need this for parallel i/o (may be ID instead)
            !       dsetName - if save, need this to specify the dataset number
            !
            !agrees with Daniel's python code when ran with constants instead of random numbers
            !
            !==================================================================================

            use normal
            implicit none

            integer                                                     :: dtIter, randIter
            integer, intent(out)                                        :: lastIter
            logical                                                     :: useFric, isOutOfBounds
            logical, intent(out)                                        :: finishedSuccessfully
            real(kind=real64)                                           :: eneg, excitationEnergy, &
                                                                           temperature, rndParam, neckVal
            real(kind=real64), dimension(2), intent(in)                 :: coordsIn
            real(kind=real64), dimension(2)                             :: coords, enegGrad, momentum, Gamma1, Gamma2, Gamma3
            real(kind=real64), dimension(3)                             :: invMet, enegHess
            real(kind=real64), dimension(2,2)                           :: randForce
            real(kind=real64), dimension(3,2)                           :: invMetGrad!, fricGrad, randForceGrad
            real(kind=real64), dimension(3,3)                           :: invMetHess
            real(kind=real64), dimension(4)                             :: randValues
            real(kind=real64), dimension(maxIters,2), intent(out)       :: allCoords, allMomenta

            !tensors used in calculation
            real(kind=real64), dimension(2)                            :: v, h, gGamma2
            real(kind=real64), dimension(2,2)                          :: dvdq, dvdp, dhdq, dhdp

            allCoords = 0.d0
            allMomenta = 0.d0

            finishedSuccessfully = .false.

            momentum = 0.d0
            rndParam = 1.d0/(2.d0*sqrt(3.d0)) !bad name, I know

            randIter=0
            coords = coordsIn

            lastIter = 1 !weird edge case in which scission happens at iteration 1

            do dtIter=1,maxIters
                call interp_wrapper(coords,&
                                    eneg,enegGrad,enegHess,&
                                    invMet,invMetGrad,invMetHess,&
                                    neckVal,isOutOfBounds)

                allCoords(dtIter,:) = coords
                allMomenta(dtIter,:) = momentum

                if (isOutOfBounds) then
                    write(*,*) 'Point (almost) out-of-bounds: ',coords
                    exit
                end if

                if (neckVal.lt.scissionNeckVal) then
                    finishedSuccessfully = .true.
                    write(*,*) 'Scission; stopping this path'
                    write(*,*) 'Neck value: ',neckVal
                    exit
                end if

                !trust PES and inverse metric and their derivatives against Jhilam's code
                
                ! excitationEnergy = startingEneg - 0.5*momentum(1)**2 * invMet(1) - &
                !     momentum(1)*momentum(2)*invMet(2) - 0.5*momentum(2)**2 * invMet(3) - eneg
                ! excitationEnergy = startingEneg + 0.5*momentum(1)**2 * invMet(1) + &
                !     momentum(1)*momentum(2)*invMet(2) + 0.5*momentum(2)**2 * invMet(3) - eneg
                excitationEnergy = startingEneg - eneg
                ! write(*,*) 'iter',dtIter, 'Potential',eneg,'exint',excitationEnergy

                !original code has option to just never use friction; that's been removed here
                if (excitationEnergy.le.0.) then
                    useFric = .false.
                else
                    useFric = .true.
                end if

                if (useFric) then
                    temperature = dsqrt(excitationEnergy/levelDensity)
                    
                    do randIter=1,4
                        randValues(randIter) = r8_normal_ab(0.d0,sqrt(2.d0)) 
                        !not sure rn what mean/std should be; probably same as in the python code
                    end do
                    ! write(*,*) 'randValues',randValues
                    
                    randForce = randomForce*sqrt(temperature)
                    ! fricGrad = fricGrad*sqrt(temperature)

                    Gamma1 = sqrt(dt)*randValues(1:2)
                    Gamma2 = dt**(1.5d0)*(0.5d0*randValues(1:2) + rndParam*randValues(3:4))
                    Gamma3 = dt**(1.5d0)*(0.5*randValues(1:2) - rndParam*randValues(3:4))
                else
                    Gamma1 = 0.d0
                    Gamma2 = 0.d0
                    Gamma3 = 0.d0
                    randForce = 0.d0
                    ! randForceGrad = 0.d0
                end if

                !for PSD matrix [[11,12],[21,22]], map 11->1, (12,21)->2, 22->3
                !eventually can make nice with matmul (or try to, at least)
                v(1) = invMet(1)*momentum(1) + invMet(2)*momentum(2)
                v(2) = invMet(2)*momentum(1) + invMet(3)*momentum(2)

                ! write(*,*) 'v',v

                h = -enegGrad - 0.5d0 * invMetGrad(1,:)*coords(1)**2 - invMetGrad(2,:)*coords(1)*coords(2) - &
                        0.5d0 * invMetGrad(3,:)*coords(2)**2 - matmul(friction,v)
                ! h(1) = h(1) - fric(1)*v(1) - fric(2)*v(2)
                ! h(2) = h(2) - fric(2)*v(1) - fric(3)*v(2)
                ! write(*,*) 'h',h

                dvdq(1,:) = invMetGrad(1,:)*momentum(1) + invMetGrad(2,:)*momentum(2)
                dvdq(2,:) = invMetGrad(2,:)*momentum(1) + invMetGrad(3,:)*momentum(2)

                ! write(*,*) 'dvdq',dvdq

                dvdp(1,1) = invMet(1)
                dvdp(1,2) = invMet(2)
                dvdp(2,1) = invMet(2)
                dvdp(2,2) = invMet(3)

                ! write(*,*) 'dvdp',dvdp

                dhdq(1,1) = -enegHess(1) - 0.5d0*invMetHess(1,1)*coords(1)**2 - invMetHess(2,1)*coords(1)*coords(2) - &
                                0.50*invMetHess(3,1)*coords(2)**2 - &
                                invMetGrad(1,1)*coords(1) - invMetGrad(2,1)*coords(2) - &
                                randForce(1,1)*dvdq(1,1) - randForce(1,2)*dvdq(2,1)
                dhdq(2,1) = -enegHess(2) - 0.5d0*invMetHess(1,2)*coords(1)**2 - invMetHess(2,2)*coords(1)*coords(2) - &
                                0.50*invMetHess(3,2)*coords(2)**2 - &
                                invMetGrad(1,2)*coords(1) - invMetGrad(2,2)*coords(2) - &
                                randForce(2,1)*dvdq(1,1) - randForce(2,2)*dvdq(2,1)
                dhdq(1,2) = -enegHess(2) - 0.5d0*invMetHess(1,2)*coords(1)**2 - invMetHess(2,2)*coords(1)*coords(2) - &
                                0.50*invMetHess(3,2)*coords(2)**2 - &
                                invMetGrad(2,1)*coords(1) - invMetGrad(3,1)*coords(2) - &
                                randForce(1,1)*dvdq(1,2) - randForce(1,2)*dvdq(2,2)
                dhdq(2,2) = -enegHess(3) - 0.5d0*invMetHess(1,3)*coords(1)**2 - invMetHess(2,3)*coords(1)*coords(2) - &
                                0.50*invMetHess(3,3)*coords(2)**2 - &
                                invMetGrad(2,2)*coords(1) - invMetGrad(3,2)*coords(2) - &
                                randForce(2,1)*dvdq(1,2) - randForce(2,2)*dvdq(2,2)

                ! write(*,*) 'dhdq',dhdq

                ! dhdp(1,:) = -fric(1)*dvdp(1,:) - fric(2)*dvdp(2,:)
                ! dhdp(2,:) = -fric(2)*dvdp(1,:) - fric(3)*dvdp(2,:)

                dhdp(1,:) = -friction(1,1)*dvdp(1,:) - friction(1,2)*dvdp(2,:)
                dhdp(2,:) = -friction(2,1)*dvdp(1,:) - friction(2,2)*dvdp(2,:)

                ! write(*,*) 'dhdp',dhdp

                gGamma2 = matmul(randForce,Gamma2)
                ! write(*,*) 'gGamma2',gGamma2
                ! gGamma2(1) = randForce(1)*Gamma2(1) + randForce(2)*Gamma2(2)
                ! gGamma2(2) = randForce(2)*Gamma2(1) + randForce(3)*Gamma2(2)

                !actually updating
                momentum = momentum + dt*h + 0.5d0*dt**2 * (matmul(dhdq,v) + matmul(dhdp,h)) + matmul(dhdp,gGamma2) + &
                                matmul(randForce,Gamma1)
                ! momentum(1) = momentum(1) + randForce(1)*Gamma1(1) + randForce(2)*Gamma1(2)! + &
                !                 ! randForceGrad(1,1)*v(1)*Gamma3(1) + randForceGrad(2,1)*v(1)*Gamma3(2) + &
                !                 ! randForceGrad(1,2)*v(2)*Gamma3(1) + randForceGrad(2,2)*v(2)*Gamma3(2)
                ! momentum(2) = momentum(2) + randForce(2)*Gamma1(1) + randForce(3)*Gamma1(2)! + &
                !                 ! randForceGrad(2,1)*v(1)*Gamma3(1) + randForceGrad(3,1)*v(1)*Gamma3(2) + &
                !                 ! randForceGrad(2,2)*v(2)*Gamma3(1) + randForceGrad(3,2)*v(2)*Gamma3(2)

                coords = coords + dt*v + 0.5d0*dt**2.d0*(matmul(dvdq,v) + matmul(dvdp,h)) + matmul(dvdp,gGamma2)

                ! write(*,*) coords, momentum
                ! stop

                lastIter = dtIter

            end do

        end subroutine time_evolve_old

        subroutine time_evolve(coordsIn,allCoords,allMomenta,lastIter,finishedSuccessfully)
            !==================================================================================
            !
            !   Time evolves one trajectory. Should take in:
            !       coords - the starting coords
            !       saveToFile - whether or not to save this run
            !       fileName - if save, need this for parallel i/o (may be ID instead)
            !       dsetName - if save, need this to specify the dataset number
            !
            !agrees with Daniel's python code when ran with constants instead of random numbers
            !
            !==================================================================================

            use normal
            implicit none

            integer                                                     :: dtIter, randIter
            integer, intent(out)                                        :: lastIter
            logical                                                     :: useFric, isOutOfBounds
            logical, intent(out)                                        :: finishedSuccessfully
            real(kind=real64)                                           :: eneg, excitationEnergy, &
                                                                           temperature, rndParam, neckVal
            real(kind=real64), dimension(2), intent(in)                 :: coordsIn
            real(kind=real64), dimension(2)                             :: coords, enegGrad, momentum, Gamma1, Gamma2, Gamma3
            real(kind=real64), dimension(3)                             :: invMet, enegHess
            real(kind=real64), dimension(2,2)                           :: randForce
            real(kind=real64), dimension(3,2)                           :: invMetGrad!, fricGrad, randForceGrad
            real(kind=real64), dimension(3,3)                           :: invMetHess
            real(kind=real64), dimension(4)                             :: randValues
            real(kind=real64), dimension(maxIters,2), intent(out)       :: allCoords, allMomenta

            !tensors used in calculation
            real(kind=real64), dimension(2)                            :: v, h, gGamma2
            real(kind=real64), dimension(2,2)                          :: dvdq, dvdp, dhdq, dhdp

            allCoords = 0.d0
            allMomenta = 0.d0

            finishedSuccessfully = .false.

            momentum = 0.d0
            rndParam = 1.d0/(2.d0*sqrt(3.d0)) !bad name, I know

            randIter=0
            coords = coordsIn

            lastIter = 1 !weird edge case in which scission happens at iteration 1

            do dtIter=1,maxIters
                call interp_wrapper(coords,&
                                    eneg,enegGrad,enegHess,&
                                    invMet,invMetGrad,invMetHess,&
                                    neckVal,isOutOfBounds)

                allCoords(dtIter,:) = coords
                allMomenta(dtIter,:) = momentum

                if (isOutOfBounds) then
                    write(*,*) 'Point (almost) out-of-bounds: ',coords
                    exit
                end if

                if (neckVal.lt.scissionNeckVal) then
                    finishedSuccessfully = .true.
                    write(*,*) 'Scission; stopping this path'
                    write(*,*) 'Neck value: ',neckVal
                    exit
                end if

                !trust PES and inverse metric and their derivatives against Jhilam's code
                
                ! excitationEnergy = startingEneg - 0.5*momentum(1)**2 * invMet(1) - &
                !     momentum(1)*momentum(2)*invMet(2) - 0.5*momentum(2)**2 * invMet(3) - eneg
                ! excitationEnergy = startingEneg + 0.5*momentum(1)**2 * invMet(1) + &
                !     momentum(1)*momentum(2)*invMet(2) + 0.5*momentum(2)**2 * invMet(3) - eneg
                excitationEnergy = startingEneg - eneg
                ! write(*,*) 'iter',dtIter, 'Potential',eneg,'exint',excitationEnergy

                !original code has option to just never use friction; that's been removed here
                if (excitationEnergy.le.0.) then
                    useFric = .false.
                else
                    useFric = .true.
                end if

                if (useFric) then
                    temperature = dsqrt(excitationEnergy/levelDensity)
                    
                    do randIter=1,4
                        randValues(randIter) = r8_normal_ab(0.d0,sqrt(2.d0)) 
                        !not sure rn what mean/std should be; probably same as in the python code
                    end do
                    ! write(*,*) 'randValues',randValues
                    
                    randForce = randomForce*sqrt(temperature)
                    ! fricGrad = fricGrad*sqrt(temperature)

                    Gamma1 = sqrt(dt)*randValues(1:2)
                    Gamma2 = dt**(1.5d0)*(0.5d0*randValues(1:2) + rndParam*randValues(3:4))
                    Gamma3 = dt**(1.5d0)*(0.5*randValues(1:2) - rndParam*randValues(3:4))
                else
                    Gamma1 = 0.d0
                    Gamma2 = 0.d0
                    Gamma3 = 0.d0
                    randForce = 0.d0
                    ! randForceGrad = 0.d0
                end if

                !for PSD matrix [[11,12],[21,22]], map 11->1, (12,21)->2, 22->3
                !eventually can make nice with matmul (or try to, at least)
                v(1) = invMet(1)*momentum(1) + invMet(2)*momentum(2)
                v(2) = invMet(2)*momentum(1) + invMet(3)*momentum(2)

                ! write(*,*) 'v',v

                h = -enegGrad - 0.5d0 * invMetGrad(1,:)*momentum(1)**2 - invMetGrad(2,:)*momentum(1)*momentum(2) - &
                        0.5d0 * invMetGrad(3,:)*momentum(2)**2 - matmul(friction,v)
                ! h(1) = h(1) - fric(1)*v(1) - fric(2)*v(2)
                ! h(2) = h(2) - fric(2)*v(1) - fric(3)*v(2)
                ! write(*,*) 'h',h

                dvdq(1,:) = invMetGrad(1,:)*momentum(1) + invMetGrad(2,:)*momentum(2)
                dvdq(2,:) = invMetGrad(2,:)*momentum(1) + invMetGrad(3,:)*momentum(2)

                ! write(*,*) 'dvdq',dvdq

                dvdp(1,1) = invMet(1)
                dvdp(1,2) = invMet(2)
                dvdp(2,1) = invMet(2)
                dvdp(2,2) = invMet(3)

                ! write(*,*) 'dvdp',dvdp

                dhdq(1,1) = -enegHess(1) - 0.5d0*invMetHess(1,1)*momentum(1)**2 - invMetHess(2,1)*momentum(1)*momentum(2) - &
                                0.50*invMetHess(3,1)*momentum(2)**2 - &
                                randForce(1,1)*dvdq(1,1) - randForce(1,2)*dvdq(2,1)
                dhdq(2,1) = -enegHess(2) - 0.5d0*invMetHess(1,2)*momentum(1)**2 - invMetHess(2,2)*momentum(1)*momentum(2) - &
                                0.50*invMetHess(3,2)*momentum(2)**2 - &
                                randForce(2,1)*dvdq(1,1) - randForce(2,2)*dvdq(2,1)
                dhdq(1,2) = -enegHess(2) - 0.5d0*invMetHess(1,2)*momentum(1)**2 - invMetHess(2,2)*momentum(1)*momentum(2) - &
                                0.50*invMetHess(3,2)*momentum(2)**2 - &
                                randForce(1,1)*dvdq(1,2) - randForce(1,2)*dvdq(2,2)
                dhdq(2,2) = -enegHess(3) - 0.5d0*invMetHess(1,3)*momentum(1)**2 - invMetHess(2,3)*momentum(1)*momentum(2) - &
                                0.50*invMetHess(3,3)*momentum(2)**2 - &
                                randForce(2,1)*dvdq(1,2) - randForce(2,2)*dvdq(2,2)

                ! write(*,*) 'dhdq',dhdq

                ! dhdp(1,:) = -fric(1)*dvdp(1,:) - fric(2)*dvdp(2,:)
                ! dhdp(2,:) = -fric(2)*dvdp(1,:) - fric(3)*dvdp(2,:)

                dhdp(1,:) = -friction(1,1)*dvdp(1,:) - friction(1,2)*dvdp(2,:)
                dhdp(2,:) = -friction(2,1)*dvdp(1,:) - friction(2,2)*dvdp(2,:)

                dhdp(:,1) = dhdp(:,1) - invMetGrad(1,:)*momentum(1) - invMetGrad(2,:)*momentum(2)
                dhdp(:,2) = dhdp(:,2) - invMetGrad(2,:)*momentum(1) - invMetGrad(3,:)*momentum(2)

                ! write(*,*) 'dhdp',dhdp

                gGamma2 = matmul(randForce,Gamma2)
                ! write(*,*) 'gGamma2',gGamma2
                ! gGamma2(1) = randForce(1)*Gamma2(1) + randForce(2)*Gamma2(2)
                ! gGamma2(2) = randForce(2)*Gamma2(1) + randForce(3)*Gamma2(2)

                !actually updating
                momentum = momentum + dt*h + 0.5d0*dt**2 * (matmul(dhdq,v) + matmul(dhdp,h)) + matmul(dhdp,gGamma2) + &
                                matmul(randForce,Gamma1)
                ! momentum(1) = momentum(1) + randForce(1)*Gamma1(1) + randForce(2)*Gamma1(2)! + &
                !                 ! randForceGrad(1,1)*v(1)*Gamma3(1) + randForceGrad(2,1)*v(1)*Gamma3(2) + &
                !                 ! randForceGrad(1,2)*v(2)*Gamma3(1) + randForceGrad(2,2)*v(2)*Gamma3(2)
                ! momentum(2) = momentum(2) + randForce(2)*Gamma1(1) + randForce(3)*Gamma1(2)! + &
                !                 ! randForceGrad(2,1)*v(1)*Gamma3(1) + randForceGrad(3,1)*v(1)*Gamma3(2) + &
                !                 ! randForceGrad(2,2)*v(2)*Gamma3(1) + randForceGrad(3,2)*v(2)*Gamma3(2)

                coords = coords + dt*v + 0.5d0*dt**2.d0*(matmul(dvdq,v) + matmul(dvdp,h)) + matmul(dvdp,gGamma2)

                ! write(*,*) coords, momentum
                ! stop

                lastIter = dtIter

            end do

        end subroutine time_evolve

        subroutine run(initialCoords,savePaths,allCoords,allMomenta,lastIter,finishedSuccessfully)
            implicit none

            integer, dimension(:), allocatable, intent(out)                     :: lastIter
            logical, intent(in)                                                 :: savePaths
            logical, dimension(:), allocatable, intent(out)                     :: finishedSuccessfully
            real(kind=real64), dimension(:,:), intent(in)                       :: initialCoords
            real(kind=real64), dimension(:,:,:), allocatable, intent(out)       :: allCoords, allMomenta
            real(kind=real64), dimension(:,:), allocatable                      :: coordsOneRun, momentaOneRun

            integer, dimension(2)                                               :: shp
            integer                                                             :: nPaths, pathIter

            shp = shape(initialCoords)
            nPaths = shp(1)

            allocate(coordsOneRun(maxIters,2))
            allocate(momentaOneRun(maxIters,2))
            allocate(lastIter(nPaths))
            allocate(finishedSuccessfully(nPaths))

            if (savePaths) then
                allocate(allCoords(nPaths,maxIters,2))
                allocate(allMomenta(nPaths,maxIters,2))
            else
                !only saves the final value
                allocate(allCoords(nPaths,1,2))
                allocate(allMomenta(nPaths,1,2))
            end if

            do pathIter=1,nPaths
                call time_evolve(initialCoords(pathIter,:),coordsOneRun,momentaOneRun,&
                                 lastIter(pathIter),finishedSuccessfully(pathIter))
                write(*,*) 'Final point: ',coordsOneRun(lastIter(pathIter),:)
                if (savePaths) then
                    allCoords(pathIter,:,:) = coordsOneRun
                    allMomenta(pathIter,:,:) = momentaOneRun
                else
                    allCoords(pathIter,1,:) = coordsOneRun(lastIter(pathIter),:)
                    allMomenta(pathIter,1,:) = momentaOneRun(lastIter(pathIter),:)
                end if
            end do

        end subroutine run

end module langevin_2d
