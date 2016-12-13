program clustermd
implicit none

type Atom
   character(4) :: symbol
   real(8) :: mass,charge,ljSigma,ljEpsilon
   !mass in uma, charge in e, ljSigma in A, ljEpsilon in Kcal/mol
   real(8),dimension(1:3) :: pos,vel
   !positions in A, velocities in A/ps 
end type Atom

type ForceInAtom
   real(8),dimension(1:3) :: total
end type ForceInAtom

real(8),parameter :: a = 11.2d0, b = 7.1d13, d = 110d0
real(8),parameter :: na = 9.26d0, da = 0.95d0
real(8),parameter :: nb = 11.42d0, db = 0.97d0
real(8),parameter :: c = 0.776d0
real(8),parameter :: forceToVelUnits = 418.4d0, kCoulomb = 332.06d0

integer :: i,j,unit1

real(8) :: f1,f2,f3,f4,fr,rah,rab,rbh,totalEnergy
real(8) :: timeStep,halfTimeStep
real(8),dimension(1:3) :: uvah,uvab,uvhb

type(Atom),dimension(1:3) :: abcomplex
type(ForceInAtom),dimension(1:3) :: force

timeStep = 0.0015d0
halfTimeStep = timeStep/2d0

!initialization
abcomplex(1)%symbol = 'O'
abcomplex(1)%mass = 93d0
abcomplex(1)%charge = 0d0
abcomplex(1)%ljSigma = 3.5d0
abcomplex(1)%ljEpsilon = 0.397d0
abcomplex(1)%pos = [1d0,1d0,1d0]
abcomplex(1)%vel = [0d0,0d0,0d0]

abcomplex(2)%symbol = 'N'
abcomplex(2)%mass = 59d0
abcomplex(2)%charge = 0d0
abcomplex(2)%ljSigma = 3.5d0
abcomplex(2)%ljEpsilon = 0.397d0
abcomplex(2)%pos = [3.7d0,1d0,1d0]
abcomplex(2)%vel = [0d0,0d0,0d0]

abcomplex(3)%symbol = 'H'
abcomplex(3)%mass = 1d0
abcomplex(3)%charge = 0.5d0
abcomplex(3)%ljSigma = 0d0
abcomplex(3)%ljEpsilon = 0d0
abcomplex(3)%pos = [2d0,1d0,1d0]
abcomplex(3)%vel = [5.0d0,0d0,0d0]

force(1)%total = [0d0,0d0,0d0]
force(2)%total = [0d0,0d0,0d0]
force(3)%total = [0d0,0d0,0d0]

!get force
rab = sqrt(sum((abcomplex(2)%pos - abcomplex(1)%pos)**2))
rah = sqrt(sum((abcomplex(3)%pos - abcomplex(1)%pos)**2))
rbh = sqrt(sum((abcomplex(3)%pos - abcomplex(2)%pos)**2))
!update charges
fr = 0.5d0*(1d0 + (rah - 1.43d0)/sqrt((rah-1.43d0)**2+0.125d0**2))
abcomplex(1)%charge = (1d0-fr)*(-0.5d0)+fr*(-1d0)
abcomplex(2)%charge = fr*(0.5d0)
!covalent forces in complex
f1 = (a*b)*exp(-a*rab)
f2 = c*d*exp(-nb*(rbh-db)**2/(2d0*rbh))*(nb*(rbh-db)/rbh)*((rbh-db)/(2d0*rbh)-1d0)
f3 =   d*exp(-na*(rah-da)**2/(2d0*rah))*(na*(rah-da)/rah)*((rah-da)/(2d0*rah)-1d0)
!electrostatic forces in complex
f4 = kCoulomb*abcomplex(1)%charge*abcomplex(2)%charge/rab**2
!vectorization
uvab = (abcomplex(2)%pos - abcomplex(1)%pos)/rab
uvah = (abcomplex(3)%pos - abcomplex(1)%pos)/rah
uvhb = (abcomplex(2)%pos - abcomplex(3)%pos)/rbh
force(1)%total = -f1*uvab - f3*uvah
force(2)%total =  f1*uvab + f2*uvhb
force(3)%total = -f2*uvhb + f3*uvah

!rab = 2.7d0
!do i = 1, 31
!   rah = rab*(i-1)/30d0
!  
!   f1 = (a*b)*exp(-a*rab)
!   f2 = c*d*exp(-nb*(rab-rah-db)**2/(2d0*(rab-rah)))*(nb*(rab-rah-db)/(rab-rah))*((rab-rah-db)/(2d0*(rab-rah))-1d0)
!   f3 = d*exp(-na*(rah-da)**2/(2d0*rah))*(na*(rah-da)/rah)*((rah-da)/(2d0*rah)-1d0)
!  
!   fvec = (abcomplex(2)%pos - abcomplex(1)%pos)/rab
!   force(3)%total = ( f2 + f3)*fvec
!   
!   write (222,'(4f26.8)') rah,force(3)%total
!   
!   totalEnergy = b*exp(-a*rab)+d*(1-exp(-(na*(rah-da)**2)/(2d0*rah)))+c*d*(1-exp(-nb*(rab-rah-db)**2/(2d0*(rab-rah))))
!  
!   write (333,'(5f26.8)') rah,f1,f2,f3,totalEnergy
!end do
!
!rah = 1d0
!do i = 1, 31
!   rab = 2.70d0*(i-1)/30d0 + 1d0
!  
!   f1 = (a*b)*exp(-a*rab)
!   f2 = c*d*exp(-nb*(rab-rah-db)**2/(2d0*(rab-rah)))*(nb*(rab-rah-db)/(rab-rah))*((rab-rah-db)/(2d0*(rab-rah))-1d0)
!   f3 = d*exp(-na*(rah-da)**2/(2d0*rah))*(na*(rah-da)/rah)*((rah-da)/(2d0*rah)-1d0)
!  
!   fvec = (abcomplex(2)%pos - abcomplex(1)%pos)/rab
!   force(3)%total = ( f2 + f3)*fvec
!   
!   write (222,'(4f26.8)') rah,force(3)%total
!   
!   totalEnergy = b*exp(-a*rab)+d*(1-exp(-(na*(rah-da)**2)/(2d0*rah)))+c*d*(1-exp(-nb*(rab-rah-db)**2/(2d0*(rab-rah))))
!  
!   write (444,'(5f26.8)') rab,f1,f2,f3,totalEnergy
!end do
!stop

do i = 1, 100000
   abcomplex(1)%vel = abcomplex(1)%vel + forceToVelUnits*halfTimeStep*force(1)%total/abcomplex(1)%mass
   abcomplex(2)%vel = abcomplex(2)%vel + forceToVelUnits*halfTimeStep*force(2)%total/abcomplex(2)%mass
   abcomplex(3)%vel = abcomplex(3)%vel + forceToVelUnits*halfTimeStep*force(3)%total/abcomplex(3)%mass

   abcomplex(1)%pos = abcomplex(1)%pos + timeStep*abcomplex(1)%vel
   abcomplex(2)%pos = abcomplex(2)%pos + timeStep*abcomplex(2)%vel
   abcomplex(3)%pos = abcomplex(3)%pos + timeStep*abcomplex(3)%vel

   !get force
   rab = sqrt(sum((abcomplex(2)%pos - abcomplex(1)%pos)**2))
   rah = sqrt(sum((abcomplex(3)%pos - abcomplex(1)%pos)**2))
   rbh = sqrt(sum((abcomplex(3)%pos - abcomplex(2)%pos)**2))
   !update charges
   fr = 0.5d0*(1d0 + (rah - 1.43d0)/sqrt((rah-1.43d0)**2+0.125d0**2))
   abcomplex(1)%charge = (1d0-fr)*(-0.5d0)+fr*(-1d0)
   abcomplex(2)%charge = fr*(0.5d0)
   !covalent forces in complex
   f1 = (a*b)*exp(-a*rab)
   f2 = c*d*exp(-nb*(rbh-db)**2/(2d0*rbh))*(nb*(rbh-db)/rbh)*((rbh-db)/(2d0*rbh)-1d0)
   f3 =   d*exp(-na*(rah-da)**2/(2d0*rah))*(na*(rah-da)/rah)*((rah-da)/(2d0*rah)-1d0)
   !electrostatic forces in complex
   f4 = kCoulomb*abcomplex(1)%charge*abcomplex(2)%charge/rab**2
   !vectorization
   uvab = (abcomplex(2)%pos - abcomplex(1)%pos)/rab
   uvah = (abcomplex(3)%pos - abcomplex(1)%pos)/rah
   uvhb = (abcomplex(2)%pos - abcomplex(3)%pos)/rbh
   force(1)%total = -f1*uvab - f3*uvah
   force(2)%total =  f1*uvab + f2*uvhb
   force(3)%total = -f2*uvhb + f3*uvah

   abcomplex(1)%vel = abcomplex(1)%vel + forceToVelUnits*halfTimeStep*force(1)%total/abcomplex(1)%mass
   abcomplex(2)%vel = abcomplex(2)%vel + forceToVelUnits*halfTimeStep*force(2)%total/abcomplex(2)%mass
   abcomplex(3)%vel = abcomplex(3)%vel + forceToVelUnits*halfTimeStep*force(3)%total/abcomplex(3)%mass

   if (mod(i,100) == 0) then
      open(newunit=unit1,file='trajectory.xyz',position='append')
      write(unit1,*) 3
      write(unit1,*) 'symbol - positions x y z'
      do j = 1, 3
         write(unit1,'(a4,3f14.8)') abcomplex(j)%symbol, abcomplex(j)%pos(1:3)
      end do
      close(unit1)
      write(222,'(i6,3f24.10)') i, abcomplex(3)%vel
      write(333,'(i6,3f24.10)') i, force(3)%total
   end if

   totalEnergy = b*exp(-a*rab)+d*(1-exp(-(na*(rah-da)**2)/(2d0*rah)))+c*d*(1-exp(-nb*(rbh-db)**2/(2d0*rbh)))
   write(111,*) i, rab, totalEnergy
end do

end program clustermd
