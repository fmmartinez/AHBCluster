module definitions
implicit none

type Atom
   character(4) :: symbol
   real(8) :: mass,charge,ljSigma,ljEpsilon
   !mass in uma, charge in e, ljSigma in A, ljEpsilon in Kcal/mol
   real(8),dimension(1:3) :: pos,vel
   !positions in A, velocities in A/ps 
end type Atom

type AtomPairData
   real(8) :: rij,qq,ljEps,ljSig
   real(8),dimension(1:3) :: vectorij
end type AtomPairData

type ForceInAtom
   real(8),dimension(1:3) :: total
end type ForceInAtom

type Forces
   type(ForceInAtom),dimension(:),allocatable :: inAtom
   real(8),dimension(:,:),allocatable :: atomPair
end type Forces

real(8),parameter :: a = 11.2d0, b = 7.1d13, d = 110d0
real(8),parameter :: na = 9.26d0, da = 0.95d0
real(8),parameter :: nb = 11.42d0, db = 0.97d0
real(8),parameter :: c = 0.776d0

real(8),parameter :: forceToVelUnits = 418.4d0, kCoulomb = 332.06d0

end module definitions
