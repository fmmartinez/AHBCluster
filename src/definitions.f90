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
   real(8) :: rij
   real(8),dimension(1:3) :: vectorij
end type AtomPairData

type ForceInAtom
   real(8),dimension(1:3) :: total
end type ForceInAtom

real(8),parameter :: forceToVelUnits = 418.4d0, kCoulomb = 332.06d0

end module definitions