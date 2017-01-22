include 'mkl_vsl.f90'

module definitions
use mkl_vsl_type
use mkl_vsl
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

type MdData
   integer :: seed, eqSteps, stepFreqEqSave, stepFreqOutTrajectory
   integer :: maxEqTries, eqPhases, eqPhaseSteps
   integer :: stepFreqVelRescale
   integer :: stepFreqCoMRemoval
   integer :: prodSteps
   integer :: nBondConstraints
   real(8) :: timeStep, halfTimeStep
   real(8) :: initialEqTempInK, targetTempInK
end type MdData

integer,parameter :: brng = VSL_BRNG_MT2203
integer,parameter :: method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
integer,parameter :: maxRattleCycles = 1000
integer,parameter :: maxShakeCycles = 10000

real(8),parameter :: a = 11.2d0, b = 7.1d13, d = 110d0
real(8),parameter :: na = 9.26d0, da = 0.95d0
real(8),parameter :: nb = 11.42d0, db = 0.97d0
real(8),parameter :: c = 0.776d0

real(8),parameter :: forceToVelUnits = 418.4d0, kCoulomb = 332.06d0
real(8),parameter :: kBoltzmann = 0.831446d0
   !Boltzmann constant in amu*(A/ps)**2*K
real(8),parameter :: KtoKcalMol = 0.00239d0
   !units of (Kcal/mol)*(ps**2/A**2)/amu
real(8),parameter :: toleranceConstraints = 1d-7
real(8),parameter :: constrainedRS = 1.781d0
end module definitions
