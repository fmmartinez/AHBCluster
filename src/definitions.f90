include 'mkl_vsl.f90'

module definitions
use mkl_vsl_type
use mkl_vsl
implicit none

!<--integration parameters
integer,parameter :: nPointsGrid = 40
real(8),parameter :: lowerLimit = 0.1d0, upperLimit = 2.4d0
real(8),parameter :: binWidth = (upperLimit-lowerLimit)/nPointsGrid
real(8),parameter :: covMinWell = 1.0d0, ionMinWell = 1.6d0
real(8),parameter :: alpha0 = 7.735d0
real(8),parameter :: alphaCov = 9.26d0
real(8),parameter :: alphaIon = 11.42d0
!-->

integer,parameter :: brng = VSL_BRNG_MT2203
integer,parameter :: method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
integer,parameter :: maxRattleCycles = 1000
integer,parameter :: maxShakeCycles = 10000

real(8),parameter :: a = 11.2d0, b = 7.1d13, d = 110d0
real(8),parameter :: na = 9.26d0, da = 0.95d0
real(8),parameter :: nb = 11.42d0, db = 0.97d0
real(8),parameter :: c = 0.776d0

real(8),parameter :: pi = 3.141592653589793d0
real(8),parameter :: pisqrt = sqrt(pi)
real(8),parameter :: forceToVelUnits = 418.4d0, kCoulomb = 332.06d0
real(8),parameter :: kBoltzmann = 0.831446d0
   !Boltzmann constant in amu*(A/ps)**2*K
real(8),parameter :: KtoKcalMol = 0.00239d0
   !units of (Kcal/mol)*(ps**2/A**2)/amu
real(8),parameter :: hbar = 0.01518d0
   !hbar constant in (Kcal/mol)*ps
real(8),parameter :: toleranceConstraints = 1d-7
real(8),parameter :: constrainedRS = 1.781d0

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

type PointData
   real(8) :: rij
   real(8),dimension(1:3) :: vectorij
end type PointData

type ForceInAtom
   real(8),dimension(1:3) :: total
end type ForceInAtom

type Forces
   type(ForceInAtom),dimension(:),allocatable :: inAtom
   real(8),dimension(:,:),allocatable :: atomPair
end type Forces

type MdData
   integer :: seed, eqSteps, nTrajectories
   integer :: stepFreqEqSave, stepFreqOutTrajectory
   integer :: stepFreqOutLog
   integer :: maxEqTries, eqPhases, eqPhaseSteps
   integer :: stepFreqVelRescale
   integer :: stepFreqCoMRemoval
   integer :: prodSteps
   integer :: nBondConstraints
   integer :: singleMap, updateLambdasOntheFly
   real(8) :: timeStep, halfTimeStep
   real(8) :: initialEqTempInK, targetTempInK
end type MdData

type BasisFunction
   real(8),dimension(1:nPointsGrid+1) :: gridPointValue
end type BasisFunction

type EvalOnGridFunction
   real(8),dimension(1:nPointsGrid+1) :: gridPointValue
end type EvalOnGridFunction

type EvalOnGridHData
   type(PointData),dimension(1:nPointsGrid+1) :: gridPoint
end type EvalOnGridHData

type MatrixList
   real(8),dimension(:,:),allocatable :: mat
end type MatrixList

type VectorForMatrix
   real(8),dimension(1:3) :: vecij
end type VectorForMatrix

type QuantumStateData
   real(8) :: hTraceN
   real(8),dimension(:),allocatable :: eigenvalues,rm,pm
   real(8),dimension(:,:),allocatable :: phiKphi, phiVsphi, SMatrix
   real(8),dimension(:,:),allocatable :: pqAp,pqBp,pAHp,pBHp,prAHp
   real(8),dimension(:,:),allocatable :: lambda,mapFactor
   real(8),dimension(:,:),allocatable :: h,vas,vbs,vhs, hs
   type(BasisFunction),dimension(:),allocatable :: phi
   type(EvalOnGridHData),dimension(:),allocatable :: gridHSolvent
   type(MatrixList),dimension(:),allocatable :: pirp,pir3p,pcr3p
!   type(MatrixList),dimension(:,:),allocatable :: pir2p
end type QuantumStateData

end module definitions
