"""
This script calculates the problem of a shock followed by a CJ deflagration, 
propagating in a tube with rear venting.

It implements thermal and chemical equilibrium calculated with temperature
dependent heat capacities.

See associated article (to be published in Combustion and Flame, 2023)
for the details on the notation, terminology, method and algorithms.

It is sectionned into 4 different part :

	USER PARAMETERS		Where the user can set the gas initial conditions
						(pressure, temperature, mixture), and the
						rear wall opening area.

	LIST OF FUNCTIONS	Where all the functions used in this script are
						defined.

	PRE-PROCESSING		Where variables derived from the initial state
						are defined.


	MAIN SCRIPT			Where the calculations to determine the flow
						properties are run.


	POST-PROCESSING		Where variables derived from the flow properties
						are calculated and printed.
	

NOTES:

	- This script runs on Python 2. Running it on Python 3 will not work as 
	  the some of the syntax is different.

	- Requires the Python thermo-chemical calculation library CANTERA to work.
	  CANTERA is called when importing library CJ3.
	
	- Requires an older version of the Shock and Detonation Toolbox (SDToolbox)
	  for Python 2 developped by Shephered et al. (2015) to work.
	  It is provided with this package and called in the library CJ3.

	- Library CJ3 contains the functions meant to compute the
	  CJ deflagration state.
	  The method is derived from the CJ detonation calculations
	  implemented in the SDToolbox.
	
	- The terminology used for the different regime of CJ deflagration
	  propagation in vented tubes is inherited from our early work on the topic.
	  Hence :
	  		Critical 		stands for	Perfectly Expanded Jet 	(PEJ)
			Supercritic		stands for	Underexpanded Jet 		(UEJ)
			Subcritic		stands for 	Subsonic Jet 			(SSJ)
	  It was left as is in order not to break the code 
	
Willstrong Rakotoarison
wrako07a@uottawa.ca
willstrong.rkt@gmail.com

Created on 		Jan. 2019
Last updated on	Nov. 2023

"""


from CJ3 import *
import numpy as np
import scipy.optimize as optim
import sys

# ==============================================================================
# ========================== USER PARAMETERS ===================================
# ==============================================================================

# Initial conditions : pressure, temperature, mixture
P1 = 101.325e3
T1 = 295
mix = "C3H8:1, O2:5, N2:%f" %(3.76*5)

# Chemcical mechanism to use
mech = "gri30.cti"

# Ratio of vent opening to tube area
A4A3 = 0.1


# ==============================================================================
# ========================= LIST OF FUNCTIONS ==================================
# ==============================================================================


# GENERAL FUNCTIONS ============================================================

def Post_Shock_CJdef(Ms):
	"""
		Returns properties of the flow after the CJ deflagration
		for a given incident shock Mach number
	"""
		
	# Shocked state
	Dshock = Ms * c1
	print "Dshock = %-12.5f  -->  Mshock = %-12.5f" %(Dshock, Ms)

	gas2 = ct.Solution(mech)
	gas2 = PostShock_fr(Dshock, P1, T1, mix, mech)
	
	rho2 = gas2.density
	P2 = gas2.P
	T2 = gas2.T

	gamma2 = gas2.cp / gas2.cv
	rgas2 = ct.gas_constant / gas2.mean_molecular_weight
	
	u2 = Dshock * rho1/rho2
	U2 = Dshock - u2	


	# CJ deflagration state
	[scj, gas3, m3] = CJ_deflagration(P2, T2, mix, mech, ERRFT = 1e-10, ERRFV = 1e-10)
	rho3 = gas3.density
	P3 = gas3.P
	T3 = gas3.T
	X3 = gas3.X
	
	c3, c3_fr = equilSoundSpeeds(gas3)
	gas3.TPX = T3, P3, X3

	Scj = U2 + scj
	u3 = scj * rho2/rho3
	U3 = Scj - u3
	M3 = U3/c3
	
	# Note : Take -U3 instead of U3 because U3 here is < 0, and U3 is used
	# only for isentropic calculations
	return gas3, -U3, -M3, Scj
	

def Print_Results():
    """
    This function is used to print the results at the end of calculations.
    """
    print("\n\n")

    print("=================================================")
    print("==================== RESULTS ====================")
    print("=================================================")

    print("\n================= FLOW CONDITIONS ===============")
    print("%-40s%-15f %s" % ("Initial pressure", P1, "Pa"))
    print("%-40s%-15f %s" % ("Initial temperature", T1, "K"))
    print("%-40s%-15s %s" % ("Initial temperature", mix, ""))
    print("%-40s%-15f %s" % ("Vent to tube area ratio A4/A3", A4A3, ""))

    print("")
    print("%-40s%-15f" % ("A4/A3 for a perfectly expanded jet is", 1.0 / rA_PEJ))
    print("%-40s%-15s %s" % ("The rear wall is a ", flow_regime, ""))

    print("\n================= SHOCKED STATE (2) =============")
    print("%-40s%-15f %s" % ("Shock speed", Ds, "m/s"))
    print("%-40s%-15f %s" % ("Shock Mach number", Ms, ""))
    print("%-40s%-15f %s" % ("Post shock pressure", gas2.P, "Pa"))
    print("%-40s%-15f %s" % ("Post shock temperature", gas2.T, "K"))

    print("\n============ CJ DEFLAGRATION STATE (3) ==========")
    print("%-40s%-15f %s" % ("CJ deflagration apparent speed", Scj, "m/s"))
    print("%-40s%-15f %s" % ("CJ deflagration relative speed", scj, "m/s"))
    print("%-40s%-15f %s" % ("CJ deflagration Mach number", MCJ, ""))
    print("%-40s%-15f %s" % ("Post CJ-deflagration pressure", gas3.P, "Pa"))
    print("%-40s%-15f %s" % ("Post CJ-deflagration temperature", gas3.T, "K"))

    print(" ")
    return

# OPEN AND SEMI-OPEN PROBLEMS ==================================================
def Open_tube(Ms):
    """
    Use Post_Shock_CJdef to determine the shock velocity to propagate a
    shock followed by a CJ deflagration in an open tube.
    P3 must be equal to P1 (to fit atmosphere condition).
    Pass this function to a root finder algorithm.
    """

    # Calculate state 3
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms)

    return gas3.P - P1

def SemiOpen_tube(Ms):
    """
    Use Shock_CJdeflagration to determine the shock velocity to propagate a
    shock followed by a CJ deflagration in a semi-open tube.
    U3 must be 0 (to fit end wall condition).
    Pass this function to a root finder algorithm.
    """

    # Calculate state 3
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms)

    return U3


# PERFECTLY EXPANDED JET (PEJ) =================================================


def P4s_minus_P1(Ms):
    """
    This function returns (P4s - P1), for a given incident shock Mach
    number.
    It must be equal to 0 to get a flow adapted to the exit conditions.
    """

    # Calculate state 3
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms)
    s3 = gas3.s
    h3 = gas3.h
    P3 = gas3.P

    # Create the gas object gas4 for state (4)
    gas4 = ct.Solution(mech)

    # Subfunction used to find P4s such that M4 = 1
    def find_P4s(P):

        # Set the state using (P, s3)
        gas4.SPX = s3, P, gas3.X
        gas4.equilibrate("SP")

        Ttemp = gas4.T
        Xtemp = gas4.X

        # Energy conservation
        U4 = (2.0 * (h3 - gas4.h + U3**2 / 2.0))**0.5

        # Sound speed
        c4, c4_fr = equilSoundSpeeds(gas4)
        gas4.TPX = Ttemp, P, Xtemp

        # Mach number
        M4 = U4 / c4

        return M4 - 1.0

    # Find P4s such that M4 = 1
    P4s = optim.newton(find_P4s, P3 / 2.0, tol=1e-5)

    return P4s - P1


def A_Astar(gas3, U3, P4s):
    """
    Returns A/A*, the ratio of sections between the channel and the
    vent opening.
    """

    # State (3) properties
    s3 = gas3.s
    h3 = gas3.h
    P3 = gas3.P

    # Gas object for state (4*)
    gas4s = ct.Solution(mech)

    # Set the state using (P4s, s3)
    gas4s.SPX = s3, P4s, gas3.X
    gas4s.equilibrate("SP")

    Ttemp = gas4s.T
    Xtemp = gas4s.X

    # Energy conservation
    U4s = (2.0 * (h3 - gas4s.h + U3**2 / 2.0))**0.5

    # Sound speed
    c4, c4_fr = equilSoundSpeeds(gas4s)
    gas4s.TPX = Ttemp, P4s, Xtemp

    # Mach number
    M4 = U4s / c4

    # Area ratio from mass conservation
    return (gas4s.density * U4s) / (gas3.density * U3)


def Find_PEJ_regime(Ms_guess=3):
    """
    Returns:
        - Ms_PEJ: minimum incident shock Mach number for PEJ regime
        - A3/A4*: area ratio for choked flow
    """

    # Find Ms for PEJ regime
    Ms_PEJ = optim.newton(P4s_minus_P1, Ms_guess, tol=1e-5)
    # Ms_PEJ = optim.brentq(P4s_minus_P1, Ms_open, Ms_semiopen, xtol=1e-5)

    # State (3) at PEJ regime
    gas3_choked, U3_choked, M3_choked, Scj_choked = Post_Shock_CJdef(Ms_PEJ)

    # Pressure P4s for choked flow
    P4smP1 = P4s_minus_P1(Ms_PEJ)
    P4s = P4smP1 + P1

    # Area ratio
    A3_A4s = A_Astar(gas3_choked, U3_choked, P4s)

    return Ms_PEJ, A3_A4s


def State_4_PEJ(Ms_PEJ):
    """
    Returns state (4) for the PEJ regime.
    """

    # State (3)
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms_PEJ)
    s3 = gas3.s
    h3 = gas3.h
    P3 = gas3.P

    # Pressure P4s
    P4smP1 = P4s_minus_P1(Ms_PEJ)
    P4s = P4smP1 + P1

    # Gas object for state (4)
    gas4 = ct.Solution(mech)

    # Set the state using (P4s, s3)
    gas4.SPX = s3, P4s, gas3.X
    gas4.equilibrate("SP")

    Ttemp = gas4.T
    Xtemp = gas4.X

    # Energy conservation
    U4s = (2.0 * (h3 - gas4.h + U3**2 / 2.0))**0.5

    # Sound speed
    c4, c4_fr = equilSoundSpeeds(gas4)
    gas4.TPX = Ttemp, P4s, Xtemp

    # Mach number
    M4 = U4s / c4

    return gas4, M4

# SUBSONIC JET (SSJ) ===========================================================

def Find_State_4_SSJ(Ms_SSJ):
    """
    This function returns (rA_iter - rA), for an arbitrary
    Ms_SSJ in the range expected for the SSJ regime.
    """

    # Calculate state 3
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms_SSJ)
    s3 = gas3.s
    h3 = gas3.h

    # Create the gas object gas4 for state (4)
    gas4 = ct.Solution(mech)

    # Set the state using (P4s, s3)
    gas4.SPX = s3, P4s, gas3.X
    gas4.equilibrate("SP")

    # dh = h3 + U3^2/2 - h4
    dh = h3 - gas4.h + U3**2 / 2.0

    # If dh < 0, impossible to compute U4
    if dh <= 0.0:
        return -1.0

    else:
        # Flow velocity at state (4)
        U4 = (2.0 * dh)**0.5

        # Area ratio to be iterated
        rA_iter = (gas4.density * U4) / (gas3.density * U3)

        # Return difference to target area ratio
        return rA_iter - rA


def State_4_SSJ(Ms_SSJ):
    """
    This function returns state (4) from the shock Mach number
    at the SSJ regime.
    """

    # Calculate state 3
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms_SSJ)
    s3 = gas3.s
    h3 = gas3.h

    # Create the gas object gas4 for state (4)
    gas4 = ct.Solution(mech)

    # Set the state using (P4s, s3)
    gas4.SPX = s3, P4s, gas3.X
    gas4.equilibrate("SP")

    Ttemp = gas4.T
    Xtemp = gas4.X

    # Flow velocity from energy conservation
    U4 = (2.0 * (h3 - gas4.h + U3**2 / 2.0))**0.5

    # Sound speed at state (4)
    c4, c4_fr = equilSoundSpeeds(gas4)
    gas4.TPX = Ttemp, P4s, Xtemp

    # Mach number at state (4)
    M4 = U4 / c4

    return gas4, M4

# UNDEREXPANDED JET (UEJ) ======================================================

def Find_State_4_UEJ(Ms_UEJ):
    """
    This function returns (rA_iter - rA), for an arbitrary Ms_UEJ
    in the range expected for the UEJ regime.
    """

    # Calculate state 3
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms_UEJ)
    s3 = gas3.s
    h3 = gas3.h
    P3 = gas3.P

    # Safety check on U3 and M3
    if U3 < 0:
        print(f"\n!!WARNING!! : U3 = {U3:e} --> Making U3 = M3 = 1e-10\n")
        U3 = 1e-10
        M3 = 1e-10

    # Create the gas object gas4 for state (4)
    gas4 = ct.Solution(mech)

    # Subfunction to find P4s such that M4 = 1
    def find_P4s(P):

        # Set the state (4) using (P, s3)
        gas4.SPX = s3, P, gas3.X
        gas4.equilibrate("SP")

        Ttemp = gas4.T
        Xtemp = gas4.X

        # Flow velocity from energy conservation
        U4 = (2.0 * (h3 - gas4.h + U3**2 / 2.0))**0.5

        # Sound speed
        c4, c4_fr = equilSoundSpeeds(gas4)
        gas4.TPX = Ttemp, P, Xtemp

        # Mach number
        M4 = U4 / c4

        return M4 - 1.0

    # Find P4s such that M4 = 1
    P4s = optim.newton(find_P4s, P3 / 1.5, maxiter=500, tol=1e-5)

    # Set the final state (4)
    gas4.SPX = s3, P4s, gas3.X
    gas4.equilibrate("SP")

    # Flow velocity
    U4 = (2.0 * (h3 - gas4.h + U3**2 / 2.0))**0.5

    # Area ratio
    rA_iter = (gas4.density * U4) / (gas3.density * U3)

    return rA_iter - rA


def State_4_UEJ(Ms_UEJ):
    """
    This function returns state (4) from the shock Mach number
    at the UEJ regime, Ms_UEJ.
    """

    # Calculate state 3
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms_UEJ)
    s3 = gas3.s
    h3 = gas3.h
    P3 = gas3.P

    # Create the gas object gas4 for state (4)
    gas4 = ct.Solution(mech)

    # Subfunction to find P4s such that M4 = 1
    def find_P4s(P):

        gas4.SPX = s3, P, gas3.X
        gas4.equilibrate("SP")

        Ttemp = gas4.T
        Xtemp = gas4.X

        # Flow velocity
        U4 = (2.0 * (h3 - gas4.h + U3**2 / 2.0))**0.5

        # Sound speed
        c4, c4_fr = equilSoundSpeeds(gas4)
        gas4.TPX = Ttemp, P, Xtemp

        # Mach number
        M4 = U4 / c4

        return M4 - 1.0

    # Find P4s such that M4 = 1
    P4s = optim.newton(find_P4s, P3 / 2.0, maxiter=500, tol=1e-5)

    # Set the final state (4)
    gas4.SPX = s3, P4s, gas3.X
    gas4.equilibrate("SP")

    Ttemp = gas4.T
    Xtemp = gas4.X

    # Flow velocity
    U4 = (2.0 * (h3 - gas4.h + U3**2 / 2.0))**0.5

    # Sound speed
    c4, c4_fr = equilSoundSpeeds(gas4)
    gas4.TPX = Ttemp, P4s, Xtemp

    # Mach number
    M4 = U4 / c4

    return gas4, M4

# ==============================================================================
# ============================ PRE-PROCESSING ==================================
# ==============================================================================

# Create gas object
gas1 = ct.Solution(mech)
gas1.TPX = T1, P1, mix

# Initial state properties
rho1 = gas1.density

gamma1 = gas1.cp / gas1.cv
rgas1 = ct.gas_constant / gas1.mean_molecular_weight
c1 = (gamma1 * rgas1 * T1)**0.5


# ==============================================================================
# ============================= MAIN SCRIPT ====================================
# ==============================================================================


# Get the PEJ regime to get the choked flow at state (4) =======================

print("\nCalculating perfectly expanded jet solution")
print("Iterating on shock speed, currently calculating for :")
Ms_PEJ, rA_PEJ = Find_PEJ_regime()

# Check if rA_PEJ is conform
if 1.0 / rA_PEJ > 1.0 or 1.0 / rA_PEJ < 0.0:
    print("\niPEJ regime not conform (< 0 or > 1) :")
    print("A4/A3_PEJ = ", 1.0 / rA_PEJ)
    print("Calculations stoped")
    sys.exit()

br_PEJ = 1.0 - 1.0 / rA_PEJ
A4A3_PEJ = 1.0 / rA_PEJ

gas4, M4 = State_4_PEJ(Ms_PEJ)
P4s = gas4.P

Ds_PEJ = Ms_PEJ * c1
gas3, U3, M3, Scj_PEJ = Post_Shock_CJdef(Ms_PEJ)

gas2_PEJ = ct.Solution(mech)
gas2_PEJ = PostShock_fr(Ds_PEJ, P1, T1, mix, mech)

P2_PEJ = gas2_PEJ.P
overP_PEJ = P2_PEJ - P1


# Calculate Semi-open tube problem =============================================

print("\nCalculating Semi-open (closed-end) tube problem (Chue, 1993)")
print("Iterating on shock speed, currently calculating for :")
Ms_semiopen = optim.newton(SemiOpen_tube, 2, maxiter=500, tol=1e-5)

Ds_semiopen = Ms_semiopen * c1
gas3, U3, M3, Scj_semiopen = Post_Shock_CJdef(Ms_semiopen)

print("U3 = ", U3)

gas2_semiopen = ct.Solution(mech)
gas2_semiopen = PostShock_fr(Ds_semiopen, P1, T1, mix, mech)

P2_semiopen = gas2_semiopen.P
overP_semiopen = P2_semiopen - P1

# Output results directly if computing for closed end tube
if A4A3 == 0.0:

    print("\nDirectly outputing closed end tube solution")
    flow_regime = "Closed end (Chue, 1993)"

    Ds = Ds_semiopen
    Ms = Ms_semiopen
    gas2 = gas2_semiopen
    P2 = P2_semiopen

    overP = P2 - P1

    # Find Scj
    gas3, U3, M3, Scj = Post_Shock_CJdef(Ms)

    # Find state 4 from the SSJ shock Mach number
    gas4, M4 = State_4_SSJ(Ms)
    P4 = gas4.P

    # Get CJ-deflagration state
    T2 = gas2.T
    scj, gas3, m3 = CJ_deflagration(
        P2, T2, mix, mech, ERRFT=1e-10, ERRFV=1e-10
    )

    # Get the CJ-deflagration Mach number
    gamma3 = gas3.cp / gas3.cv
    c3 = (gamma3 * gas3.P / gas3.density)**0.5
    MCJ = scj / c3

    # Print results
    Print_Results()

    sys.exit()

# Otherwise, set rA = 1/A4A3
else:
    rA = 1.0 / A4A3
    

##### FIN DE CONVERSION py2 -> py3 avec CHATITO #####


# Find subsonic jet (SSJ) regime : for rA < rA_PEJ =================================

if A4A3 > A4A3_PEJ:

	print "\nCalculating subsonic-jet solution"
	print "Iterating on shock speed, currently calculating for :"
		
	# Find the incident shock Mach number for the SSJ regime
	Ms = optim.newton(Find_State_4_SSJ, Ms_PEJ, tol = 1e-5)
	Ds = Ms * c1
	
	# Find incident shock overpressure	
	gas2 = ct.Solution(mech)
	gas2 = PostShock_fr(Ds, P1, T1, mix, mech)
	P2 = gas2.P
	overP = P2 - P1
	
	# Find Scj
	gas3, U3, M3, Scj = Post_Shock_CJdef(Ms)
	
	# Find state 4 from the SSJ shock Mach number
	gas4, M4 = State_4_SSJ(Ms)
	P4 = gas4.P
	
	flow_regime = "Subsonic jet (SSJ)"



# Find underexpanded jet (UEJ) regime : for rA > rA_PEJ ========================

if A4A3 <= A4A3_PEJ:

	print "\nCalculating underexpanded jet solution"
	print "Iterating on shock speed, currently calculating for :"
	
	# Find the incident shock Mach number for the UEJ regime
	Ms = optim.brentq(Find_State_4_UEJ, Ms_PEJ, Ms_semiopen, xtol = 1e-5)
	Ds = Ms * c1
	
	# Find incident shock overpressure
	gas2 = ct.Solution(mech)
	gas2 = PostShock_fr(Ds, P1, T1, mix, mech)

	P2 = gas2.P
	overP = P2 - P1	
	
	# Find Scj
	gas3, U3, M3, Scj = Post_Shock_CJdef(Ms)
	
	# Find state 4 from the super shock Mach number
	gas4, M4 = State_4_UEJ(Ms)

	flow_regime = "Underexpanded (sonic) jet (UEJ)"


# ==============================================================================
# =========================== POST-PROCESSING ==================================
# ==============================================================================

# Get CJ-deflagration state
T2 = gas2.T
[scj, gas3, m3] = CJ_deflagration(P2, T2, mix, mech, ERRFT = 1e-10, ERRFV = 1e-10)

# Get the CJ-deflagration Mach number
gamma3 = gas3.cp / gas3.cv
c3 = (gamma3 * gas3.P / gas3.density)**0.5
MCJ = scj/c3

Print_Results()















