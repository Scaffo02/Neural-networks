import numpy as np
from scipy.linalg import logm, expm
import matplotlib.pyplot as plt
from scipy.io import savemat

# ---------------------------------------------------------
# 1. Definitions of Basis States (from SQ.m)
# ---------------------------------------------------------

# Computational Basis (used for Q and T)
zero = np.array([[1], [0]])
uno = np.array([[0], [1]])

# Base R (using a = pi/5 as in SQ.m)
a = np.pi / 5
R_p = 1/np.sqrt(2) * (zero + 1j * np.exp(1j * a) * uno)
R_m = 1/np.sqrt(2) * (zero - 1j * np.exp(1j * a) * uno)

# Base S
S_p = 1/np.sqrt(2) * (zero + uno)
S_m = 1/np.sqrt(2) * (zero - uno)

# Kron products for projectors
# We will construct the correlation calculation logic in the next steps.

# ---------------------------------------------------------
# 2. Helper Functions
# ---------------------------------------------------------

def entropy(rho):
    """Calculates Von Neumann entropy of a density matrix."""
    # Ensure hermitian
    if np.allclose(rho, rho.conj().T):
        evals = np.linalg.eigvalsh(rho)
    else:
        # Should be hermitian but maybe numerical noise
        evals = np.linalg.eigvals(rho)
        evals = np.real(evals) # Take real part if complex due to noise

    # Remove negative or zero eigenvalues for log
    evals = evals[evals > 1e-12]
    # Normalize if needed (should be trace 1)
    # evals = evals / np.sum(evals)

    if len(evals) == 0:
        return 0.0

    return -np.sum(evals * np.log(evals))

def get_projector(v1, v2):
    """Returns projector |v1><v1| kron |v2><v2|"""
    # v1, v2 are column vectors
    # kron(v1, v2) is column vector
    psi = np.kron(v1, v2)
    return psi @ psi.conj().T

def calculate_diagonal_rho(rho, basis1, basis2):
    """
    Calculates the diagonal density matrix (post-measurement state)
    in the tensor product basis of basis1 and basis2.
    basis1: list of vectors [b1_0, b1_1]
    basis2: list of vectors [b2_0, b2_1]
    """
    rho_d = np.zeros_like(rho, dtype=complex)

    for b1 in basis1:
        for b2 in basis2:
            P = get_projector(b1, b2)
            # Term: (P rho P) ?? No, formula in SQ.m is:
            # ((kron...)' * rho * kron...) .* kron... * kron...'
            # This means scalar = <psi|rho|psi>. Then scalar * |psi><psi|.
            # This is exactly P rho P if P is rank-1 projector.
            # But formula says: (psi' * rho * psi) is a scalar.

            # Let psi = kron(b1, b2)
            psi = np.kron(b1, b2)
            scalar = (psi.conj().T @ rho @ psi).item() # <psi|rho|psi>
            rho_d += scalar * (psi @ psi.conj().T)

    return rho_d

def calculate_correlations(rho):
    """Calculates C_QS, C_RS, C_RT, C_QT."""

    # Bases lists
    # Q and T use computational basis (zero, uno)
    Basis_Q = [zero, uno]
    Basis_T = [zero, uno]
    Basis_R = [R_p, R_m]
    Basis_S = [S_p, S_m]

    # Entropies
    S_rho = entropy(rho)

    # C_QS
    rho_d_QS = calculate_diagonal_rho(rho, Basis_Q, Basis_S)
    C_QS = entropy(rho_d_QS) - S_rho

    # C_RS
    rho_d_RS = calculate_diagonal_rho(rho, Basis_R, Basis_S)
    C_RS = entropy(rho_d_RS) - S_rho

    # C_RT
    rho_d_RT = calculate_diagonal_rho(rho, Basis_R, Basis_T)
    C_RT = entropy(rho_d_RT) - S_rho

    # C_QT
    rho_d_QT = calculate_diagonal_rho(rho, Basis_Q, Basis_T)
    C_QT = entropy(rho_d_QT) - S_rho

    return C_QS, C_RS, C_RT, C_QT

# ---------------------------------------------------------
# 3. Hamiltonians and Evolution
# ---------------------------------------------------------

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])
identity = np.eye(2)

# Hamiltonian 1: Heisenberg Interaction
# H = J (X.X + Y.Y + Z.Z)
J1 = 1.0
H1 = J1 * (np.kron(sigma_x, sigma_x) + np.kron(sigma_y, sigma_y) + np.kron(sigma_z, sigma_z))

# Hamiltonian 2: Ising Interaction (Z.Z) with transverse field
# H = J Z.Z + h (X.I + I.X)
J2 = 1.0
h2 = 0.5
H2 = J2 * np.kron(sigma_z, sigma_z) + h2 * (np.kron(sigma_x, identity) + np.kron(identity, sigma_x))

# Initial State: Separable state |00>
psi0 = np.kron(zero, zero)
rho0 = psi0 @ psi0.conj().T

# Or maybe an entangled state?
# psi0 = 1/sqrt(2) * (|01> - |10>) (Singlet)
# If we start with eigenstate of H1, nothing happens.
# Let's start with |00> which is not eigenstate of Singlet Hamiltonian?
# Singlet H1 eigs: -3, 1, 1, 1. |00> is in the triplet subspace (eigenvalue 1).
# So |00> is eigenstate of H1. It won't evolve.
# Let's verify H1 |00>.
# XX |00> = |11>
# YY |00> = -|11>
# ZZ |00> = |00>
# (XX+YY+ZZ)|00> = |00>. Yes.
# So I should use a different initial state or different H1.
# Let's use psi0 = |0> x (|0> + |1>)/sqrt(2) = (|00> + |01>)/sqrt(2)
psi0 = np.kron(zero, (zero + uno)/np.sqrt(2))
rho0 = psi0 @ psi0.conj().T


# Time steps
times = np.linspace(0, 10, 100)
results_H1 = {'C_QS': [], 'C_RS': [], 'C_RT': [], 'C_QT': []}
results_H2 = {'C_QS': [], 'C_RS': [], 'C_RT': [], 'C_QT': []}

for t in times:
    # Evolution H1
    U1 = expm(-1j * H1 * t)
    rho_t_1 = U1 @ rho0 @ U1.conj().T

    c1 = calculate_correlations(rho_t_1)
    results_H1['C_QS'].append(c1[0])
    results_H1['C_RS'].append(c1[1])
    results_H1['C_RT'].append(c1[2])
    results_H1['C_QT'].append(c1[3])

    # Evolution H2
    U2 = expm(-1j * H2 * t)
    rho_t_2 = U2 @ rho0 @ U2.conj().T

    c2 = calculate_correlations(rho_t_2)
    results_H2['C_QS'].append(c2[0])
    results_H2['C_RS'].append(c2[1])
    results_H2['C_RT'].append(c2[2])
    results_H2['C_QT'].append(c2[3])


# ---------------------------------------------------------
# 4. Save and Plot
# ---------------------------------------------------------

# Save to .mat
# We need to restructure dict for savemat to be nice in MATLAB if needed, but dict is fine.
savemat('correlazioni_evolution.mat', {
    'times': times,
    'H1_C_QS': results_H1['C_QS'],
    'H1_C_RS': results_H1['C_RS'],
    'H1_C_RT': results_H1['C_RT'],
    'H1_C_QT': results_H1['C_QT'],
    'H2_C_QS': results_H2['C_QS'],
    'H2_C_RS': results_H2['C_RS'],
    'H2_C_RT': results_H2['C_RT'],
    'H2_C_QT': results_H2['C_QT']
})

# Calculate Bell Parameter (S-like quantity used in class_bell.ipynb)
# sum = C_QS + C_RS + C_RT - C_QT
def calc_S(res):
    return np.array(res['C_QS']) + np.array(res['C_RS']) + np.array(res['C_RT']) - np.array(res['C_QT'])

S1 = calc_S(results_H1)
S2 = calc_S(results_H2)

# Plotting
plt.figure(figsize=(12, 10))

# Plot Correlations H1
plt.subplot(2, 2, 1)
plt.plot(times, results_H1['C_QS'], label='C_QS')
plt.plot(times, results_H1['C_RS'], label='C_RS')
plt.plot(times, results_H1['C_RT'], label='C_RT')
plt.plot(times, results_H1['C_QT'], label='C_QT')
plt.title('Correlations H1 (Heisenberg)')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Entropy Diff')

# Plot Correlations H2
plt.subplot(2, 2, 2)
plt.plot(times, results_H2['C_QS'], label='C_QS')
plt.plot(times, results_H2['C_RS'], label='C_RS')
plt.plot(times, results_H2['C_RT'], label='C_RT')
plt.plot(times, results_H2['C_QT'], label='C_QT')
plt.title('Correlations H2 (Ising + Transverse)')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Entropy Diff')

# Plot Bell Parameter S
plt.subplot(2, 1, 2)
plt.plot(times, S1, label='S H1 (Heisenberg)')
plt.plot(times, S2, label='S H2 (Ising)')
plt.title('Bell Parameter Sum (C_QS + C_RS + C_RT - C_QT)')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Sum')
plt.grid(True)

plt.tight_layout()
plt.savefig('confronto_correlazioni_plot.png')
print("Simulation complete. Results saved to 'correlazioni_evolution.mat' and 'confronto_correlazioni_plot.png'.")
