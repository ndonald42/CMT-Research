import numpy as np
from scipy import sparse
from scipy.linalg import eig
from itertools import permutations, combinations
import numpy.linalg as la

# Helper function for Kronecker product
def fold_kronecker(matrices):
    result = matrices[0]
    for m in matrices[1:]:
        result = np.kron(result, m)
    return result

# Combinatorial function for binomial coefficient
def comb(n, k):
    from math import factorial
    return factorial(n) // (factorial(k) * factorial(n - k))

# Initialize lattice parameters and couplings
Latticedimensions = [2, 2]  # 2x2 lattice
partition = [[1, 2], [3, 4], [5, 6]]  # Lattice partition
outside = [7, 8]  # Outside sites
particles = []
Signum1 = []
J = 1
K = 1
h = 1

# Build Signum1 for entropy sign calculation
for k in range(1, len(partition) + 1):
    Signum1.append(np.full(comb(len(partition), k), (-1)**(k + 1)))
Signum2 = np.concatenate(Signum1)

# Define chains and set used to determine sign on v.n entropy
for i in range(Latticedimensions[1]):
    particles.append(2 * Latticedimensions[0])
L1 = len(particles)
L2 = sum(particles)
L3 = sum(particles[:-1])

# Define Pauli matricies and define magnetic field strength/direction
vec1 = (h / np.sqrt(3)) * np.array([1, 1, 1])
sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])
identity = np.eye(2)
sigma = [sigma_x, sigma_y, sigma_z]
U = vec1[0]*sigma[0] + vec1[1]*sigma[1] + vec1[2]*sigma[2]

# Define g function which finds bond type in Kitaev chain
def g(x, k):
    if x == 2 and k % 2 == 1:
        return sigma_x
    elif x == 2 and k % 2 == 0:
        return sigma_y
    elif x == 1:
        return identity
    return identity

# Function to build chain Kronecker product list for one chain
def chain(n):
    Permutation1 = list(permutations([2] + [1] * (n - 2) + [2]))
    Info = []
    for k in range(len(Permutation1)):
        row = []
        for j in range(1, n - 1):
            segment = [Permutation1[k][i-1:i] for i in range(j, j + 3)]
            row.append(int(segment == [[1], [2], [1]]))
        Info.append(row)
    Info = np.array(Info).T
    Positions = [i for i, x in enumerate(np.any(Info == 1, axis=0)) if x]
    Permutation1 = [Permutation1[i] for i in range(len(Permutation1)) if i not in Positions]
    result = [[g(Permutation1[y][l], y + 1) for l in range(n)] for y in range(len(Permutation1))]
    return result

# Helper for site indexing. Finds site number before and after site j
Helper = [[sum(particles[:j-1]), sum(particles[j:])] for j in range(1, L1 + 1)]

# Build chains and lists the kronecker products for each chain with appropriate length. Also joins the chains together
Chains = [chain(p) for p in particles]
Lattice1 = []
for i in range(L1):
    chain_i = Chains[i]
    before = Helper[i][0]
    after = Helper[i][1]
    lattice_chain = []
    for config in chain_i:
        lattice_chain.append([identity] * before + config + [identity] * after)
    Lattice1.append(lattice_chain)
Lattice2 = [item for sublist in Lattice1 for item in sublist]

# Z-bond interactions
Zbonds0 = [[(i, j) for j in range(1, particles[i-1] + 1)] for i in range(1, L1 + 1)]
Zbonds1 = []
for i in range(1, L1):
    row = []
    for j in range(1, particles[i-1] + 1):
        if (i+1, j-1) in Zbonds0[i] and Zbonds0[i-1][j-1][1] % 2 == 0:
            row.append([Zbonds0[i-1][j-1], Zbonds0[i][j-2]])
    Zbonds1.append([x for x in row if x])
Zbonds2 = []
for i in range(len(Zbonds1)):
    for j in range(len(Zbonds1[i])):
        site1 = Zbonds1[i][j][0]
        site2 = Zbonds1[i][j][1]
        idx1 = sum(particles[:site1[0]-1]) + site1[1]
        idx2 = sum(particles[:site2[0]-1]) + site2[1]
        Zbonds2.append([idx1, idx2])
Zbonds3 = []
for i in range(len(Zbonds2)):
    config = [identity] * L2
    config[Zbonds2[i][0] - 1] = sigma_z
    config[Zbonds2[i][1] - 1] = sigma_z
    Zbonds3.append(config)
Zbonds4 = Zbonds3

# Twisted boundary conditions
TBC1 = [[i, L3 + i + 1] for i in range(1, particles[0] + 1) if i % 2 == 1]
TBC2 = []
for i in range(len(TBC1)):
    config = [identity] * L2
    config[TBC1[i][0] - 1] = sigma_z
    config[TBC1[i][1] - 1] = sigma_z
    TBC2.append(config)

# Combine Kronecker products
KroneckerProducts1 = Lattice2 + Zbonds4 + TBC2

# Construct Hamiltonian H1 for the lattice
H1 = sparse.csr_matrix((2**L2, 2**L2), dtype=complex)
for kp in KroneckerProducts1:
    H1 += (J / 2) * sparse.csr_matrix(fold_kronecker(kp))

# Construct Hamiltonian H2 for the magnetic field interaction
KroneckerProducts2 = []
for i in range(L2):
    config = [identity] * L2
    config[i] = U
    KroneckerProducts2.append(config)
H2 = sparse.csr_matrix((2**L2, 2**L2), dtype=complex)
for kp in KroneckerProducts2:
    H2 += h * sparse.csr_matrix(fold_kronecker(kp))

# Total Hamiltonian
H = (H1 + H2).toarray()

# Compute ground state numerically
eigenvalues, eigenvectors = eig(-H)
lambda1 = -np.max(np.real(eigenvalues))
groundstate = eigenvectors[:, np.argmax(np.real(eigenvalues))]

# Puts kronecker product basis in binary
HamiltonianBasis = [np.binary_repr(j, width=L2) for j in range(2**L2)]
HamiltonianBasis = [[int(d) for d in x] for x in HamiltonianBasis]

# Compute all viable subsets of the partition
subsets = list(combinations(partition, r) for r in range(1, len(partition) + 1))
subsets = [item for sublist in subsets for item in sublist]
regions = [np.concatenate(s) for s in subsets]

# Lists all regions in lattice for the selected partition
regionbases1 = []
for r in regions:
    basis = [np.binary_repr(j, width=len(r)) for j in range(2**len(r))]
    regionbases1.append([[int(d) for d in x] for x in basis])

# Lists indices of hamiltonian basis vectors needed for each region
projectedspaces = []
for k in range(len(regionbases1)):
    region_k = regions[k]
    basis_k = regionbases1[k]
    space_k = []
    for l in range(len(basis_k)):
        indices = []
        for i in range(2**L2):
            matches = all(HamiltonianBasis[i][r-1] == basis_k[l][regions[k].tolist().index(r)] for r in region_k)
            if matches:
                indices.append(i + 1)
        space_k.append(indices)
    projectedspaces.append(space_k)

# Decomposes ground state into direct product of subspaces
projections = []
for i in range(len(projectedspaces)):
    proj_i = []
    for j in range(len(projectedspaces[i])):
        proj_ij = [groundstate[k-1] for k in projectedspaces[i][j] if k > 0]
        proj_i.append(proj_ij)
    projections.append(proj_i)

# Forms matrix representation of inner product bilinear form on each decomposition
BilinearForm = []
for i in range(len(projections)):
    form_i = []
    for j in range(len(projections[i])):
        row = []
        for k in range(len(projections[i])):
            row.append(np.dot(np.conj(projections[i][j]), projections[i][k]))
        form_i.append(row)
    BilinearForm.append(np.array(form_i))

# Diagonalize bilinear forms
Topo1 = [np.real(la.eigvals(bf)) for bf in BilinearForm]

# Partial von Neumann entropies
Topo2 = []
for i in range(len(Topo1)):
    entropy = -sum(x * np.log(x) for x in Topo1[i] if x > 1e-10)
    Topo2.append(np.real(entropy))
print("The Partial Von Neumann Entropies are:")
print(np.array(Topo2))

# Topological entanglement entropy
TopologicalEntropy = sum(Topo2[i] * Signum2[i] for i in range(len(Topo2)))
print("The Topological Entanglement Entropy is:")
print(TopologicalEntropy)
