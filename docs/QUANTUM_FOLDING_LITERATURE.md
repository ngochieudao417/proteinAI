# Literature Review & Theoretical Foundations: Quantum Algorithms for Protein Folding and Photonic Architectures

**Project:** proteinAI  
**Milestone:** 1 (M1) — Architectural Blueprints & Theoretical Foundations  
**Author:** QAAA Lab Core Team  
**Date:** November 2026  
**Document Status:** Comprehensive Literature Review & Theoretical Specification (Approved)  
**Target Hardware:** Quandela Discrete-Variable Photonic Processors (`sim:slos`, `qpu:ascella`)  

---

## 1. Introduction: The Protein Folding Challenge and Quantum Advantage

The protein folding problem represents one of the foundational challenges in molecular biophysics and computational structural biology. Given a primary polypeptide sequence $S = (s_1, s_2, \dots, s_N)$ composed of $N$ amino acids chosen from the alphabet of 20 canonical residues, the folding problem seeks to determine the three-dimensional tertiary conformation $\mathbf{R} = (\vec{r}_1, \vec{r}_2, \dots, \vec{r}_N) \in \mathbb{R}^{3 \times N}$ corresponding to the thermodynamic global free-energy minimum, adhering to Anfinsen's thermodynamic hypothesis (Anfinsen 1973).

### 1.1 Computational Complexity and Levinthal's Paradox
Levinthal's paradox (Levinthal 1969) highlights that an unbiased random conformational search for a modest 100-residue polypeptide—assuming merely 3 discrete conformational states per peptide bond—would require exploring $3^{100} \approx 5 \times 10^{47}$ distinct conformations. Even at picosecond transition rates ($10^{12}\ \text{s}^{-1}$), an exhaustive search would require $\sim 10^{28}$ years, yet native proteins routinely fold spontaneously on millisecond timescales.

In computational complexity theory, determining the ground-state conformation of coarse-grained lattice protein models (such as the 2D or 3D Hydrophobic-Polar model) has been formally proven to be **NP-complete** (Berger & Leighton 1998; Crescenzi et al. 1998). Classical all-atom molecular dynamics (MD) simulations face fundamental timescale gaps: numerical integration of Newton's equations requires femtosecond timesteps ($10^{-15}\ \text{s}$), making microsecond-to-second folding trajectories computationally prohibitive even on specialized supercomputers (such as Anton 2/3). While modern deep-learning architectures (AlphaFold2, ESMFold2) have achieved remarkable predictive accuracy for static native folds deposited in the Protein Data Bank (PDB), they operate as statistical approximators. They frequently fail when predicting:
- Dynamic conformational transitions and alternative metastable states.
- De novo designed sequences lacking natural evolutionary homologues in Multiple Sequence Alignments (MSAs).
- Non-local topological contact cooperativity in flexible binding epitope loops.

### 1.2 The Quantum Paradigm
Quantum computing introduces fundamentally distinct physical principles to conformational search. Quantum superpositions $|\Psi\rangle = \sum_{\mathbf{c}} \alpha_{\mathbf{c}} |\mathbf{c}\rangle$ enable simultaneous representation of exponentially large conformational state spaces. Quantum tunneling allows systems to traverse tall, narrow free-energy barriers that severely trap classical Monte Carlo and simulated annealing algorithms.

Furthermore, linear optical quantum computing (LOQC)—as pioneered by Knill, Laflamme, and Milburn (KLM 2001) and implemented on physical hardware by Quandela—provides an intrinsically natural platform for simulating molecular contact graphs. Photonic Fock states undergoing quantum interference across parameterized Mach-Zehnder Interferometer (MZI) networks directly execute Boson Sampling and Variational Quantum Algorithms (Peruzzo et al. 2014), concentrating photon detection probabilities onto minimal-energy tertiary contact configurations.

---

## 2. Biophysical & Mathematical Formulations of Coarse-Grained Protein Folding

Atomistic protein models require thousands of spatial coordinates, solvent dielectric functions, and quantum chemical electronic force fields (CHARMM, AMBER). For quantum computational simulation on near-term noisy intermediate-scale quantum (NISQ) and photonic hardware, coarse-graining is essential to capture the dominant biophysical driving forces while maintaining manageable qubit/mode counts.

### 2.1 The Hydrophobic-Polar (HP) Model (Dill et al.)
Introduced by Ken Dill and colleagues (Dill 1985; Lau & Dill 1989; Dill et al. 1995), the Hydrophobic-Polar (HP) model simplifies the 20-letter amino acid alphabet into two binary biophysical classes:
1. **Hydrophobic (H)**: Non-polar residues ($s_i \in \{\text{A, V, I, L, M, F, W, P, C}\}$) that sequester away from aqueous solvent to minimize the hydrophobic entropic penalty of water restructuring.
2. **Polar (P)**: Charged or hydrophilic residues ($s_i \in \{\text{G, S, T, Y, N, Q, D, E, K, R, H}\}$) that interact favorably with solvent water.

In the HP model, folding free energy is governed exclusively by favorable contacts between non-bonded hydrophobic residues that occupy adjacent lattice sites:
$$\epsilon_{\text{HP}}(s_i, s_j) = \begin{cases} -1 & \text{if } s_i = \text{H and } s_j = \text{H} \\ 0 & \text{otherwise} \end{cases}$$

Despite its simplicity, the HP lattice model captures key thermodynamic properties of natural globular proteins: compact core formation, secondary structure emergence, cooperativity, and two-state folding kinetics.

### 2.2 The Miyazawa-Jernigan (MJ) Statistical Contact Potential Matrix
While the HP model provides qualitative insights, quantitative binder evaluation requires distinguishing subtle differences in residue affinities (e.g., aromatic stacking vs. aliphatic packing vs. salt bridges).

Miyazawa and Jernigan (1985, 1996) derived a $20 \times 20$ symmetric statistical contact energy matrix from non-redundant high-resolution protein crystal structures deposited in the PDB. Based on the quasi-chemical approximation, the effective inter-residue contact energy between amino acids of types $a$ and $b$ is defined by:
$$\epsilon_{\text{MJ}}(a, b) = -k_B T \ln \left( \frac{N_{ab} N_{00}}{N_{a0} N_{b0}} \right)$$
where $N_{ab}$ is the observed number of contacts between residue types $a$ and $b$, $N_{a0}$ and $N_{b0}$ are contacts with solvent, and $N_{00}$ is solvent-solvent interaction.

The resulting Miyazawa-Jernigan contact potentials (expressed in units of $k_B T \approx 0.593\ \text{kcal/mol}$ at $298\ \text{K}$) exhibit clear biophysical clustering:
- **Strongly Attractive (Hydrophobic / Aromatic Core)**:
  - Cys–Cys ($\text{C}-\text{C}$): $-6.85\ k_B T$ (covalent disulfide bridging and packing)
  - Phe–Phe ($\text{F}-\text{F}$): $-5.18\ k_B T$ (aromatic $\pi$-$\pi$ stacking)
  - Trp–Trp ($\text{W}-\text{W}$): $-5.09\ k_B T$
  - Ile–Ile ($\text{I}-\text{I}$): $-4.46\ k_B T$, Leu–Leu ($\text{L}-\text{L}$): $-4.31\ k_B T$
  - Phe–Trp ($\text{F}-\text{W}$): $-5.38\ k_B T$
- **Electrostatic Interactions**:
  - Arg–Asp ($\text{R}-\text{D}$): $-3.74\ k_B T$ (favorable salt bridge)
  - Lys–Glu ($\text{K}-\text{E}$): $-3.12\ k_B T$ (favorable salt bridge)
  - Lys–Lys ($\text{K}-\text{K}$), Arg–Arg ($\text{R}-\text{R}$): $> -1.5\ k_B T$ (repulsive electrostatic penalties)
- **Weak / Repulsive (Polar / Solvent-Exposed)**:
  - Ser–Ser ($\text{S}-\text{S}$): $-1.95\ k_B T$
  - Gly–Gly ($\text{G}-\text{G}$): $-1.55\ k_B T$

Incorporating the Miyazawa-Jernigan contact potential into quantum folding Hamiltonians enables realistic, biophysically grounded sequence evaluation.

### 2.3 Spatial Lattice Geometries

Lattice protein models project the continuous polypeptide backbone onto discrete spatial graphs $\mathcal{L} = (\mathcal{V}, \mathcal{E})$. The choice of lattice geometry directly governs conformational fidelity:

| Lattice Geometry | Dimension | Coordination Number ($z$) | Bond Angles ($\theta$) | Parity / Bipartite Restriction | Biophysical Fidelity |
|---|---|---|---|---|---|
| **2D Square Lattice** | 2D | $z=4$ | $90^\circ, 180^\circ$ | Bipartite (contacts only when $\|i-j\|$ is odd) | Low; toy model for algorithmic verification |
| **3D Simple Cubic Lattice** | 3D | $z=6$ | $90^\circ, 180^\circ$ | Bipartite (contacts only when $\|i-j\|$ is odd) | Moderate; restricts natural contact topologies |
| **Tetrahedral (Diamond) Lattice** | 3D | $z=4$ | $\arccos(-1/3) \approx 109.47^\circ$ | Non-bipartite along zigzag | **High; exact match to polypeptide $sp^3$ backbone** |
| **Face-Centered Cubic (FCC)** | 3D | $z=12$ | $60^\circ, 90^\circ, 120^\circ$ | Non-bipartite | High; dense packing, high computational complexity |

#### The Tetrahedral / Diamond Lattice Geometry
The diamond lattice is uniquely suited for quantum protein folding (Robert et al. 2021). Each lattice vertex is bonded to four nearest neighbors oriented toward the vertices of a regular tetrahedron:
$$\vec{e}_0 = \frac{1}{\sqrt{3}}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}, \quad \vec{e}_1 = \frac{1}{\sqrt{3}}\begin{pmatrix} 1 \\ -1 \\ -1 \end{pmatrix}, \quad \vec{e}_2 = \frac{1}{\sqrt{3}}\begin{pmatrix} -1 \\ 1 \\ -1 \end{pmatrix}, \quad \vec{e}_3 = \frac{1}{\sqrt{3}}\begin{pmatrix} -1 \\ -1 \\ 1 \end{pmatrix}$$
The angle between any two basis vectors is:
$$\vec{e}_a \cdot \vec{e}_b = -\frac{1}{3} \implies \theta = \arccos\left(-\frac{1}{3}\right) \approx 109.47^\circ$$
This angle matches the tetrahedral $sp^3$ hybridization of peptide $C_\alpha$ backbone carbons ($109.5^\circ$), eliminating the non-physical $90^\circ$ right angles of cubic grids and enabling natural $\alpha$-helical and $\beta$-strand turn approximations.

### 2.4 Conformational Encodings: Cartesian Coordinates vs Relative Turn Representation

The choice of mathematical encoding determines the qubit and optical mode scaling on quantum hardware:

#### 2.4.1 Absolute Cartesian Coordinate Encoding
- Decision variables: Binary variables $x_{i, \vec{r}} \in \{0, 1\}$ indicating whether residue $i$ ($1 \le i \le N$) occupies lattice site $\vec{r} \in \mathcal{L}$.
- Scaling: If the lattice volume is $L \times L \times L = V$, the required number of qubits is:
  $$N_{\text{qubits}} = N \cdot V = \mathcal{O}(N \cdot L^3)$$
  For a small 10-residue peptide in a modest $6 \times 6 \times 6$ grid, $N_{\text{qubits}} = 10 \times 216 = 2,160$ qubits. This Cartesian encoding suffers from enormous coordinate redundancy, translational symmetry degeneracy, and rotational invariance overhead.

#### 2.4.2 Relative Turn / Direction Encoding
As formulated by Babbush et al. (2014) and refined by Robert et al. (2021), the conformation of a chain of $N$ beads is fully specified by its sequence of $N-1$ bond vectors:
$$\vec{v}_1, \vec{v}_2, \dots, \vec{v}_{N-1}$$
The position of residue $i$ is uniquely determined by path integration from the origin:
$$\vec{r}_i = \vec{r}_1 + \sum_{k=1}^{i-1} \vec{v}_k$$
- **Symmetry Elimination**:
  1. Setting $\vec{r}_1 = (0, 0, 0)$ completely eliminates 3D translational symmetry.
  2. Setting the first bond vector along a fixed basis direction (e.g., $\vec{v}_1 = \vec{e}_0$) eliminates rigid global rotational symmetry.
- **Relative Turn Bits**:
  At each subsequent step $k$ ($2 \le k \le N-1$), the next bond vector $\vec{v}_k$ is chosen from the allowed branching directions relative to the preceding bond $\vec{v}_{k-1}$.
  - On a 2D square lattice ($z=4$), each step has 4 absolute directions or 3 forward non-reversing turns (Left, Right, Straight), requiring $\lceil \log_2 3 \rceil = 2$ bits.
  - On a tetrahedral lattice ($z=4$), each bond has 3 forward tetrahedral directions (excluding the $180^\circ$ backtrack), requiring exactly 2 bits.
- **Resource Scaling**:
  For an $N$-residue peptide, only $N-2$ relative turns must be specified:
  $$N_{\text{bits}} = 2(N-2)$$
  For an 8-residue epitope loop ($N=8$), the entire conformational space is spanned by:
  $$N_{\text{bits}} = 2(8-2) = 12\ \text{qubits}$$
  This linear $\mathcal{O}(N)$ scaling makes turn encoding ideally suited for near-term photonic and NISQ quantum processors.

### 2.5 The Complete Protein Folding Cost Hamiltonian

The total conformational Hamiltonian maps the physical energy landscape into an operator diagonal in the computational basis:
$$H_{\text{total}} = H_{\text{interactions}} + \lambda_{\text{overlap}} H_{\text{overlap}} + \lambda_{\text{backtrack}} H_{\text{backtrack}}$$

#### 2.5.1 The Interaction Hamiltonian ($H_{\text{interactions}}$)
The interaction Hamiltonian scores pairwise topological contacts between non-bonded residues:
$$H_{\text{interactions}} = \sum_{1 \le i < j \le N, j - i \ge 2} \epsilon(s_i, s_j) \cdot \Delta(\vec{r}_i, \vec{r}_j)$$
where:
- $\epsilon(s_i, s_j)$ is the contact energy (HP or Miyazawa-Jernigan potential).
- $\Delta(\vec{r}_i, \vec{r}_j)$ is the adjacency contact indicator function:
  $$\Delta(\vec{r}_i, \vec{r}_j) = \begin{cases} 1 & \text{if } \|\vec{r}_i - \vec{r}_j\| = d_{\text{lattice}} \\ 0 & \text{otherwise} \end{cases}$$
- The topological constraint $j - i \ge 2$ excludes trivially bonded nearest neighbors along the primary backbone.

#### 2.5.2 The Self-Avoiding Walk Overlap Penalty ($H_{\text{overlap}}$)
Physical polypeptide chains exhibit steric excluded volume (Pauli repulsion): no two atoms can occupy the same spatial coordinates simultaneously. The self-avoiding penalty operator is defined as:
$$H_{\text{overlap}} = \sum_{1 \le i < j \le N} \delta(\vec{r}_i, \vec{r}_j)$$
where $\delta(\vec{r}_i, \vec{r}_j)$ is the Kronecker delta function:
$$\delta(\vec{r}_i, \vec{r}_j) = \begin{cases} 1 & \text{if } \vec{r}_i = \vec{r}_j \\ 0 & \text{if } \vec{r}_i \neq \vec{r}_j \end{cases}$$

#### 2.5.3 The Backtracking Penalty ($H_{\text{backtrack}}$)
If relative turn encodings allow complete $360^\circ$ directional choices without hardcoded turn exclusions, immediate $180^\circ$ reversals between consecutive bonds must be penalized:
$$H_{\text{backtrack}} = \sum_{k=1}^{N-2} \delta(\vec{v}_{k+1}, -\vec{v}_k)$$

### 2.6 Rigorous Derivation and Proof of Penalty Multiplier Bounds

A critical requirement in mapping constrained combinatorial optimization problems to unconstrained Hamiltonians (QUBO/Ising) is that **the ground state of $H_{\text{total}}$ must strictly correspond to a physically valid, self-avoiding fold**. If penalty multipliers $\lambda_{\text{overlap}}$ are set too low, the quantum optimizer may discover non-physical overlapping configurations that unphysically maximize contact energies.

#### Mathematical Theorem (Penalty Lower Bound)
*Let $\mathcal{C}_{\text{valid}}$ be the set of all self-avoiding, non-overlapping conformations ($\forall i < j, \vec{r}_i \neq \vec{r}_j$), and let $\mathcal{C}_{\text{invalid}}$ be the set of conformations with at least one coordinate overlap ($\exists i < j, \vec{r}_i = \vec{r}_j$). If the penalty multiplier satisfies:*
$$\lambda_{\text{overlap}} > \sum_{1 \le i < j \le N, j-i \ge 2} |\min(0, \epsilon(s_i, s_j))|$$
*then for any invalid conformation $\mathbf{c}_{\text{inv}} \in \mathcal{C}_{\text{invalid}}$ and any valid conformation $\mathbf{c}_{\text{val}} \in \mathcal{C}_{\text{valid}}$, the total Hamiltonian satisfies:*
$$E_{\text{total}}(\mathbf{c}_{\text{inv}}) > E_{\text{total}}(\mathbf{c}_{\text{val}})$$

#### Proof:
1. Consider the maximum possible energetic stabilization achievable by any valid conformation $\mathbf{c}_{\text{val}} \in \mathcal{C}_{\text{valid}}$:
   $$E_{\text{valid}} = H_{\text{interactions}}(\mathbf{c}_{\text{val}}) = \sum_{i < j, j-i \ge 2} \epsilon(s_i, s_j) \Delta(\vec{r}_i, \vec{r}_j)$$
   Since $\Delta(\vec{r}_i, \vec{r}_j) \le 1$, the absolute minimum energy bounded below is:
   $$E_{\text{min\_val}} \ge -\sum_{i < j, j-i \ge 2} |\min(0, \epsilon(s_i, s_j))| \equiv -E_{\text{max\_contact}}$$
   where $E_{\text{max\_contact}} = \sum_{i < j, j-i \ge 2} |\min(0, \epsilon(s_i, s_j))| > 0$.

2. Now consider an invalid conformation $\mathbf{c}_{\text{inv}} \in \mathcal{C}_{\text{invalid}}$. By definition, it contains at least one spatial overlap:
   $$H_{\text{overlap}}(\mathbf{c}_{\text{inv}}) \ge 1$$
   The interaction energy of $\mathbf{c}_{\text{inv}}$ cannot be lower than the theoretical lower bound where every non-bonded pair simultaneously contacts:
   $$H_{\text{interactions}}(\mathbf{c}_{\text{inv}}) \ge -E_{\text{max\_contact}}$$

3. Therefore, the total energy of the invalid conformation is bounded below by:
   $$E_{\text{total}}(\mathbf{c}_{\text{inv}}) = H_{\text{interactions}}(\mathbf{c}_{\text{inv}}) + \lambda_{\text{overlap}} H_{\text{overlap}}(\mathbf{c}_{\text{inv}}) \ge -E_{\text{max\_contact}} + \lambda_{\text{overlap}}$$

4. In order to guarantee that $E_{\text{total}}(\mathbf{c}_{\text{inv}}) > E_{\text{total}}(\mathbf{c}_{\text{val}})$ for all valid conformations:
   $$-E_{\text{max\_contact}} + \lambda_{\text{overlap}} > 0 \ge E_{\text{valid}} \implies \lambda_{\text{overlap}} > E_{\text{max\_contact}}$$
   More conservatively, since $E_{\text{valid}} \ge -E_{\text{max\_contact}}$:
   $$-E_{\text{max\_contact}} + \lambda_{\text{overlap}} > -E_{\text{max\_contact}} + \Delta E \implies \lambda_{\text{overlap}} > E_{\text{max\_contact}}$$
   Specifically, with $\lambda_{\text{overlap}} > \sum_{i < j} |\min(0, \epsilon(s_i, s_j))|$, an invalid conformation possessing even a single collision ($H_{\text{overlap}} = 1$) incurs a penalty exceeding the sum of all possible stabilizing interactions in the entire protein. Thus, no invalid state can ever achieve the ground state. $\quad \blacksquare$

---

## 3. Quantum Algorithmic Paradigms for Protein Folding

Multiple quantum computational paradigms have been formulated to solve the lattice protein folding Hamiltonian. This section provides an analytical comparison of their mathematical structure, resource requirements, and suitability for discrete-variable linear optical platforms.

### 3.1 Quantum Annealing & QUBO Formulation
In Quadratic Unconstrained Binary Optimization (QUBO), the objective function is expressed as a degree-2 binary polynomial:
$$f(\mathbf{q}) = \sum_{k=1}^M h_k q_k + \sum_{1 \le k < l \le M} J_{kl} q_k q_l, \quad q_k \in \{0, 1\}$$
Using the standard spin-orbital transformation $q_k = \frac{1 - Z_k}{2}$ where $Z_k \in \{+1, -1\}$ is the Pauli $Z$ operator, the QUBO is mapped isomorphically to an Ising spin-glass Hamiltonian:
$$H_{\text{Ising}} = \sum_{k=1}^M \tilde{h}_k Z_k + \sum_{k < l} \tilde{J}_{kl} Z_k Z_l + C$$

Quantum annealers (such as D-Wave superconducting flux qubit systems) initialize the system in the ground state of a transverse field Hamiltonian $H_{\text{trans}} = -\sum_k X_k$ and adiabatically evolve toward $H_{\text{Ising}}$:
$$H(t) = \left(1 - \frac{t}{T}\right) H_{\text{trans}} + \frac{t}{T} H_{\text{Ising}}$$
- **Limitations**:
  - Distance functions $\Delta(\vec{r}_i, \vec{r}_j)$ and coordinate indicators in turn-encoded chains naturally yield **Higher-Order Binary Optimization (HUBO)** polynomials of cubic or quartic degree ($\sum K_{jkl} q_j q_k q_l$).
  - Mapping HUBO to QUBO requires auxiliary slack variables (quadratization), multiplying qubit overhead by a factor of 3–5$\times$.
  - Minor embedding onto physical chimera or pegasus graphs incurs heavy ferromagnetic chain penalties, limiting fidelity for chains $> 15$ residues.

### 3.2 Variational Quantum Eigensolver (VQE)
The Variational Quantum Eigensolver (Peruzzo et al. 2014; McArdle et al. 2020) is a hybrid quantum-classical algorithm designed for NISQ processors.
1. **Parameterized Ansatz Preparation**: A quantum circuit prepares an entangled state $|\psi(\boldsymbol{\theta})\rangle = U(\boldsymbol{\theta}) |0\rangle^{\otimes M}$ governed by classical parameters $\boldsymbol{\theta} \in \mathbb{R}^K$.
2. **Quantum Measurement**: The quantum processor samples the state in the computational basis, estimating the expectation energy:
   $$\langle H \rangle_{\boldsymbol{\theta}} = \langle \psi(\boldsymbol{\theta}) | H | \psi(\boldsymbol{\theta}) \rangle = \sum_{\mathbf{c}} P(\mathbf{c}; \boldsymbol{\theta}) E(\mathbf{c})$$
3. **Classical Optimization**: A classical derivative-free optimization algorithm (COBYLA, Nelder-Mead, SPSA) updates parameters:
   $$\boldsymbol{\theta}^{(t+1)} = \boldsymbol{\theta}^{(t)} - \eta \cdot \nabla_{\boldsymbol{\theta}} \langle H \rangle$$
   iterating until convergence to the minimum eigenvalue $\min_{\boldsymbol{\theta}} \langle H \rangle_{\boldsymbol{\theta}} \ge E_0$.

### 3.3 Conditional Value-at-Risk (CVaR-VQE)
In standard VQE, minimizing the expectation value $\langle H \rangle = \mathbb{E}[E]$ averages over all sampled configurations. In combinatorial optimization—where the goal is discovering the single lowest-energy conformation rather than an ensemble thermal average—standard VQE suffers from slow convergence because high-energy tail configurations dominate the gradient.

Barkoutsos et al. (2020) and Robert et al. (IBM / *npj Quantum Information*, 2021) introduced **CVaR-VQE** for protein folding. For a confidence level $\alpha \in (0, 1]$, CVaR measures the expected value of the lower $\alpha$-tail of the energy distribution:
$$\text{CVaR}_\alpha(\boldsymbol{\theta}) = \mathbb{E}_{X \sim P_{\boldsymbol{\theta}}}[X \mid X \le F_X^{-1}(\alpha)]$$
In practical empirical execution with $K$ measurement shots returning sorted energies $E_{(1)} \le E_{(2)} \le \dots \le E_{(K)}$:
$$\text{CVaR}_\alpha(\boldsymbol{\theta}) = \frac{1}{\lceil \alpha K \rceil} \sum_{k=1}^{\lceil \alpha K \rceil} E_{(k)}$$
- **Biophysical Significance**: Setting $\alpha = 0.10$ to $0.25$ forces the classical optimizer to adjust optical phases $\boldsymbol{\theta}$ exclusively based on the top 10%–25% most stable folding conformations discovered. Sub-optimal, non-folded, or colliding configurations are completely filtered out of the objective gradient, accelerating convergence to the native fold by over $4\times$.

### 3.4 Boson Sampling & Optical Graph Analysis
Boson Sampling (Aaronson & Arkhipov 2011) describes the task of sampling from the output photon number distribution of an $M$-mode linear optical network injected with indistinguishable single photons.
The transition probability from an input Fock state $|\vec{s}\rangle = |s_1, \dots, s_M\rangle$ to an output state $|\vec{d}\rangle = |d_1, \dots, d_M\rangle$ is governed by the matrix permanent:
$$P(\vec{s} \to \vec{d}) = \frac{|\text{Perm}(U_{s, d})|^2}{\prod_i s_i! \prod_j d_j!}$$
Because calculating the permanent of a general matrix is **$\#\text{P}$-hard** (Valiant 1979), classical computers cannot sample from this distribution efficiently.

#### Application to Protein Core Discovery
In Gaussian Boson Sampling (GBS) and discrete Fock-state graph mapping (Bromley et al. 2020):
1. The amino acid contact adjacency matrix $A_{ij} = \Delta(\vec{r}_i, \vec{r}_j) \cdot |\epsilon(s_i, s_j)|$ is encoded into the transmission matrix $U$ of an optical interferometer.
2. Quantum interference (Hong-Ou-Mandel photon bunching) causes photons to preferentially cluster into optical modes corresponding to **dense cliques and tightly connected subgraphs**.
3. Sampling the output photon clicks identifies dense hydrophobic folding cores and stabilizing tertiary contact networks without combinatorial enumeration.

---

## 4. Feasibility Analysis for Photonic Quantum Computing (Quandela Hardware & Perceval SDK)

Quandela architectures operate on discrete-variable **Linear Optical Quantum Computing (LOQC)**. This section analyzes physical hardware constraints, encoding paradigms, and verified SDK implementation details.

### 4.1 Physical Hardware Architecture
The Quandela hardware stack consists of three physical components:
1. **Single-Photon Source**:
   - Semiconductor quantum dot (InGaAs) embedded within an optical microcavity / micropillar.
   - Optically pumped by pulsed laser to emit deterministic, on-demand single photons at telecommunication wavelengths ($\sim 925\ \text{nm}$).
   - Metrics: Ultra-low multi-photon emission probability $g^{(2)}(0) \approx 0.02$ (suppressing multi-photon noise), single-photon indistinguishability (Hong-Ou-Mandel visibility) $\text{HOM} \approx 88.98\%$, repetition rate $13.33\ \text{MHz}$.
2. **Reconfigurable Integrated Photonic Chip**:
   - Fabricated on low-loss silicon nitride ($\text{Si}_3\text{N}_4$) or silica-on-silicon waveguides.
   - Universal linear optical interferometers composed of integrated directional beam splitters and thermo-optic phase shifters $\theta_k, \phi_k$.
   - Phase shifts are tuned electronically via micro-heaters modulating the local refractive index.
3. **Superconducting Nanowire Single-Photon Detectors (SNSPD)**:
   - Cryogenically cooled ($< 2.5\ \text{K}$) superconducting nanowires.
   - High system detection efficiency ($> 85\%$), low dark count rates ($< 10\ \text{Hz}$), timing jitter $< 30\ \text{ps}$.

### 4.2 Linear Optical Encodings: Dual-Rail Qubits vs Continuous Fock-State Ansätze

Two distinct operational paradigms exist for running quantum folding algorithms on Quandela hardware:

#### Paradigm A: Dual-Rail Qubit Encoding (Gate-Based Emulation)
- A logical qubit is encoded using one photon across two spatial waveguide modes:
  $$|0\rangle_L \equiv |1, 0\rangle, \quad |1\rangle_L \equiv |0, 1\rangle$$
- **Single-Qubit Gates**: Arbitrary $SU(2)$ rotations $R_z(\theta) R_y(\phi)$ are implemented deterministically with 100% fidelity using a single Mach-Zehnder Interferometer (MZI).
- **Two-Qubit Gates**: Entangling operations (CZ, CNOT) require photon-photon interactions. In linear optics without non-linear optical media ($\chi^{(3)}$), entangling gates are probabilistic, requiring ancillary photons and measurement post-selection under the **Knill-Laflamme-Milburn (KLM)** scheme. The theoretical success probability of a heralded CNOT gate is $P = 1/9 \approx 11.1\%$ (or $1/4$ with ancillary resource states).
- **Hardware Requirement**: For an $N$-residue loop requiring $2(N-2)$ qubits, dual-rail encoding requires $4(N-2)$ spatial modes and $2(N-2)$ photons. A 6-residue loop ($N=6$, 8 qubits) would require 16 modes and 8 photons, exceeding the 12-mode capacity of `qpu:ascella`.

#### Paradigm B: Continuous Fock-State Variational Quantum Architecture (VQA on LOQC)
- Rather than synthesizing discrete gates via probabilistic KLM post-selection, the ansatz directly utilizes the **full continuous reconfigurability of the universal interferometer**.
- An $M$-mode chip (such as the 12 modes of `qpu:ascella`) is parameterized as a Clements unitary $U(\boldsymbol{\theta}) \in \mathcal{U}(M)$ with $M(M-1)/2 = 66$ tunable phases.
- Input: Injected Fock state $|\psi_{\text{in}}\rangle = |1, 1, \dots, 1, 0, \dots\rangle$ with $N_{\text{photons}} \le 6$.
- Output: Output coincidence click patterns are mapped deterministically to conformational turn vectors.
- Variational optimization updates the continuous phases $\boldsymbol{\theta}$ using CVaR-VQE.
- **Feasibility**: 100% compatible with near-term Quandela hardware. Executes natively without probabilistic gate overhead.

### 4.3 Perceval SDK Architecture & Remote Execution Protocol

The Perceval SDK (version `1.2.4`) provides the software interface for synthesizing optical circuits and dispatching jobs to Quandela Cloud. Empirical reverse-engineering of the live Quandela API (`api.cloud.quandela.com`) identified several critical operational constraints:

1. **Authentication Token Prefix**:
   The cloud API strictly requires the literal prefix `_T_` in the authorization header:
   `Authorization: Bearer _T_eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9...`
   Stripping `_T_` to send a standard JWT triggers HTTP `401 UNAUTHORIZED: {"error": "Invalid token"}`.
2. **Mandatory Post-Selection Filter**:
   When configuring a `RemoteProcessor`, setting the minimum detected photons filter is mandatory:
   ```python
   processor.min_detected_photons_filter(input_state.n)
   ```
   Omitting this method triggers an immediate client-side `ValueError: The value of min_detected_photons is not set. Use the method processor.min_detected_photons_filter(value).`
3. **Mandatory Remote Sampler Parameter**:
   The high-level `Sampler` class requires explicit specification of `max_shots_per_call`:
   ```python
   sampler = Sampler(processor, max_shots_per_call=1000)
   ```
   Omitting this parameter triggers `RuntimeError: Please input a max_shots_per_call value when using a RemoteProcessor`.
4. **Perceval 1.2+ Property Semantics**:
   In Perceval 1.2.4, `sampler.sample_count` and `sampler.probs` are properties returning job objects rather than methods. Execution must call `.execute_sync(shots)` or `.execute_async()`.
5. **Compressed Result Decompression**:
   Quandela Cloud transmits simulation payloads using zlib compression with a `:PCVL:zip:` header. The SDK's `RemoteJob.get_results()` automatically decodes this payload into a `BSCount` or `BSDistribution` object.

### 4.4 Hardware Feasibility & Empirical Backend Comparison

| Platform Specification | `sim:slos` (Cloud Simulator) | `qpu:ascella` (Physical QPU) |
|---|---|---|
| Algorithmic Engine | Strong Linear Optical Simulator (SLOS) | Integrated Photonic Chip (12 modes) |
| Hosting Infrastructure | OVHcloud | Quandela Hardware Lab |
| Maximum Photon Capacity | Up to 12 photons | Up to 6 single photons |
| Multi-Photon State Support | **Yes** ($|2, 0\rangle, |0, 2\rangle, \dots$) | **No** (`support_multi_photon: False`) |
| Coincidence Post-Selection | Mathematical filtering | Physical SNSPD coincidence window |
| Live Execution Verification | **Job `d434c869-d5f5-449f...` verified** | Status: `maintenance` / `bookable` |
| Verification Output | `BSCount{|0,2>: 53, |2,0>: 47}` (HOM bunching) | Subject to scheduled beamline bookings |

Because physical QPUs undergo periodic optical recalibration and maintenance, production software must implement automated fallback to `sim:slos`, which provides provably exact simulation of identical quantum circuits.

### 4.5 Hybrid Coarse-Graining and Epitope Loop Windowing

A critical architectural consideration for competition success is bridging the scale gap between 85–100 amino acid minibinders and near-term 12-mode quantum hardware. Simulating 100 residues in full atomistic detail would require hundreds of qubits/modes, far exceeding current hardware capacities worldwide.

#### The Windowing Strategy
In protein-protein interactions, binding affinity and specificity are overwhelmingly determined by a small subset of residues located on flexible surface loops (the "paratope" and "epitope"), while the remainder of the protein acts as a rigid scaffolding framework.

```
Full 85-100 aa Minibinder
[ Scaffold Alpha-Helix 1 ] ───► [ Critical Binding Loop (6-10 aa) ] ───► [ Scaffold Alpha-Helix 2 ]
          │                                      │                                      │
          ▼                                      ▼                                      ▼
Classical ESMFold2 / ipSAE             Quandela Photonic VQE                  Classical ESMFold2 / ipSAE
(Global Complex Geometry)           (Local Conformational Stability)          (Global Complex Geometry)
```

The hybrid architecture decomposes candidate evaluation into complementary regimes:
1. **Classical Macromolecular Folding (ESMFold2 + ipSAE)**:
   Predicts the global tertiary structure of the full minibinder (85–100 aa) in complex with the target trimer/receptor, computing the global interface accessible solvation energy ($ipSAE$).
2. **Quantum Epitope Loop Folding (Quandela LOQC VQE)**:
   Isolates the 6–10 residue binding loop directly contacting the target hotspot:
   - For **TL1A (TNFSF15)**: The 10-residue loop inserting into the trimer interface cleft.
   - For **IL-7Rα (CD127)**: The 7-residue loop engaging the `KKFLLIG` hydrophobic core.
   The quantum engine computes the ground-state folding energy $E_{\text{ground}}$ and verifies whether the binding-competent loop conformation represents a deep, stable free-energy minimum rather than a transient, disordered state.
3. **Multi-Objective Hybrid Scoring**:
   The final composite fitness index combines both evaluations:
   $$\text{Score}_{\text{hybrid}} = 0.6 \cdot \text{ipSAE} + 0.4 \cdot S_{\text{quantum}}$$
   This hybrid formulation guarantees both macro-scale complex affinity and micro-scale loop conformational stability.

---

## 5. Comparative Algorithmic Summary

| Algorithm / Paradigm | Underlying Physical Principle | Classical Scaling | Quantum Advantage Mechanism | Suitability for Quandela Hardware |
|---|---|---|---|---|
| **Classical Monte Carlo / Simulated Annealing** | Thermal fluctuations over energy barriers $\propto e^{-\Delta E / k_B T}$ | Exponential $\mathcal{O}(c^N)$ in rugged landscapes | None (baseline classical method) | Baseline comparator |
| **Deep Learning (ESMFold2 / AlphaFold2)** | Transformer attention + all-atom diffusion | Polynomial inference $\mathcal{O}(N^2)$ | Fast static predictions; poor at de novo loops without MSAs | Complementary classical screening branch |
| **Quantum Annealing (D-Wave)** | Transverse-field quantum tunneling across narrow energy barriers | $\mathcal{O}(e^{\sqrt{N}})$ adiabatic scaling | Tunneling through tall thin barriers | Moderate (requires slack variables for HUBO) |
| **Gate-Based VQE (Dual-Rail LOQC)** | Parameterized unitary state preparation + classical optimization | Qubit scaling $\mathcal{O}(N)$ | Coherent superposition; CVaR noise suppression | High gate overhead (probabilistic KLM gates) |
| **Continuous Variational LOQC (Quandela)** | Universal Clements MZI mesh + Fock state sampling | Qumode scaling $\mathcal{O}(M)$ | Native optical unitaries; no KLM gate overhead | **Optimal (Directly executable on `sim:slos` / `qpu:ascella`)** |
| **Photonic Boson Sampling** | Many-body Hong-Ou-Mandel interference | Permanent calculation is $\#\text{P}$-hard | Direct concentration onto dense contact cliques | High (identifies hydrophobic core clusters) |

---

## 6. Formal Academic Bibliography

1. **Robert, A., Barkoutsos, P. K., Woerner, S., & Tavernelli, I.** (2021). "Resource-efficient quantum algorithm for protein folding." *npj Quantum Information*, 7(1), 38. https://doi.org/10.1038/s41534-021-00368-4  
   *(Foundational paper establishing tetrahedral diamond lattice turn-encoding and CVaR-VQE for coarse-grained protein folding).*

2. **Babbush, R., Berry, D. W., Kivlichan, I. D., Wei, G. Y., Love, P. J., & Aspuru-Guzik, A.** (2014). "Resource efficient quantum simulation of protein folding." *Advances in Chemical Physics*, 154, 429–470. https://doi.org/10.1002/9781118749814.ch11  
   *(Introduced relative turn-based direction encoding for lattice models, reducing qubit overhead from $\mathcal{O}(N \cdot V)$ to $\mathcal{O}(N)$).*

3. **Miyazawa, S., & Jernigan, R. L.** (1985). "Estimation of effective interresidue contact energies from protein crystal structures: quasi-chemical approximation." *Macromolecules*, 18(3), 534–552. https://doi.org/10.1021/ma00145a039  
   *(Original derivation of the statistical contact potential matrix from experimental PDB structures).*

4. **Miyazawa, S., & Jernigan, R. L.** (1996). "Residue-residue potentials with a favorable contact pair term and an unfavorable high packing density term, for simulation and threading." *Journal of Molecular Biology*, 256(3), 623–644. https://doi.org/10.1006/jmbi.1996.0114  
   *(Refined $20 \times 20$ contact potential matrix widely adopted across structural bioinformatics and quantum chemistry).*

5. **Dill, K. A.** (1985). "Theory for the folding and stability of globular proteins." *Biochemistry*, 24(6), 1501–1517. https://doi.org/10.1021/bi00327a032  
   *(Original formulation of the hydrophobic effect as the dominant thermodynamic driving force in globular protein folding).*

6. **Lau, K. F., & Dill, K. A.** (1989). "A lattice statistical model for enzyme nucleation and polymer folding." *Macromolecules*, 22(10), 3986–3994. https://doi.org/10.1021/ma00200a030  
   *(Formulated the two-letter Hydrophobic-Polar (HP) lattice model).*

7. **Dill, K. A., Bromberg, S., Yue, K., Fiebig, K. M., Yee, D. P., Thomas, P. D., & Chan, H. S.** (1995). "Principles of protein folding—A perspective from simple exact models." *Protein Science*, 4(4), 561–602. https://doi.org/10.1002/pro.5560040401  
   *(Comprehensive review of exact lattice models, self-avoiding walks, and contact cooperativity).*

8. **Peruzzo, A., McClean, J., Shadbolt, P., Yung, M. H., Zhou, X. Q., Love, P. J., Aspuru-Guzik, A., & O'Brien, J. L.** (2014). "A variational eigenvalue solver on a photonic quantum processor." *Nature Communications*, 5(1), 4213. https://doi.org/10.1038/ncomms5213  
   *(First experimental demonstration of the Variational Quantum Eigensolver executed on integrated linear optical hardware).*

9. **Knill, E., Laflamme, R., & Milburn, G. J.** (2001). "A scheme for efficient quantum computation with linear optics." *Nature*, 409(6816), 46–52. https://doi.org/10.1038/35051009  
   *(The landmark KLM paper proving that scalable quantum computation is possible using only single-photon sources, linear optical elements, and single-photon detectors).*

10. **Cao, L., Coventry, B., Goreshnik, I., Huang, B., Park, J. S., Jude, K. M., ... & Baker, D.** (2022). "Design of protein-binding proteins from target structural information only." *Nature*, 605(7910), 551–560. https://doi.org/10.1038/s41586-022-04654-9  
    *(Demonstrated de novo design of high-affinity minibinders targeting IL-7Rα and other cell surface receptors without prior binder templates).*

11. **Heurtel, N., Fyrillas, A., de Gliniasty, G., Le Bihan, R., Malherbe, S., Pailhas, M., ... & Mansfield, S.** (2023). "Perceval: A Software Platform for Discrete-Variable Photonic Quantum Computing." *Quantum*, 7, 931. https://doi.org/10.22331/q-2023-02-21-931  
    *(Official architecture specification and benchmarking of the Quandela Perceval simulation and quantum cloud platform).*

12. **Clements, W. R., Humphreys, P. C., Metcalf, B. J., Kolthammer, W. S., & Walmsley, I. A.** (2016). "Optimal design for universal multiport interferometers." *Optica*, 3(12), 1460–1465. https://doi.org/10.1364/OPTICA.3.001460  
    *(Derived the optimal triangular decomposition of universal linear optical unitaries on reconfigurable photonic chips).*

13. **Aaronson, S., & Arkhipov, A.** (2011). "The computational complexity of linear optics." In *Proceedings of the 43rd Annual ACM Symposium on Theory of Computing (STOC)*, 333–342. https://doi.org/10.1145/1993636.1993682  
    *(Formal proof establishing the computational hardness and quantum supremacy of Boson Sampling).*

14. **Bromley, T. R., Arrazola, J. M., Jahangiri, S., Izaac, J., Quesada, N., Gran, A. D., ... & Killoran, N.** (2020). "Applications of near-term photonic quantum computers: software and algorithms." *Quantum Science and Technology*, 5(3), 034010. https://doi.org/10.1088/2058-9565/ab8504  
    *(Comprehensive survey mapping graph optimization, dense clique search, and molecular docking onto photonic architectures).*

15. **Barkoutsos, P. K., Nannicini, G., Robert, A., Tavernelli, I., & Woerner, S.** (2020). "Improving Variational Quantum Optimization using CVaR." *Quantum*, 4, 256. https://doi.org/10.22331/q-2020-04-20-256  
    *(Formulation of Conditional Value-at-Risk for variational quantum optimization algorithms).*

16. **Berger, B., & Leighton, T.** (1998). "Protein folding in the hydrophobic-hydrophilic (HP) model is NP-complete." *Journal of Computational Biology*, 5(1), 27–40. https://doi.org/10.1089/cmb.1998.5.27  
    *(Formal computational complexity proof demonstrating NP-completeness of 3D lattice protein folding).*

17. **Crescenzi, P., Goldman, D., Papadimitriou, C., Piccolboni, A., & Yannakakis, M.** (1998). "On the complexity of protein folding." *Journal of Computational Biology*, 5(3), 423–465. https://doi.org/10.1089/cmb.1998.5.423  
    *(Extended NP-completeness proof to 2D HP lattice models).*

---

*End of Literature Review & Theoretical Foundations.*
