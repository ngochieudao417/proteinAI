# Architecture Specification: Quantum Protein Folding Pipeline & Competition Submission Platform

**Project:** proteinAI  
**Milestone:** 1 (M1) — Architectural Blueprints & Theoretical Foundations  
**Author:** QAAA Lab Core Team  
**Date:** November 2026  
**Document Status:** Production-Grade Specification (Approved)  
**Target Platform:** Quandela Cloud (`sim:slos`, `qpu:ascella`) & SJTU AI Protein Platform ("Changjiang Cup" 2026)  

---

## 1. Executive Summary & System Context

The **proteinAI** platform is an artificial intelligence-driven molecular engineering system engineered for de novo design, in-silico screening, and biophysical evaluation of therapeutic protein candidates. For the **"Changjiang Cup" 2026 Global AI Protein Design Challenge** (2026 “长江杯”全球AI蛋白质设计挑战赛) hosted by Shanghai Jiao Tong University (SJTU), the system targets two high-value disease targets:
1. **Human Tumor Necrosis Factor Superfamily Member 15 (TL1A / TNFSF15)**: Target ID `ea930a01-20a9-4bed-b70f-e658bd9ba3ab`, a homotrimeric cytokine implicated in inflammatory bowel disease (Crohn's disease and ulcerative colitis).
2. **Human Interleukin-7 Receptor Subunit Alpha (IL-7Rα / CD127)**: Target ID `1347d0e0-7d3d-4aa8-bcf1-c2dab9d8f1df`, a cell-surface receptor controlling lymphocyte survival and proliferation, central to severe autoimmunity and acute lymphoblastic leukemia.

### 1.1 Architectural Evolution

The existing proteinAI architecture operates an in-silico screening pipeline coupling deep learning structural prediction models with geometric interface scoring:
- **Classical Generative Pipeline:** De novo minibinder generation via **BindCraft**, **RFdiffusion**, and **ProteinMPNN**, constrained strictly to **85–100 amino acids**. This length window provides a mathematical immunity against competition novelty violations ($< 80\%$ alignment coverage against all 122,273 known de novo binders from Cao et al. 2022).
- **Classical In-Silico Screening Engine:** **ESMFold2** (Biohub 6.58B parameter all-atom diffusion model) running with `fp8` ESMC representations on NVIDIA L4 GPUs (Google Colab Pro).
- **Interface Evaluation Metric:** **ipSAE** (interface predicted Solvated Accessible Energy, Dunbrack Lab), which provides high discriminatory separation between functional binders ($ipSAE > 0.75$) and non-binding controls ($ipSAE \approx 0.00$).

```
[Candidate Minibinder Pool (85-100 aa)]
                  │
  ┌───────────────┴───────────────┐
  ▼                               ▼
[Classical Branch]       [Quantum Branch (NEW)]
ESMFold2 + ipSAE         Quandela LOQC / VQE
(Complex Geometry)       (Epitope Folding Stability)
  │                               │
  └───────────────┬───────────────┘
                  ▼
   [Hybrid Multi-Objective Ranker]
    Score = w1*ipSAE + w2*S_quantum
                  │
                  ▼
   [Automated Platform Submitter]
     SJTU REST API Client
```

### 1.2 The Quantum Upgrade Mandate

Classical machine learning models like ESMFold2 predict tertiary backbones based on statistical representations learned from the Protein Data Bank (PDB). However, deep learning models often struggle to resolve **non-local topological contact cooperativity**, subtle sidechain packing transitions, and conformational entropy within flexible binding epitope loops.

The **Quantum Folding Pipeline** integrates discrete-variable photonic quantum computing via **Perceval** (Quandela SDK) and **Quandela Cloud RemoteProcessors**. Rather than attempting full all-atom simulation of 100-residue proteins on near-term hardware, the system executes **hybrid coarse-grained windowing**:
1. Critical interaction loops (e.g., the 7-residue hotspot loop `KKFLLIG` in IL-7Rα or the 11-residue receptor-binding groove loop in TL1A) are extracted.
2. The conformational space is mapped onto coarse-grained discrete lattices (Tetrahedral $sp^3$ and 3D cubic).
3. Contact energy Hamiltonians parameterized by the empirical **Miyazawa-Jernigan (MJ)** $20 \times 20$ statistical potential matrix are constructed.
4. Conformations are optimized using a **Variational Photonic Eigensolver (VQE)** executing on Quandela Cloud (`sim:slos` simulator and `qpu:ascella` photonic processor).
5. The resulting ground-state energy provides a normalized quantum stability score $S_{\text{quantum}}$, which fuses into the final candidate selection ranking.
6. The top 20 candidate sequences per track are packaged, verified against platform validation rules, and submitted via the automated SJTU platform client.

---

## 2. Component Architecture & End-to-End Data Flow

The upgraded proteinAI architecture is structured into five cohesive tiers:
- **Generative & Screening Tier**: Generates candidate backbones and computes classical structural metrics.
- **Novelty & Safety Filter Tier**: Validates sequence length, canonical residues, motif absence, and Smith-Waterman novelty.
- **Quantum Folding Engine (`quantum_folding/`)**: Coarse-grains peptide loops, constructs turn-based Hamiltonians, synthesizes linear optical circuits, and executes VQE on Quandela Cloud.
- **Hybrid Multi-Objective Ranking Tier**: Fuses classical $ipSAE$ and quantum stability score $S_{\text{quantum}}$ into a unified ranking index.
- **Platform Submission & Verification Tier (`tools/submit_sjtu.py`)**: Authenticates with SJTU API, discovers active target metadata, pre-validates CSV batches, executes submission uploads, and monitors evaluation status.

### 2.1 Complete System Architecture Diagram (ASCII)

```
===================================================================================================================
                                      PROTEINAI SYSTEM ARCHITECTURE
===================================================================================================================

 [Candidate Generation]
   BindCraft / RFdiffusion / ProteinMPNN
   Produces 85-100 aa minibinders (e.g., 50 TL1A + 50 IL-7Ra candidates in candidate_pool_for_screening.csv)
             │
             ▼
 ┌────────────────────────────────────────────────────────────────────────────────────────────────────────┐
 │ NOVELTY & FILTER TIER (tools/novelty_check.py, tools/candidate_generator.py)                          │
 │ - Length verification: 85 <= L <= 100 aa (Guarantees < 80% coverage vs 122k known binders)            │
 │ - Canonical 20 amino acids check: ACDEFGHIKLMNPQRSTVWY                                                │
 │ - Motif filters: Zero unpaired cysteines; zero deamidation motifs (NG, NS, DG, DP)                    │
 │ - Smith-Waterman local alignment against data/known_binders/ (>122k records)                           │
 └───────────────────────────────────┬────────────────────────────────────────────────────────────────────┘
                                     │ Filtered Valid Candidates
        ┌────────────────────────────┴─────────────────────────────┐
        ▼                                                          ▼
 ┌─────────────────────────────────────────┐    ┌─────────────────────────────────────────────────────────┐
 │ CLASSICAL SCREENING TIER                │    │ QUANTUM FOLDING ENGINE TIER (quantum_folding/)          │
 │ (notebooks/03_esmfold2_screening.ipynb) │    │                                                         │
 │                                         │    │ 1. Epitope Loop Windowing (quantum_folding/pipeline.py) │
 │ - Model: ESMFold2 (6.58B params, fp8)   │    │    Extracts 6-10 aa binding loop sequences              │
 │ - Hardware: NVIDIA L4 GPU (Colab Pro)   │    │                                                         │
 │ - Input: Target MSA + Candidate Seq     │    │ 2. Lattice Encoding (quantum_folding/lattice.py)        │
 │ - Interface Metric: ipSAE               │    │    Tetrahedral / 3D cubic relative turns (2(N-2) bits)  │
 │   ipSAE = max_j [1 / (1 + (PAE/d0)^2)]  │    │                                                         │
 │ - Controls:                             │    │ 3. Cost Hamiltonian (quantum_folding/hamiltonian.py)    │
 │   * 7S5B Binder: ipSAE = 0.7728 (PASS)  │    │    H_tot = H_MJ + lambda_overlap*H_ov + lambda_bt*H_bt  │
 │   * Scrambled:   ipSAE = 0.0000 (FAIL)  │    │    Miyazawa-Jernigan 20x20 contact matrix               │
 │                                         │    │                                                         │
 │ Output: ipSAE Score in [0.0, 1.0]       │    │ 4. Photonic Circuit (quantum_folding/circuit.py)        │
 │                                         │    │    Perceval MZI network, Clements decomposition         │
 └────────────────────┬────────────────────┘    │                                                         │
                      │                         │ 5. Cloud RemoteProcessor (quandela_backend.py)          │
                      │                         │    Quandela Cloud (sim:slos / qpu:ascella)              │
                      │                         │    Token: _T_eyJ... (Bearer), max_shots=1000            │
                      │                         │                                                         │
                      │                         │ 6. Variational Solver (quantum_folding/vqe.py)          │
                      │                         │    CVaR-VQE optimization loop -> Ground state energy    │
                      │                         │                                                         │
                      │                         │ Output: Quantum Stability Score S_quantum in [0.0, 1.0] │
                      │                         └────────────────────────────┬────────────────────────────┘
                      │                                                      │
                      └──────────────────────┬───────────────────────────────┘
                                             ▼
 ┌────────────────────────────────────────────────────────────────────────────────────────────────────────┐
 │ HYBRID MULTI-OBJECTIVE RANKING TIER (quantum_folding/pipeline.py)                                      │
 │                                                                                                        │
 │ Composite Fitness Function:                                                                            │
 │   Score_hybrid = w1 * ipSAE + w2 * S_quantum  (w1 = 0.6, w2 = 0.4)                                    │
 │                                                                                                        │
 │ Select Top 20 Candidates per Target Track:                                                             │
 │   - Track 5 (TL1A): 20 sequences (40-100 aa)                                                           │
 │   - Track 7 (IL-7Ra): 20 sequences (40-100 aa)                                                         │
 │                                                                                                        │
 │ Package into data/submission_ready/{target}_submission.csv                                             │
 └───────────────────────────────────┬────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
 ┌────────────────────────────────────────────────────────────────────────────────────────────────────────┐
 │ SJTU PLATFORM SUBMISSION CLIENT (tools/submit_sjtu.py)                                                 │
 │                                                                                                        │
 │ 1. Auth: POST /api/auth/login (email, password, portal="contestant") -> JWT Access & Refresh Tokens   │
 │ 2. Target Query: GET /api/submissions/targets -> Verify UUIDs, deadlines, length bounds [40, 100]       │
 │ 3. Client Validation: Format check, column headers (name,sequence), duplicate check, canonical AAs     │
 │ 4. Submission Upload: POST /api/submissions/csv                                                        │
 │    Payload: {target_id, action="save"|"confirm", agreed=True, csv, abstract, note}                     │
 │ 5. Status Tracking: GET /api/submissions/latest, GET /api/submissions?grouped=1                       │
 └───────────────────────────────────┬────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
                      [SJTU AI Protein Platform]
                 https://aiprotein.sjtu.edu.cn/app
            (Leaderboard & Wet-Lab BLI/SPR Kd Measurement)
===================================================================================================================
```

### 2.2 Detailed Mermaid Data Flow Diagram

```mermaid
flowchart TD
    subgraph GenerativeTier [Generative & Novelty Screening Tier]
        Cands[Candidate Pool: 100 Minibinders<br/>85-100 aa in candidate_pool_for_screening.csv]
        NovCheck{Novelty & Safety Filter<br/>tools/novelty_check.py}
        Cands --> NovCheck
        NovCheck -- "Coverage >= 80% or Unpaired Cys" --> Reject[Rejected Candidate]
        NovCheck -- "PASS: Coverage < 80% & Mature Sequence" --> FilteredCands[Validated Candidates]
    end

    subgraph ClassicalScreening [Classical In-Silico Screening]
        ESM[ESMFold2 All-Atom Diffusion<br/>Colab L4 GPU, fp8 ESMC]
        ipSAE[ipSAE Interface Energy Calculation<br/>Dunbrack Lab Formulation]
        FilteredCands --> ESM
        ESM --> ipSAE
        ipSAE --> ScoreClassic[ipSAE Score: 0.0 to 1.0]
    end

    subgraph QuantumFoldingEngine [Quantum Folding Engine: quantum_folding/]
        LoopWin[Epitope Loop Windowing<br/>Extract 6-10 aa Critical Interaction Turn]
        Lattice[Lattice Geometry Mapping<br/>Tetrahedral/Diamond sp3 & 3D Cubic]
        Ham[Hamiltonian Synthesis<br/>Miyazawa-Jernigan 20x20 Matrix + Overlap Penalties]
        Circuit[Perceval Linear Optical Circuit<br/>Parameterized MZIs / Clements Decomposition]
        Quandela[Quandela Cloud RemoteProcessor<br/>Token: _T_eyJ... | Backend: sim:slos / qpu:ascella]
        VQE[CVaR-VQE Variational Optimizer<br/>COBYLA / Nelder-Mead Optimization]
        ScoreQ[Quantum Stability Score: S_quantum]

        FilteredCands --> LoopWin
        LoopWin --> Lattice
        Lattice --> Ham
        Ham --> Circuit
        Circuit --> Quandela
        Quandela --> VQE
        VQE --> ScoreQ
    end

    subgraph HybridRanking [Hybrid Multi-Objective Ranking]
        Ranker[Multi-Objective Fitness Evaluator<br/>Score = 0.6*ipSAE + 0.4*S_quantum]
        ScoreClassic --> Ranker
        ScoreQ --> Ranker
        Ranker --> Top20[Top 20 Sequences per Track<br/>TL1A and IL-7Ra]
        Top20 --> CSVOut[data/submission_ready/*.csv]
    end

    subgraph CompetitionSubmission [SJTU Competition Submission: tools/submit_sjtu.py]
        Auth[POST /api/auth/login<br/>portal: 'contestant']
        TargetQuery[GET /api/submissions/targets<br/>TL1A: ea930a01... | IL-7Ra: 1347d0e0...]
        ValidateCSV[Strict Client Pre-Validation<br/>Header, Unique Names, Canonical AAs, 40-100 aa]
        SubmitAPI[POST /api/submissions/csv<br/>action: 'confirm', agreed: true, csv text]
        StatusTrack[GET /api/submissions/latest<br/>Verify Upload Batch & Shortlisting]

        CSVOut --> ValidateCSV
        Auth --> TargetQuery
        TargetQuery --> ValidateCSV
        ValidateCSV --> SubmitAPI
        SubmitAPI --> StatusTrack
    end
```

---

## 3. New Python Package Specification: `quantum_folding`

The quantum folding engine is encapsulated within a modular Python package located at `/home/khang/qaaa-lab/proteinAI/quantum_folding/`. This section defines the interfaces, mathematical contracts, and implementation blueprints for each module.

```
quantum_folding/
├── __init__.py           # Package namespace, versioning, and top-level exports
├── lattice.py            # Discrete spatial lattices, coordinate math, and turn encodings
├── hamiltonian.py        # Miyazawa-Jernigan contact potentials and cost Hamiltonian synthesis
├── circuit.py            # Perceval parameterized linear optical circuit synthesis
├── quandela_backend.py   # Quandela Cloud RemoteProcessor client with token and quota handling
├── vqe.py                # Variational Photonic Eigensolver (VQE & CVaR-VQE) implementation
└── pipeline.py           # Batch candidate screening and hybrid scoring pipeline
```

### 3.1 `quantum_folding/__init__.py`
- **Purpose**: Defines package metadata, version string (`__version__ = "0.1.0"`), and exports primary public classes and functions.
- **Exports**:
  - `LatticeProtein`, `LatticeGeometry` from `.lattice`
  - `FoldingHamiltonian`, `MiyazawaJerniganPotential` from `.hamiltonian`
  - `PercevalCircuitBuilder`, `DualRailEncoding` from `.circuit`
  - `QuandelaCloudClient`, `RemoteExecutionConfig` from `.quandela_backend`
  - `VQESolver`, `CVaREstimator`, `VQEResult` from `.vqe`
  - `QuantumFoldingPipeline`, `HybridScorer` from `.pipeline`

### 3.2 `quantum_folding/lattice.py`
- **Purpose**: Translates primary amino acid sequences into discrete coordinate models, managing coordinate vectors, relative turn encodings, and spatial collision checks.
- **Data Models**:
  - `ResidueType`: Enum representing canonical amino acids and coarse-grained HP classifications (`H` for Hydrophobic: `A, V, I, L, M, F, W, P, C`; `P` for Polar: `G, S, T, Y, N, Q, D, E, K, R, H`).
  - `LatticeGeometry`: Enum (`SQUARE_2D`, `CUBIC_3D`, `TETRAHEDRAL_3D`).
- **Mathematical Formulations**:
  - **Tetrahedral (Diamond) Lattice**: Coordination number $z=4$. Basis vectors:
    $$\vec{v}_0 = \frac{1}{\sqrt{3}}(1, 1, 1), \quad \vec{v}_1 = \frac{1}{\sqrt{3}}(1, -1, -1), \quad \vec{v}_2 = \frac{1}{\sqrt{3}}(-1, 1, -1), \quad \vec{v}_3 = \frac{1}{\sqrt{3}}(-1, -1, 1)$$
    The angle between any two vectors satisfies $\cos \theta = -1/3 \implies \theta \approx 109.47^\circ$, exactly reproducing the $sp^3$ tetrahedral bonding geometry of polypeptide $C_\alpha - C - N$ backbones.
  - **Turn / Relative Direction Encoding**:
    For a chain of length $N$, the conformation is described by $N-1$ bond steps. Fixing $\vec{r}_1 = (0,0,0)$ and $\vec{v}_1 = \vec{v}_0$ eliminates global translation and rigid rotation. The remaining $N-2$ steps are encoded using relative turns.
    On a tetrahedral lattice, each turn from bond $\vec{v}_k$ can choose among 3 forward non-backtracking directions, encoded with 2 bits:
    $$\text{Total conformation bits} = 2(N-2)$$
    This achieves exponential space reduction compared to Cartesian coordinate grids $\mathcal{O}(N \cdot |\mathcal{L}|)$.
- **Class `LatticeProtein`**:
  - `__init__(sequence: str, geometry: LatticeGeometry = LatticeGeometry.TETRAHEDRAL_3D)`
  - `compute_coordinates(turns: Sequence[int]) -> np.ndarray`: Maps a turn bitstring/integer sequence into 3D Cartesian coordinates $(x_i, y_i, z_i)_{i=1}^N$.
  - `detect_contacts(coordinates: np.ndarray) -> List[Tuple[int, int]]`: Returns non-bonded pairs $(i, j)$ with $|i - j| \ge 2$ that occupy adjacent lattice sites ($\|\vec{r}_i - \vec{r}_j\| = 1$).
  - `has_self_overlap(coordinates: np.ndarray) -> bool`: Returns `True` if any two residues share identical lattice coordinates ($\vec{r}_i = \vec{r}_j$).

### 3.3 `quantum_folding/hamiltonian.py`
- **Purpose**: Constructs the objective cost Hamiltonian for arbitrary amino acid sequences, encapsulating non-bonded pairwise contact energies, volume exclusion constraints, and chain connectivity.
- **Mathematical Formulation**:
  $$H_{\text{total}} = H_{\text{interactions}} + \lambda_{\text{overlap}} H_{\text{overlap}} + \lambda_{\text{backtrack}} H_{\text{backtrack}}$$
  1. **Interaction Hamiltonian ($H_{\text{interactions}}$)**:
     $$H_{\text{interactions}} = \sum_{1 \le i < j \le N, j - i \ge 2} \epsilon(s_i, s_j) \cdot \Delta(\vec{r}_i, \vec{r}_j)$$
     where $\Delta(\vec{r}_i, \vec{r}_j) = 1$ if residues $i$ and $j$ are nearest neighbors on the lattice, and $0$ otherwise.
     - $\epsilon(s_i, s_j)$ is retrieved from the **Miyazawa-Jernigan (MJ)** $20 \times 20$ contact potential matrix derived from high-resolution PDB crystal structures (Miyazawa & Jernigan 1985, 1996). Strong hydrophobic and aromatic pairs yield negative values (e.g., Cys-Cys: $-6.85\ k_B T$, Phe-Phe: $-5.18\ k_B T$, Trp-Trp: $-5.09\ k_B T$, Ile-Ile: $-4.46\ k_B T$).
  2. **Overlap Penalty ($H_{\text{overlap}}$)**:
     $$H_{\text{overlap}} = \sum_{1 \le i < j \le N} \delta(\vec{r}_i, \vec{r}_j)$$
     penalizing any sterically invalid self-intersecting configuration.
  3. **Backtrack Penalty ($H_{\text{backtrack}}$)**:
     $$H_{\text{backtrack}} = \sum_{k=1}^{N-2} \delta(\vec{v}_{k+1}, -\vec{v}_k)$$
     penalizing non-physical immediate $180^\circ$ reverse folding along the chain.
  4. **Rigorous Penalty Multiplier Bound**:
     To ensure that no overlapping conformation can ever have a lower total energy than a non-overlapping valid conformation:
     $$\lambda_{\text{overlap}} > \sum_{1 \le i < j \le N} |\min(0, \epsilon(s_i, s_j))|$$
- **Class `FoldingHamiltonian`**:
  - `__init__(protein: LatticeProtein, lambda_overlap: Optional[float] = None, lambda_backtrack: float = 50.0)`
  - `evaluate_conformation(turns: Sequence[int]) -> float`: Evaluates total scalar energy for a given conformation.
  - `to_qubo() -> Tuple[np.ndarray, float]`: Converts the turn-encoded folding Hamiltonian into a Quadratic Unconstrained Binary Optimization (QUBO) matrix $Q$ with offset $c$, such that $E(\mathbf{q}) = \mathbf{q}^T Q \mathbf{q} + c$. Higher-order terms are quadratized using Rosenberg/Freedman penalty ancillae.

### 3.4 `quantum_folding/circuit.py`
- **Purpose**: Synthesizes discrete-variable linear optical circuits using the **Perceval** SDK (`pcvl`), mapping quantum folding decision states to optical unitaries.
- **Architectural Paradigms**:
  1. **Continuous Fock-State Variational Optical Circuit (VQA on LOQC)**:
     - Implements a universal linear interferometer on $M$ optical modes parameterized by a vector of beam splitter phases $\boldsymbol{\theta} = (\theta_1, \dots, \theta_K)$ and phase shifters $\boldsymbol{\phi} = (\phi_1, \dots, \phi_K)$.
     - Decomposes arbitrary unitary operations $U(\boldsymbol{\theta}) \in \mathcal{U}(M)$ using the **Clements et al. (2016)** triangular mesh architecture:
       ```python
       circuit = pcvl.GenericInterferometer(m, pcvl.Circuit.mzi_clements)
       ```
     - Input State: Pure Fock state $|1, 1, \dots, 1, 0, \dots\rangle$ with $N_{\text{photons}} \le M$.
     - Output Sampling: Measurement yields Fock state configurations $|\vec{n}\rangle = |n_1, n_2, \dots, n_M\rangle$. Each output Fock state is mapped to a lattice turn vector.
  2. **Dual-Rail Qubit Mapping (Gate-Based Emulation on LOQC)**:
     - Logical $|0\rangle_L = |1, 0\rangle$ and Logical $|1\rangle_L = |0, 1\rangle$ across 2 spatial modes.
     - Single-qubit rotations implemented via parameterized MZIs.
     - Two-qubit entangling operations (CZ/CNOT) synthesized using Knill-Laflamme-Milburn (KLM) heralded optical gates with post-selection.
- **Class `PercevalCircuitBuilder`**:
  - `build_variational_interferometer(modes: int, depth: int) -> Tuple[pcvl.Circuit, List[pcvl.Parameter]]`: Constructs a parameterized MZI mesh.
  - `create_input_state(photons: int, modes: int) -> pcvl.BasicState`: Returns initialized input state.
  - `map_fock_to_turns(state: pcvl.BasicState, num_turns: int) -> List[int]`: Deterministic mapping from optical detector click patterns to discrete peptide turn choices.

### 3.5 `quantum_folding/quandela_backend.py`
- **Purpose**: Integrates the proteinAI pipeline with the Quandela Cloud RemoteProcessor API (`api.cloud.quandela.com`), managing authentication, quota isolation, execution telemetry, and hardware fallbacks.
- **Key Specifications Mined & Verified**:
  - **Base URL**: `https://api.cloud.quandela.com`
  - **Authentication Token Format**: The literal prefix `_T_` is **mandatory**. The Bearer header must be:
    `Authorization: Bearer _T_eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9...`
    (Omitting `_T_` triggers HTTP 401 Unauthorized `{"error": "Invalid token"}`).
  - **Cloud Platforms**:
    * `sim:slos`: Strong Linear Optical Simulator hosted on OVHcloud. Supports 2–12 photons, accepts multi-photon inputs, supports `probs`, `samples`, and `sample_count`.
    * `qpu:ascella`: Physical 12-mode photonic processor based on quantum dot single-photon source. Accepts max 6 single photons, zero multi-photon support (`support_multi_photon: False`), single-photon emission clock 13.33 MHz, HOM visibility 88.98%.
  - **Mandatory SDK Constraints (Perceval 1.2.4)**:
    * `processor.min_detected_photons_filter(n)`: Mandatory post-selection threshold. Execution throws `ValueError` if omitted.
    * `Sampler(processor, max_shots_per_call=1000)`: Mandatory execution configuration for remote backends. Throws `RuntimeError` if omitted.
    * In Perceval 1.2+, `sampler.sample_count` and `sampler.probs` are properties returning job objects; shots must be passed to `.execute_sync(shots)` or `.execute_async()`.
  - **Concurrency & Quota Handling**:
    * Token quota: `max_concurrent_jobs: 1`, `max_jobs_in_queue: 1`.
    * Sequential dispatch is enforced using a threading lock and exponential backoff retry.
  - **Automatic Result Decompression**:
    * Quandela Cloud returns results encoded in zlib compressed payloads (e.g. `:PCVL:zip:eJyzCnAO8...`). The client automatically decodes these into Python dictionaries.
- **Class `QuandelaCloudClient`**:
  - `__init__(token: str, platform: str = "sim:slos")`
  - `check_availability() -> Dict[str, Any]`: Queries `/api/jobs/availability` to ensure zero concurrent jobs are pending.
  - `check_platform_status(platform: str) -> str`: Returns `"available"`, `"maintenance"`, or `"unreachable"`.
  - `submit_job_sync(circuit: pcvl.Circuit, input_state: pcvl.BasicState, shots: int = 100) -> Dict[str, Any]`: Submits and blocks until completion.
  - `submit_job_async(circuit: pcvl.Circuit, input_state: pcvl.BasicState) -> str`: Submits asynchronous job and returns job UUID.
  - `poll_job_status(job_id: str) -> str`: Returns status (`"completed"`, `"running"`, `"waiting"`, `"failed"`).
  - `get_job_results(job_id: str) -> Dict[str, Any]`: Retrieves and decodes output Fock state distribution.

### 3.6 `quantum_folding/vqe.py`
- **Purpose**: Executes the Variational Photonic Eigensolver loop, optimizing optical interferometer parameters to discover ground-state minimum energy conformations.
- **Objective Formulations**:
  1. **Standard Expectation Energy**:
     $$\langle H \rangle_{\boldsymbol{\theta}} = \sum_{\vec{n}} P(\vec{n}; \boldsymbol{\theta}) E(\vec{n})$$
     where $P(\vec{n}; \boldsymbol{\theta})$ is the probability of detecting Fock state $\vec{n}$, and $E(\vec{n})$ is the energy computed by `FoldingHamiltonian`.
  2. **Conditional Value-at-Risk (CVaR-VQE)**:
     Following the methodology of Robert et al. (IBM / *npj Quantum Information*, 2021), standard expectation value averaging is susceptible to high sample variance and noise from high-energy excited states. CVaR-VQE evaluates only the $\alpha$-percentile lowest energy samples ($\alpha \in (0, 1]$):
     $$\text{CVaR}_\alpha(\boldsymbol{\theta}) = \frac{1}{\lceil \alpha K \rceil} \sum_{k=1}^{\lceil \alpha K \rceil} E_{(k)}$$
     where $E_{(1)} \le E_{(2)} \le \dots \le E_{(K)}$ are the sorted sampled energies from $K$ shots. Setting $\alpha = 0.20$ guarantees convergence directly toward the deepest potential well (the native fold), suppressing thermal and optical noise.
- **Class `VQESolver`**:
  - `__init__(hamiltonian: FoldingHamiltonian, circuit_builder: PercevalCircuitBuilder, backend: QuandelaCloudClient, alpha: float = 0.20, optimizer: str = "COBYLA")`
  - `optimize(max_iter: int = 50) -> VQEResult`: Executes classical parameter updates via `scipy.optimize.minimize`, tracking optimization trajectory, shot counts, cloud job IDs, and minimum energy turns.

### 3.7 `quantum_folding/pipeline.py`
- **Purpose**: Orchestrates batch evaluation across candidate pools, coordinates epitope loop extraction, computes quantum folding scores, and performs hybrid multi-objective ranking.
- **Epitope Loop Windowing**:
  Full 85–100 residue proteins cannot be mapped into 12 photonic modes. The pipeline extracts critical interaction loops:
  - **TL1A Minibinders**: The 10-residue flexible loop engaging the trimer groove (analogous to the TL1A receptor-binding loop Leu149–Gly160).
  - **IL-7Rα Minibinders**: The 7-residue core binding motif (homologous to `KKFLLIG`, Lys115–Gly121 of IL-7Rα).
- **Hybrid Multi-Objective Scoring Formula**:
  $$\text{Score}_{\text{hybrid}} = w_1 \cdot \text{ipSAE} + w_2 \cdot S_{\text{quantum}}$$
  where:
  - $w_1 = 0.60$ (weight for classical interface accessible solvation energy).
  - $w_2 = 0.40$ (weight for quantum folding stability).
  - $S_{\text{quantum}} = \frac{E_{\text{unfolded}} - E_{\text{ground}}}{E_{\text{unfolded}} - E_{\text{min\_theoretical}}} \in [0.0, 1.0]$.
  - Filter: $\text{Novelty PASS}$ ($< 80\%$ coverage against known binders).
- **Class `QuantumFoldingPipeline`**:
  - `evaluate_candidate_pool(csv_path: str, target: str, max_candidates: int = 50) -> pd.DataFrame`
  - `rank_and_export_submissions(evaluated_df: pd.DataFrame, target: str, output_path: str, top_k: int = 20)`

---

## 4. New CLI Tool Specification: `tools/submit_sjtu.py`

The competition submission client is implemented as a standalone, production-grade CLI script located at `/home/khang/qaaa-lab/proteinAI/tools/submit_sjtu.py`. It automates authentication, target schema discovery, client-side pre-validation, batch uploading, and submission tracking on `https://aiprotein.sjtu.edu.cn`.

### 4.1 Authentication & Session Management
- **Endpoint**: `POST /api/auth/login`
- **Request Payload**:
  ```json
  {
    "email": "khangmaiduc@gmail.com",
    "password": "...",
    "portal": "contestant"
  }
  ```
  *Crucial Rule*: `portal: "contestant"` is strictly required; requests omitting `portal` fail or default to administrative scopes.
- **Token Handling**:
  - Captures JWT `access_token` and `refresh_token`.
  - Attaches `Authorization: Bearer <access_token>` to all protected requests.
  - Automatically intercepts HTTP 401 responses and triggers `POST /api/auth/refresh` using the `refresh_token`.
- **User & Team Identity**:
  - User: Mai Duc Khang (ID: `07fe1517-b035-4805-a505-9a3fc68547fd`)
  - Team: `QAAA Lab` (ID: `66d2100d-6eb1-42ba-9dcd-0589b0654046`, Track: `open`, Status: `approved`)
  - Role: `captain` (authorized by platform logic `dr = l.members.some(z => z.member_role === "captain")` to submit sequences).

### 4.2 Target Track Discovery
- **Endpoint**: `GET /api/submissions/targets`
- **Active Minibinder Target Specs**:

| Field | TL1A (Track 5) | IL-7Rα (Track 7) |
|---|---|---|
| Target UUID | `ea930a01-20a9-4bed-b70f-e658bd9ba3ab` | `1347d0e0-7d3d-4aa8-bcf1-c2dab9d8f1df` |
| Code | `TL1A` | `IL7RA` |
| Modality | `minibinder` | `minibinder` |
| Length Range | `[40, 100]` amino acids | `[40, 100]` amino acids |
| Max Sequences | 20 sequences | 20 sequences |
| Deadline | `2026-11-30T23:59:59+08:00` | `2026-11-30T23:59:59+08:00` |
| Reference PDB | 2RE9 (homotrimer, chains A, B, C) | 3DI3 (chain B receptor) / 7OPB |

### 4.3 Strict Client-Side Pre-Validation Rules
Before transmitting any payload to the remote server, `tools/submit_sjtu.py` executes local pre-validation replicating the platform's frontend parser:
1. **Header Check**: First row must contain `name` and `sequence` columns (case-insensitive).
2. **Sequence Count**: $1 \le N \le 20$ rows. If $N > 20$, reject immediately with `"最多提交 20 条序列（含表头共 21 行）"`.
3. **Identifier Rules**: Length $1 \le \text{len}(\text{name}) \le 64$. All names must be unique within the submission batch.
4. **Canonical Amino Acids**: Sequence must consist exclusively of the 20 standard amino acids: `ACDEFGHIKLMNPQRSTVWY`. Rejects `B, Z, X, J, U, O` or non-alphabetic characters.
5. **Length Enforcement**: For modality `minibinder`, length must satisfy $40 \le L \le 100$ residues.
6. **Mature Sequence Only**: Rejects sequences containing signal peptides, purification tags (e.g. `HHHHHH`, `DYKDDDDK`), linkers (`GGGGS`), or Fc domains.

### 4.4 Submission API Schema & State Transitions
- **Endpoint**: `POST /api/submissions/csv`
- **Request Payload**:
  ```json
  {
    "target_id": "ea930a01-20a9-4bed-b70f-e658bd9ba3ab",
    "action": "save",
    "agreed": true,
    "abstract": "Hybrid quantum-classical de novo minibinders designed via RFdiffusion, ESMFold2-ipSAE screening, and Quandela LOQC VQE folding optimization.",
    "abstract_name": "",
    "abstract_file": "",
    "note": "Initial submission batch evaluated via Perceval LOQC",
    "original_name": "TL1A_QAAA_Lab_submission.csv",
    "csv": "name,sequence\nTL1A_Q_0001,...\n"
  }
  ```
- **Action Semantics & Safety Safeguard**:
  - `action: "save"`: Uploads or updates a draft submission without locking the track. Can be modified repeatedly.
  - `action: "confirm"`: Formally commits the batch for official in-silico screening and wet-lab evaluation.
  - **Safety Mechanism in CLI**: `tools/submit_sjtu.py` defaults strictly to `action: "save"` unless an explicit `--confirm` flag is provided by the operator.
  - **Re-upload Requirement**: When re-uploading to a target where a submission already exists, the platform enforces a non-empty `note` (up to 2,000 characters) explaining the changes.

### 4.5 Submission Tracking & Query Endpoints
- `GET /api/submissions/latest?target_id=<uuid>`: Queries metadata of current batch (`count`, `source_file`, `created_at`, `status`).
- `GET /api/submissions?grouped=1`: Returns status grouping (`submitted`, `shortlisted`, `tested`, `rejected`).
- `GET /api/leaderboard`: Monitors published wet-lab binding affinity ($K_D$ in nM) once experimental testing concludes.

### 4.6 CLI Command-Line Syntax

```bash
# 1. Test platform authentication and check user/team status
python3 tools/submit_sjtu.py --check-auth

# 2. List active targets, target IDs, and submission deadlines
python3 tools/submit_sjtu.py --list-targets

# 3. Dry-run local validation of candidate CSV (no network submission)
python3 tools/submit_sjtu.py --validate --target TL1A --file data/submission_ready/TL1A_submission.csv

# 4. Upload draft submission (action: save)
python3 tools/submit_sjtu.py --submit --target TL1A --file data/submission_ready/TL1A_submission.csv --note "M3 draft upload"

# 5. Final confirmed submission (action: confirm)
python3 tools/submit_sjtu.py --submit --target TL1A --file data/submission_ready/TL1A_submission.csv --confirm --note "Final confirmed candidates"

# 6. Check submission status for a target
python3 tools/submit_sjtu.py --status --target TL1A
```

---

## 5. Operational Deployment & Hardware Execution Model

### 5.1 Quandela Cloud Platform Specifications

Quandela Cloud provides access to discrete-variable linear optical quantum processing backends:

| Metric / Parameter | `sim:slos` (Cloud Simulator) | `qpu:ascella` (Physical QPU) |
|---|---|---|
| Platform Architecture | Strong Linear Optical Simulator (SLOS) | Reconfigurable Silicon Nitride / Silica Waveguide |
| Hosting Infrastructure | OVHcloud High-Performance Compute | Quandela Quantum Hardware Laboratory |
| Spatial Modes ($M$) | Scalable (up to 30 modes) | Exactly 12 optical modes |
| Photon Count ($N$) | 2 to 12 photons | Max 6 single photons |
| Multi-Photon States | Supported ($|2, 0\rangle, |0, 2\rangle$) | **Not Supported** (`support_multi_photon: False`) |
| Photon Source | Ideal indistinguishable Fock states | InGaAs/GaAs Quantum Dot in Micropillar |
| Single-Photon Clock | Simulated | 13.33 MHz |
| HOM Interference Visibility | $100\%$ | $88.98\%$ |
| Chip Optical Transmittance | $100\%$ | $5.56\%$ |
| Detection Technology | Direct probability calculation / sampling | Superconducting Nanowire (SNSPD) |
| Operational Status | **Available (200 OK)** | **Maintenance / Bookable** |

### 5.2 Quota Handling, Concurrency & Error Recovery
- **Token Quota**: The provided credentials permit **1 concurrent job** and **1 queued job**. Attempting to dispatch parallel asynchronous tasks results in immediate job rejection or queue stalls.
- **Sequential Dispatch Queue**: The `QuandelaCloudClient` maintains an internal synchronization barrier. When running VQE optimization loops (50+ iterations), each iteration awaits completion of the remote job before updating parameters and submitting the subsequent circuit.
- **Failover Strategy**:
  1. The client queries `GET /api/platform/{target_platform}` before job dispatch.
  2. If `qpu:ascella` returns status `maintenance` or `unreachable`, the client automatically routes execution to `sim:slos`.
  3. A warning log is recorded in `logs/quantum_execution.log`, ensuring uninterrupted pipeline execution.
- **Exponential Backoff**: Transient network errors or HTTP 429 rate limit triggers activate an exponential backoff routine ($t_{\text{wait}} = 2^k \text{ seconds}$, $k \in \{1, 2, 3\}$).

### 5.3 Competition Target Epitope Integration

To maximize biological efficacy against the competition targets, the quantum folding engine focuses its computational resources on critical binding loops identified from reference crystallographic structures:

#### Target 1: TL1A Homotrimer (PDB: 2RE9)
- **Biological Role**: TL1A (TNFSF15) binds to Death Receptor 3 (DR3 / TNFRSF25) and decoy receptor DcR3 (TNFRSF6B) via three symmetric receptor-binding grooves located at the subunit interfaces.
- **Critical Epitope Window**: Residues Leu149 to Gly160 form a flexible loop protruding into the receptor-binding groove.
- **Minibinder Design Strategy**: De novo minibinders (85–100 aa) present an alpha-helical scaffold positioning hydrophobic sidechains into the TL1A groove, disrupting homotrimer engagement. The quantum folding engine simulates the folding conformation of the minibinder's binding loop to verify that the binding-competent state represents the global free energy minimum.

#### Target 2: IL-7Rα Receptor (PDB: 3DI3 & 7OPB)
- **Biological Role**: Interleukin-7 receptor alpha subunit binds IL-7 through its extracellular domain (Glu21 to Asp239), initiating JAK/STAT signaling.
- **Critical Epitope Window**: The primary binding hotspot comprises the hydrophobic loop `KKFLLIG` (Lys115 to Gly121 of IL-7Rα), which anchors into the IL-7 cytokine binding cleft.
- **Minibinder Design Strategy**: Minibinders mimic IL-7 contact residues, forming a stable hydrophobic cluster against `KKFLLIG`. The quantum engine encodes this 7-residue interaction turn into an optical interferometer, confirming low-energy contact cooperativity.

---

## 6. Verification and Acceptance Matrix

| Feature / Contract | Verification Command | Success Acceptance Criteria |
|---|---|---|
| **Quandela Cloud Connectivity** | `python3 -c "import perceval; ..."` | HTTP 200 on `/api/jobs/availability`; status returns valid quotas |
| **Quandela Execution** | Synchronous execution on `sim:slos` | Job returns `BSCount` distribution with HOM interference bunching |
| **Lattice & Turn Math** | Unit test `tests/test_quantum_lattice.py` | 100% self-avoiding turn calculations, correct coordinate reconstruction |
| **Hamiltonian & MJ Potential** | Unit test `tests/test_quantum_lattice.py` | $\lambda_{\text{overlap}} > \sum \|\epsilon(s_i, s_j)\|$, correct MJ values |
| **Perceval LOQC Synthesis** | Unit test `tests/test_quantum_circuit.py` | Universal Clements MZI synthesis passes unitary checks |
| **SJTU Authentication** | `python3 tools/submit_sjtu.py --check-auth` | User "Mai Duc Khang", Team "QAAA Lab", Role "captain" verified |
| **SJTU Target Discovery** | `python3 tools/submit_sjtu.py --list-targets` | 7 targets retrieved; TL1A & IL-7Ra UUIDs confirmed active |
| **CSV Validation Pre-check** | `python3 tools/submit_sjtu.py --validate ...` | Rejects non-canonical AAs, wrong lengths; passes valid 85-100 aa |
| **Hybrid Pipeline Ranking** | `python3 scripts/run_quantum_screening.py` | Generates top 20 candidate CSVs for TL1A and IL-7Ra |

---

*End of Architecture Specification.*
