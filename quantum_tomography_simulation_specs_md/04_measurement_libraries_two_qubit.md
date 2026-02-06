# Two-qubit measurement libraries (Libraries A and B)

This file specifies the two two-qubit measurement libraries used in the simulation.

System: \(q=2\) qubits ⇒ \(N=4\), \(d=15\).

We assume the two-qubit Bloch basis is the scaled Pauli-product basis
\(\sigma_{(\alpha,\beta)}=\tfrac12(\sigma_\alpha\otimes\sigma_\beta)\) for \((\alpha,\beta)\neq(I,I)\)
as specified in `02_state_basis_one_two_qubit.md`.

---

## Shared helpers

### One-qubit primitives
Use `build_pauli_basis_1q()` to obtain:
- \(I, X, Y, Z\)

### Kronecker product
Use `kronecker(A,B)`.

### Computational basis vectors |00>,|01>,|10>,|11>
Represent as standard basis vectors in \(\mathbb C^4\):
- \(|00\rangle = (1,0,0,0)^\top\)
- \(|01\rangle = (0,1,0,0)^\top\)
- \(|10\rangle = (0,0,1,0)^\top\)
- \(|11\rangle = (0,0,0,1)^\top\)

### Rank-1 projector from a state vector
```r
proj_rank1 <- function(psi) {
  psi <- matrix(psi, ncol = 1)
  psi %*% Conj(t(psi))
}
```

---

## Library A (binary-optimal): 15 two-qubit Pauli-parity PVMs

### Settings
Define the 15 non-identity two-qubit Pauli operators:
\[
\mathcal P_2 := \{\sigma_\alpha\otimes\sigma_\beta:\alpha,\beta\in\{I,X,Y,Z\},(\alpha,\beta)\neq(I,I)\}.
\]
So \(|\mathcal P_2|=15\).

### Measurement at setting P
Each \(P\in\mathcal P_2\) satisfies \(P^2=I_4\) and has spectrum \(\{+1,+1,-1,-1\}\).
The associated binary PVM is:
\[
Q_{P,+}=\frac{I_4+P}{2},\qquad Q_{P,-}=\frac{I_4-P}{2}.
\]

### Builder
We define a canonical ordering for the 15 settings; it should match the basis ordering where possible.

Suggested ordering (same as in `build_pauli_product_basis_2q()`):
`IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ`

Builder:

```r
build_library_2q_A <- function(order = c("IX","IY","IZ","XI","XX","XY","XZ","YI","YX","YY","YZ","ZI","ZX","ZY","ZZ")) {
  pauli <- build_pauli_basis_1q()
  single <- list(I=pauli$I, X=pauli$sigmas$X, Y=pauli$sigmas$Y, Z=pauli$sigmas$Z)
  kron <- function(A,B) kronecker(A,B)
  I4 <- diag(4)

  P_list <- list()
  Q_list <- list()
  for (lab in order) {
    a <- substr(lab,1,1); b <- substr(lab,2,2)
    P <- kron(single[[a]], single[[b]])
    P_list[[lab]] <- P
    Qp <- 0.5*(I4 + P)
    Qm <- 0.5*(I4 - P)
    Q_list[[lab]] <- list(Qp, Qm)
  }
  build_ab_indexing(Q_list, setting_labels = names(Q_list), N = 4, extra = list(P_list = P_list))
}
```

---

## Library B (rank-1 PVM-optimal): 5 MUB bases (4 outcomes each)

For \(N=4\) (prime power), there exist \(N+1=5\) mutually unbiased bases (MUBs).
Measuring all 5 bases is informationally complete and meets the counting minimum \(k=5\) for \(r=4\).

### One-qubit states
Define:
\[
|+\rangle = \frac{|0\rangle+|1\rangle}{\sqrt2},\quad
|-\rangle = \frac{|0\rangle-|1\rangle}{\sqrt2},
\]
\[
|+i\rangle = \frac{|0\rangle+i|1\rangle}{\sqrt2},\quad
|-i\rangle = \frac{|0\rangle-i|1\rangle}{\sqrt2}.
\]

### The 5 explicit bases

Let the computational basis be \(\{|00\rangle,|01\rangle,|10\rangle,|11\rangle\}\).

1. **Basis 1 (Z⊗Z / computational):**
   \[
   B_Z=\{|00\rangle,|01\rangle,|10\rangle,|11\rangle\}.
   \]

2. **Basis 2 (X⊗X product basis):**
   \[
   B_X=\{|++\rangle,|+-\rangle,|-+\rangle,|--\rangle\}.
   \]

3. **Basis 3 (Y⊗Y product basis):**
   \[
   B_Y=\{|+i,+i\rangle,|+i,-i\rangle,|-i,+i\rangle,|-i,-i\rangle\}.
   \]

4. **Basis 4 (entangled MUB):**
   \[
   B_{E1}=\left\{
   \tfrac12(|00\rangle+|01\rangle+i|10\rangle-i|11\rangle),
   \tfrac12(|00\rangle+|01\rangle-i|10\rangle+i|11\rangle),
   \tfrac12(|00\rangle-|01\rangle+i|10\rangle+i|11\rangle),
   \tfrac12(|00\rangle-|01\rangle-i|10\rangle-i|11\rangle)
   \right\}.
   \]

5. **Basis 5 (entangled MUB):**
   \[
   B_{E2}=\left\{
   \tfrac12(|00\rangle+i|01\rangle-|10\rangle+i|11\rangle),
   \tfrac12(|00\rangle-i|01\rangle+|10\rangle+i|11\rangle),
   \tfrac12(|00\rangle+i|01\rangle+|10\rangle-i|11\rangle),
   \tfrac12(|00\rangle-i|01\rangle-|10\rangle-i|11\rangle)
   \right\}.
   \]

These 5 bases are mutually unbiased: for \(a\neq a'\) and \(|\psi\rangle\in B_a,|\phi\rangle\in B_{a'}\),
\(|\langle\psi|\phi\rangle|^2=1/4\).

### Turning bases into a measurement (POVM)
For each basis \(B_a=\{|\psi_{a,0}\rangle,\dots,|\psi_{a,3}\rangle\}\), define rank-1 projectors:
\[
\Pi_{a,b}=|\psi_{a,b}\rangle\langle\psi_{a,b}|,
\quad b\in\{0,1,2,3\}.
\]
The PVM is \(\{\Pi_{a,b}\}_{b=0}^3\).

### Optional: assign eigenvalues to form an observable
You can define \(O_a=\sum_{b=0}^3\lambda_b\Pi_{a,b}\) with distinct eigenvalues (e.g. 3,1,-1,-3).
This is useful only for physics interpretation; the tomography code only needs the projectors.

### Builder
```r
build_library_2q_B <- function() {
  pauli <- build_pauli_basis_1q()
  # |0>, |1>
  ket0 <- matrix(c(1,0), ncol=1) + 0i
  ket1 <- matrix(c(0,1), ncol=1) + 0i

  ketp  <- (ket0 + ket1)/sqrt(2)
  ketm  <- (ket0 - ket1)/sqrt(2)
  ketpi <- (ket0 + 1i*ket1)/sqrt(2)
  ketmi <- (ket0 - 1i*ket1)/sqrt(2)

  kron <- function(A,B) kronecker(A,B)

  # computational basis in C^4
  ket00 <- kron(ket0, ket0)
  ket01 <- kron(ket0, ket1)
  ket10 <- kron(ket1, ket0)
  ket11 <- kron(ket1, ket1)

  BZ  <- list(ket00, ket01, ket10, ket11)
  BX  <- list(kron(ketp,ketp), kron(ketp,ketm), kron(ketm,ketp), kron(ketm,ketm))
  BY  <- list(kron(ketpi,ketpi), kron(ketpi,ketmi), kron(ketmi,ketpi), kron(ketmi,ketmi))

  BE1 <- list(
    0.5*(ket00 + ket01 + 1i*ket10 - 1i*ket11),
    0.5*(ket00 + ket01 - 1i*ket10 + 1i*ket11),
    0.5*(ket00 - ket01 + 1i*ket10 + 1i*ket11),
    0.5*(ket00 - ket01 - 1i*ket10 - 1i*ket11)
  )

  BE2 <- list(
    0.5*(ket00 + 1i*ket01 - ket10 + 1i*ket11),
    0.5*(ket00 - 1i*ket01 + ket10 + 1i*ket11),
    0.5*(ket00 + 1i*ket01 + ket10 - 1i*ket11),
    0.5*(ket00 - 1i*ket01 - ket10 - 1i*ket11)
  )

  bases <- list(BZ=BZ, BX=BX, BY=BY, BE1=BE1, BE2=BE2)

  # Q_list: each setting a is a list of 4 projectors
  Q_list <- lapply(bases, function(B) lapply(B, proj_rank1))

  build_ab_indexing(Q_list, setting_labels = names(Q_list), N = 4, extra = list(bases = bases))
}
```

---

## Informational completeness (counting sanity check)

- Library A: \(k=15\) settings, each has \(r=2\) outcomes ⇒ total constraints \(15(2-1)=15=d\).
- Library B: \(k=5\) settings, each has \(r=4\) outcomes ⇒ total constraints \(5(4-1)=15=d\).

Both achieve the counting minimum for their outcome structures.

