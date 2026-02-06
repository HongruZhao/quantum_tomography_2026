# Bloch bases and state parametrizations (one qubit and two qubits)

This file specifies the Bloch-coordinate representation used in the simulation study.

## 1) Common Bloch model

For Hilbert-space dimension \(N\), let \(d=N^2-1\).
Fix a traceless Hermitian basis \(\{\sigma_j\}_{j=1}^d\) of \(\mathrm{Herm}_0(\mathbb C^N)\).

We represent a state by its Bloch vector \(\theta\in\mathbb R^d\) via
\[
\rho(\theta) = \frac{I}{N} + \frac12\sum_{j=1}^d \theta_j \sigma_j.
\]

In R, this is already implemented by `rho_of_theta(theta, sigmas, N)`.

### Basis normalization
The PDF’s SU(N) construction yields \(\mathrm{tr}(\sigma_j\sigma_k)=2\delta_{jk}\).
For Pauli-product bases we may use a different scaling, but then:
- the mapping \(\theta \mapsto \rho(\theta)\) is still correct,
- but Frobenius / observable metric constants change by a scalar factor (usually immaterial for argmin/argmax rules).

To avoid confusion, we explicitly define the basis scaling we will use for each system.

---

## 2) One qubit (N = 2)

### Pauli matrices
Use the standard Pauli matrices:
\[
X=\begin{bmatrix}0&1\\1&0\end{bmatrix},\quad
Y=\begin{bmatrix}0&-i\\ i&0\end{bmatrix},\quad
Z=\begin{bmatrix}1&0\\0&-1\end{bmatrix}.
\]
They satisfy \(\mathrm{tr}(X^2)=\mathrm{tr}(Y^2)=\mathrm{tr}(Z^2)=2\) and are orthogonal.

### Parametrization
With `sigmas = list(X, Y, Z)`, the Bloch map becomes:
\[
\rho(\theta)=\frac12\left(I + \theta_1 X + \theta_2 Y + \theta_3 Z\right).
\]
Full rank corresponds to \(\|\theta\|_2<1\).

### Required R builder

Implement a function that returns these matrices:

```r
build_pauli_basis_1q <- function() {
  I2 <- diag(2)
  X <- matrix(c(0,1,1,0), 2,2) + 0i
  Y <- matrix(c(0,-1i,1i,0), 2,2)
  Z <- matrix(c(1,0,0,-1), 2,2) + 0i
  list(N = 2, sigmas = list(X=X, Y=Y, Z=Z), I = I2)
}
```

---

## 3) Two qubits (N = 4)

### Pauli-product operators
Define the 16 two-qubit Pauli products:
\[
\sigma_\alpha\otimes\sigma_\beta,\quad \alpha,\beta\in\{I,X,Y,Z\},
\]
excluding \(I\otimes I\) to obtain \(d=15\) traceless operators.

### Scaling choice consistent with the PDF simulation section
The PDF writes:
\[
\rho(\theta)=\frac14 I\otimes I + \frac14 \sum_{(\alpha,\beta)\neq(I,I)} \theta_{\alpha\beta}\,\sigma_\alpha\otimes\sigma_\beta.
\]

To match the existing `rho_of_theta` convention
\(\rho(\theta)=I/N + \tfrac12\sum_j \theta_j\sigma_j\),
we choose the Bloch basis elements as:
\[
\sigma_{(\alpha,\beta)} := \frac12(\sigma_\alpha\otimes\sigma_\beta),\quad (\alpha,\beta)\neq(I,I).
\]
Then:
\[
\frac12\theta_{(\alpha,\beta)}\,\sigma_{(\alpha,\beta)}
= \frac12 \theta_{(\alpha,\beta)} \cdot \frac12(\sigma_\alpha\otimes\sigma_\beta)
= \frac14\theta_{(\alpha,\beta)}(\sigma_\alpha\otimes\sigma_\beta),
\]
so `rho_of_theta` reproduces the PDF’s coefficient \(1/4\).

### Required R builder

```r
build_pauli_product_basis_2q <- function(order = c("IX","IY","IZ","XI","XX","XY","XZ","YI","YX","YY","YZ","ZI","ZX","ZY","ZZ")) {
  pauli <- build_pauli_basis_1q()
  I <- pauli$I; X <- pauli$sigmas$X; Y <- pauli$sigmas$Y; Z <- pauli$sigmas$Z

  # helper: kron
  kron <- function(A,B) kronecker(A,B)

  # map label -> matrix
  single <- list(I=I, X=X, Y=Y, Z=Z)

  # build the 15 traceless products with scaling 1/2
  sigmas <- list()
  for (lab in order) {
    a <- substr(lab,1,1); b <- substr(lab,2,2)
    sigmas[[lab]] <- 0.5 * kron(single[[a]], single[[b]])
  }
  list(N = 4, sigmas = sigmas, I = diag(4))
}
```

### Notes on ordering
- We intentionally exclude `"II"`.
- The chosen order defines the coordinate convention for \(\theta\).
- All downstream code must treat this ordering as canonical when plotting or comparing estimates.

---

## 4) Converting between rho and theta

### `theta_from_rho(rho, sigmas)` (already in CVXR1121.R)
Computes \(\theta_j = \mathrm{Re}\,\mathrm{tr}(\rho\sigma_j)\).
This is consistent with the basis scaling above because:
\[
\mathrm{tr}(\rho\sigma_j) = \mathrm{tr}\Big(\big(\tfrac{I}{N}+\tfrac12\sum_k\theta_k\sigma_k\big)\sigma_j\Big)
= \tfrac12\sum_k \theta_k\,\mathrm{tr}(\sigma_k\sigma_j),
\]
so if the basis is orthogonal but not necessarily normalized to 2, `theta_from_rho` returns a scaled coordinate.
Therefore:

**Project requirement:** Always use `rho_of_theta` and `theta_from_rho` with the *same* `sigmas` definition within a run.

If you need a strictly orthonormal coordinate system (e.g., for GF = 1/2 I), implement an optional basis
normalization step and transform G and FI by congruence accordingly.

