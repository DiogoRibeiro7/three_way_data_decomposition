# Tucker3 centroid-mode extension

## Motivation

The validated S3 baseline uses a Tucker2 structure for the component means. Rocci, Vichi, and Ranalli explicitly note that a Tucker3 mean model could additionally reduce the centroid/group mode when the number of mixture components is large.

The extension is

\[
\mu_{gjk}
=
\mu_{jk}
+
\sum_{p=1}^{P}
\sum_{q=1}^{Q}
\sum_{r=1}^{R}
a_{gp} b_{jq} c_{kr} \eta_{pqr}.
\]

The new matrix \(A \in \mathbb{R}^{G\times P}\) reduces the centroid mode from \(G\) component means to \(P\) centroid components.

## Identification used in this implementation

Let \(\pi\) denote the vector of mixture probabilities. The decomposition separates an explicit grand mean from group deviations by imposing

\[
\pi^\top A = 0.
\]

Therefore

\[
\sum_{g=1}^{G} \pi_g \mu_g = \mu.
\]

The current model layer additionally represents \(A\) by an orthonormal basis. As with any Tucker model, basis rotations are not individually identifiable; the estimable object is the corresponding multilinear subspace and reconstructed mean tensor.

## Effective parameter count

The centroid subspace lies in a \((G-1)\)-dimensional contrast space. A rank-\(P\) centroid subspace therefore contributes

\[
P(G-1-P)
\]

degrees of freedom.

For ranks \((P,Q,R)\), the full S3-Tucker3 model count implemented by `scr_tucker3_parameter_count()` is

\[
\begin{aligned}
\nu ={}&
(G-1)
+ JK
+ PQR \\
&+ P(G-1-P)
+ Q(J-Q)
+ R(K-R) \\
&+ \frac{J(J+1)}{2}
+ \frac{K(K+1)}{2}
-1.
\end{aligned}
\]

Setting \(P=G-1\) removes the centroid-subspace term and gives

\[
PQR=(G-1)QR,
\]

so the count reduces exactly to the published S3/Tucker2 parameter count.

## Current scope

This PR implements only the statistical model layer:

- centered centroid-basis construction;
- Tucker3 mean reconstruction;
- dimension and identifiability checks;
- effective parameter counting.

It deliberately does **not** implement an optimizer yet.

The next methodological step is to derive a block update for \(A\) and the core tensor under the S3 likelihood while preserving the covariance-weighted updates for \(B\) and \(C\). That optimizer should be validated first on simulated data where the true centroid rank is known.

## Reference

Rocci, R., Vichi, M. & Ranalli, M. (2025). *Mixture models for simultaneous classification and reduction of three-way data*. Computational Statistics, 40, 469–507. Equation (25) and concluding remarks.
