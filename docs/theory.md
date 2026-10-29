# Theoretical background

This chapter defines the models that the joint fit, the sequential fit, and the two single-peak notebooks optimize. Implementation steps and default numbers are in [Workflows](workflows.md). Column names are in [Outputs](outputs.md).

## Difference density as the data

An omit or difference map stores, at each voxel, a density residual: structure-factor coefficients that show density the current atomic model does not account for. A positive blob at limited resolution is locally similar to a Gaussian. That similarity is the reason to fit one.

The similarity is also the reason to be careful. The observed width mixes the atomic scattering factor, the atomic displacement, the resolution cutoff, the map sharpening, the solvent model, and any neighbors that fall inside the same window. The notebooks state the working rule, and the library follows it. A principal width is a property of the fitted contour. It is not an atomic displacement parameter. The formal Gaussian integral is a property of the amplitude and the covariance. It is not a count of electrons. A center is a position in the map's Cartesian frame. It becomes an atom only after that frame is matched to a structural model and the surrounding positive and negative features have been inspected together.

The library models positive peaks. Amplitudes are bounded below by zero, and a constant background, which may be negative, carries the local baseline of the difference map. Negative peaks are fit by the standalone notebooks, where the amplitude is constrained to the sign of the extremum.

Voxel noise in a crystallographic map is correlated, and neighboring voxels are not independent samples. Every \(R^2\), standard error, and information criterion below is computed as if the residuals were independent. Those scores rank models and expose a fit that hit a bound. They do not, by themselves, prove that two centers are physically resolved. The joint-fit log prints that limitation on purpose.

## Cartesian coordinates on a CCP4 grid

CCP4 stores the grid in column, row, section order. The header records which of those stored axes corresponds to fractional \(x\), \(y\), and \(z\) (`MAPC`, `MAPR`, `MAPS`), where the sampled block starts (`NCSTART`, `NRSTART`, `NSSTART`), and how many intervals span the full cell (`NX`, `NY`, `NZ`). The first array index is therefore not crystallographic \(x\) unless the header says so.

gemmi reads the file with symmetry expansion turned off (`setup=False`). For a stored index \(\mathbf q = (i, j, k)\) the fractional coordinate is

\[
\mathbf f = (\mathbf q_{xyz} + \mathbf s_{xyz}) / \mathbf n_{xyz},
\]

where \(\mathbf q_{xyz}\) is \(\mathbf q\) reordered onto fractional \(x,y,z\), \(\mathbf s_{xyz}\) is the matching start, and \(\mathbf n_{xyz} = (N_X, N_Y, N_Z)\). The Cartesian coordinate in Å is \(\mathbf r = O\mathbf f\), with \(O\) the unit-cell orthogonalization matrix. In code this product is the row vector `fractional @ orth.T`. The same product handles a non-orthogonal cell.

The inverse map, Cartesian Å back to a stored index, is used only when a notebook samples a plane through the grid. The Gaussian itself is always evaluated in Cartesian Å. The sampled block is kept as a finite box. Opposite faces are not identified, and undefined voxels are not invented to fill a unit cell.

The joint fit reports the shortest Cartesian length of the three stored-index steps and uses \(1.8\) times that length as the neighborhood for a local maximum. The sequential fit uses a \(3\times 3\times 3\) index neighborhood for the same idea.

## The Gaussian family

All four workflows are special cases of

\[
\rho(\mathbf r)
  = b
  + \sum_{j=1}^{K}
      A_j
      \exp\!\left(
        -\frac12
        (\mathbf r - \mathbf r_j)^\top
        \Sigma_j^{-1}
        (\mathbf r - \mathbf r_j)
      \right).
\]

\(A_j\) is the amplitude in map units, \(\mathbf r_j\) is the center in Å, \(\Sigma_j\) is a \(3\times 3\) covariance in Å\(^2\), and \(b\) is a constant background in map units. The level set

\[
(\mathbf r - \mathbf r_j)^\top \Sigma_j^{-1} (\mathbf r - \mathbf r_j) = n^2
\]

is the \(n\sigma\) ellipsoid of component \(j\).

### Isotropic component

The library's default shape is one width \(\sigma_j > 0\):

\[
\Sigma_j = \sigma_j^2 I,
\qquad
(\mathbf r - \mathbf r_j)^\top \Sigma_j^{-1} (\mathbf r - \mathbf r_j)
  = \frac{\lVert \mathbf r - \mathbf r_j \rVert^2}{\sigma_j^2}.
\]

The joint model uses one background for the whole sum. A sequential step uses \(K = 1\) and its own background. In either case the isotropic width is bounded to the configured interval, `[min_width, max_width]` or `[sigma_min, sigma_max]`.

### Axis-aligned anisotropic component

`notebooks/xfel_gaussian_fitting_tutorial.ipynb` allows three different widths and keeps the principal axes on the Cartesian directions:

\[
\Sigma = \mathrm{diag}(s_x^2, s_y^2, s_z^2).
\]

The center is a free parameter, the amplitude carries the sign of the extremum, and \(K = 1\). This model is not a switch in the library. A feature whose long axis is tilted relative to \(x,y,z\) is outside its family.

### Fully rotated ellipsoid

A symmetric \(\Sigma\) has six independent entries and is a covariance only when it is positive definite. Both library workflows and `xfel_gaussian_full_covariance_tutorial.ipynb` enforce that by optimizing the Cholesky factor with a positive diagonal,

\[
\Sigma = LL^\top,
\qquad
L =
\begin{pmatrix}
\ell_{00} & 0 & 0 \\
\ell_{10} & \ell_{11} & 0 \\
\ell_{20} & \ell_{21} & \ell_{22}
\end{pmatrix},
\qquad
\ell_{ii} > 0.
\]

The quadratic form is evaluated as a squared Euclidean norm in the Cholesky coordinates, without forming \(\Sigma^{-1}\):

\[
(\mathbf r - \mathbf r_0)^\top \Sigma^{-1} (\mathbf r - \mathbf r_0)
  = \lVert L^{-1}(\mathbf r - \mathbf r_0) \rVert^2.
\]

`scipy.linalg.solve_triangular` applies \(L^{-1}\). The six stored parameters are \((\ell_{00}, \ell_{10}, \ell_{20}, \ell_{11}, \ell_{21}, \ell_{22})\). Off-diagonal \(\ell_{i0}\) are what rotate the ellipsoid. Euler angles are not parameters: they are undefined for a sphere and interchangeable when two widths match. Orientation is reported after the fit by diagonalizing \(\Sigma\).

Any real \(L\) with a positive diagonal gives a positive-definite \(\Sigma\). The optimizer keeps the diagonal at or above the configured minimum width, and it keeps every Cholesky entry inside \(\pm\) the maximum width. Those are bounds on \(L\). A principal width, which is a square root of an eigenvalue of \(\Sigma\), can lie outside that interval.

The initial factor is a multiple of the identity (\(0.8\,I\) in the joint fit, clipped `fit_radius / 3` in the sequential fit, \(0.8\,I\) in the full-covariance notebook). The first residual is therefore the isotropic model. Anisotropy and rotation appear only when they reduce the residual.

### How many numbers each component adds

The joint fit adds one shared background to the count below. A sequential fit is one component plus its own background. The full-covariance notebook fixes the center and fits eight numbers: signed amplitude, six Cholesky entries, and a background. The axis-aligned notebook fits eight numbers as well: signed amplitude, three center coordinates, three widths, and a background.

| Shape | Center | Free parameters per component |
| --- | --- | --- |
| Isotropic | refined | 5: \(A\), three center coordinates, \(\sigma\) |
| Isotropic | fixed | 2: \(A\), \(\sigma\) |
| Full covariance | refined | 10: \(A\), three center coordinates, six Cholesky entries |
| Full covariance | fixed | 7: \(A\), six Cholesky entries |

Fixed centers are removed from the parameter vector. They are not parameters with a tight bound. The BIC penalty, defined below, therefore changes when `fix_centers` changes, and so does the set of models being compared.

## Shape descriptors

After a fit, \(\Sigma = LL^\top\) is diagonalized. Eigenvalues are sorted from largest to smallest, and the principal widths are their positive square roots \(s_1 \ge s_2 \ge s_3\). Each eigenvector is flipped so that its largest-magnitude component is positive. If that produces a left-handed frame, the minor axis is flipped. The result is a right-handed frame. When two widths are nearly equal, the eigenvectors in their plane are weakly determined. The covariance itself remains the stable description; quoting \(\Sigma\) is then more honest than quoting a direction.

The equivalent isotropic width is the geometric mean

\[
\sigma_{\mathrm{eq}} = (s_1 s_2 s_3)^{1/3}.
\]

Along principal axis \(i\), and for the equivalent width,

\[
\mathrm{FWHM} = 2\sqrt{2\ln 2}\; s \approx 2.355\, s.
\]

The formal integral of one component, in map-unit·Å\(^3\), is

\[
\int A \exp\!\left(-\frac12 \mathbf u^\top\Sigma^{-1}\mathbf u\right) \,dV
  = A\,(2\pi)^{3/2}\sqrt{\det\Sigma}
  = A\,(2\pi)^{3/2}\,\ell_{00}\ell_{11}\ell_{22}.
\]

The product of the Cholesky diagonals equals \(\sqrt{\det\Sigma}\) because those diagonals are positive. For an isotropic component the same expression is \(A (2\pi)^{3/2} \sigma^3\). The integral extends over all space. A fit window that clips the tails does not change the formula; it changes how well \(A\) and \(\Sigma\) are known.

The condition number is \(\lambda_{\max}/\lambda_{\min}\) of \(\Sigma\). A value near 1 is a round contour. A very large value is a needle or a pancake, often a peak cut by the edge of a crop or two densities joined inside one window. An isotropic fit has condition number 1, equal principal widths, a diagonal covariance, and Cartesian principal axes up to the sign convention above. The CSV schema is the same for both shapes, so those columns are filled rather than omitted.

A plane through the center with orthonormal directions \(\mathbf u\) and \(\mathbf v\) cuts the \(n\sigma\) ellipsoid in an ellipse. In plane coordinates \((\alpha, \beta)\), with \(P = [\mathbf u\;\mathbf v]\),

\[
\begin{pmatrix}\alpha & \beta\end{pmatrix}
P^\top \Sigma^{-1} P
\begin{pmatrix}\alpha \\ \beta\end{pmatrix}
  = n^2.
\]

The notebooks and `xfel_gaussian_module_calculation.ipynb` draw the \(n = 1\) curve as a solid line and the \(n = 2\) curve as a dashed line on the Cartesian \(xy\), \(xz\), and \(yz\) planes. A circle on one of those planes means the two widths in that plane are equal. A tilted ellipse is the rotation the full covariance adds.

## Nonlinear least squares

Every workflow minimizes the sum of squared voxel residuals

\[
\mathrm{RSS}(\theta) = \sum_{i=1}^{n} \bigl(\rho(\mathbf r_i;\theta) - y_i\bigr)^2
\]

with `scipy.optimize.least_squares`, method `trf` (trust-region reflective), linear loss, and `x_scale='jac'`. The bounds in the previous sections are the feasible box. The residual RMS is

\[
\mathrm{RMSE} = \sqrt{\mathrm{RSS}/n},
\]

and, when the voxel values are not all equal,

\[
R^2 = 1 - \frac{\mathrm{RSS}}{\sum_i (y_i - \bar y)^2}.
\]

\(R^2\) and RMSE are computed inside the fit window only. They inherit the independence assumption above. An \(R^2\) near 1 means the Gaussian-plus-baseline family matches those voxels. It does not mean the feature is one atom.

`n_at_bound` counts parameters whose bound is active, using the solver's active-mask entries of \(\pm 1\). A width or a Cholesky entry sitting on its upper bound usually means the window is too small for that axis, or the blob is not a single component. The parameter is then known only as "at least this large" inside the chosen bounds.

The joint fit tries several starts and keeps the lowest RSS. The sequential fit and the notebooks use one start. Repeating a sequential or notebook fit from a mildly anisotropic start is a useful check that \(\Sigma\) is not an accident of \(L = \sigma I\); the full-covariance notebook lists that check.

## Joint decomposition and model order

Inside one sphere the joint model is the sum at the top of this chapter, with a single \(b\) shared by every component. Components are ordered by the seed scan, from the tallest seed downward, not by the refined amplitude.

### Why a shared background

One constant is enough to remove the local floor of the difference map. A separate baseline per component would trade against the amplitudes wherever the Gaussians overlap, and the overlap is the usual reason to fit them jointly.

### Information criteria

Assume, for ranking only, that the residuals are independent and normal with a common unknown variance. The maximum-likelihood variance on \(n\) voxels is \(\hat\sigma^2 = \mathrm{RSS}/n\), and the profile log-likelihood satisfies

\[
-2\ln L = n\ln(\mathrm{RSS}/n) + n\ln(2\pi) + n.
\]

The Bayesian information criterion and Akaike information criterion then differ across models on the same voxels only through \(\mathrm{RSS}\) and the parameter count \(p\):

\[
\mathrm{BIC} = n\ln(\mathrm{RSS}/n) + p\ln n,
\qquad
\mathrm{AIC} = n\ln(\mathrm{RSS}/n) + 2p.
\]

The code drops the additive constant \(n\ln(2\pi) + n\), and it floors \(\mathrm{RSS}/n\) at \(10^{-14}\) before the logarithm. Differences of BIC are the meaningful quantities. The selected model is the minimum BIC over \(K = 1, \ldots, K_{\max}\). AIC is written and is not used for the choice.

Here \(p\) is the count in the table above, times \(K\), plus one background. Moving from an isotropic component to a full covariance, or from a fixed center to a free center, increases \(p\) and is accepted by BIC only when the drop in \(\mathrm{RSS}\) pays for that increase. Because the criteria use \(\ln(\mathrm{RSS}/n)\), a fit that drives the residual to the floor can post a large negative BIC. Compare values from the same ROI. A BIC from a different radius, or from the synthetic grid, is a different \(n\).

The criterion does not know about correlated map noise, so it can prefer a larger \(K\) than a crystallographer would accept. Pair separations and the residual figure are the companion evidence. The printed note from the joint fit is the operational form of this paragraph: \(K\), BIC, and a covariance standard error do not alone prove physical resolvability.

### Parameter uncertainty

For the selected joint model the program forms the Gauss–Newton approximation

\[
\mathrm{Cov}(\hat\theta)
  \approx
  \hat\sigma^2_{\nu}\,(J^\top J)^{+},
\qquad
\hat\sigma^2_{\nu} = \mathrm{RSS}/(n - p),
\]

where \(J\) is the Jacobian of the residuals at the solution and \(^{+}\) is a pseudoinverse. The exported standard errors are the square roots of the center diagonal entries. They are NaN when the center is fixed, because a fixed center has no column in \(J\). Amplitude and Cholesky uncertainties are computed internally only as part of this vector and are not written. A direction in which \(J\) has no rank contributes nothing reliable; a tiny standard error on a poorly constrained parameter is a property of the pseudoinverse, not a measurement of precision. The approximation again treats voxels as independent.

### Separation of the selected centers

For every unordered pair the joint fit reports the Euclidean distance \(d\), the ratio \(d/(s^{(a)}_{\mathrm{eq}} + s^{(b)}_{\mathrm{eq}})\), and the ratio of \(d\) to the average of the two equivalent full widths at half maximum. A ratio well below 1 means the \(1\sigma\) ellipsoids, replaced here by equivalent spheres, overlap strongly. The fit can still split them, especially when the centers are fixed on two local maxima. The ratios are diagnostics, not a resolution theorem.

### Two-peak synthetic benchmark

The synthetic mode asks a narrower question than a map fit. Two isotropic Gaussians of known width and known separation, plus independent normal noise, are sampled on a regular Cartesian grid. BIC is then allowed to choose \(K = 1\) or \(K = 2\). As the true separation grows, a useful model starts to prefer \(K = 2\) and the recovered distance approaches the true distance. The column \(\mathrm{BIC}(K=1) - \mathrm{BIC}(K=2)\) is positive when the two-component model has the lower BIC.

This is a check that the optimizer and the criterion recover a pair that really is two isotropic blobs. It is not a calibration of omit-map noise, and the truth stays isotropic even if `full_covariance` is on.

## Sequential subtraction

The sequential model is greedy. At step \(k\) the data are the current residual, and the model is a single term of the sum above plus a baseline that belongs to that step only:

\[
\rho_k(\mathbf r) = b_k + A_k \exp\!\left(-\frac12 (\mathbf r-\mathbf r_k)^\top\Sigma_k^{-1}(\mathbf r-\mathbf r_k)\right).
\]

The fit uses voxels inside a sphere about the current seed. After a successful fit with \(A_k > 0\), the residual is replaced by

\[
\rho_{k+1}(\mathbf r) = \rho_k(\mathbf r) - A_k \exp(\cdots)
\]

inside a sphere of radius five times the longest principal width. The baseline \(b_k\) is not removed. Each local baseline would otherwise be subtracted again at the next site and walk the residual downward across the map.

This is a matching-pursuit style of decomposition. It is cheaper than a joint fit over the whole map, and it lets a strong peak reveal a weaker neighbor. The price is order dependence. A broad feature can be carved into several Gaussians. A weak peak under a strong one is invisible until the strong one is removed, and it is never seen if the strong fit absorbed it. Joint refinement of all accepted centers is not a later step.

### Robust threshold

The stopping height is a multiple of a robust scale of the original map. With \(m\) the median of the finite voxels,

\[
s = 1.4826 \cdot \mathrm{median}_i \lvert y_i - m \rvert,
\]

and the search stops when the strongest remaining local maximum of the residual is below `threshold_sigma` \(\times s\). The constant \(1.4826 \approx 1/\Phi^{-1}(3/4)\) converts a median absolute deviation into a standard deviation for normal data. A few very strong peaks move \(s\) less than they move a root-mean-square scale. The scale is still a property of the whole map, including the peaks and any solvent region in the file. On a tight crop that contains only the site, \(s\) is larger and the same `threshold_sigma` is a higher bar. The notebook detectors use the ordinary map standard deviation for an analogous screen, and they say so: a global scale is a screening heuristic.

The threshold is compared with the residual height at the chosen voxel, before the Gaussian is fit. The fitted amplitude can differ from that height because the baseline and the width absorb part of the peak.

### Identifiability at one site

Holding the center fixed answers a different question from refining it. A fixed center reports the shape of the density around the detected voxel, and a grid maximum that sits on a merged feature stays on that voxel. A free center can slide toward the better placement, up to `max_shift` in the sequential fit and anywhere in the joint ROI when centers are free, including onto a neighbor. Use fixed centers when the detected maxima are the objects of interest. Use free centers when the voxel maximum is only a start.

## Reading a result

A joint or sequential table is ready to interpret when the following are true together.

- The fit window surrounds the feature. A principal width near the subtraction or ROI scale, or `n_at_bound` greater than zero, means a bound or the window edge is doing part of the work.
- The residual section, from the module notebook or from the joint diagnostic scatter, has no leftover blob at the same site. A leftover blob is a missing component or a shape the chosen family cannot represent.
- Pair distances, when \(K > 1\), are large enough compared with the equivalent widths that the components are not a split of one contour. Overlapping sites are a reason to use the joint model rather than one ellipsoid.
- The Cartesian frame has been aligned with the structural model before any atom name is attached. Positions may lie outside the conventional unit cell because the map's stored origin was respected. Symmetry copies are a later, explicit choice.
- The PDB names `MN1`–`MN4` and `O1`, `O2`, ... have not been read as elemental identities. They are labels so a graphics program can show the centers.

Repeat the fit at a smaller and a larger radius when a number will be quoted. Agreement of \(\Sigma\) across those windows, and across an isotropic start versus the optimized factor, is the practical statement of stability that the information criterion does not supply.
