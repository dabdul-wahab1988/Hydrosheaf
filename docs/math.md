# Hydrosheaf Math Notes

This document expands the proposal mathematics in a compact, implementation-aligned form.

## Notation

- Graph: $G = (V, E)$ with directed edges $e = (u \to v)$.
- Ion state at node $v$: $x_v \in \mathbb{R}^8$ in fixed order $(\mathrm{Ca}, \mathrm{Mg}, \mathrm{Na}, \mathrm{HCO_3}, \mathrm{Cl}, \mathrm{SO_4}, \mathrm{NO_3}, \mathrm{F})$.
- Auxiliary observations: $y_v = (\mathrm{EC}, \mathrm{TDS}, \mathrm{pH}, \ldots)$.
- Reaction dictionary: $S \in \mathbb{R}^{8 \times m}$ and extents $z \in \mathbb{R}^m$.
- Weights: $W = \mathrm{diag}(w_1, \ldots, w_8)$.

## Units and Conversions

For ion $i$ with molar mass $M_i$ in g/mol:
$$x_i\ (\mathrm{mmol/L}) = \frac{C_i\ (\mathrm{mg/L})}{M_i}$$
Charge equivalents:
$$x_i\ (\mathrm{meq/L}) = |q_i| \cdot x_i\ (\mathrm{mmol/L})$$

## Transport + Reaction Residual

Edge residual for parameters $\theta_e$ and reaction extents $z_e$:
$$r_e(\theta_e, z_e) = x_v - (A(\theta_e) x_u + b(\theta_e) + S z_e)$$

The per-edge objective is:
$$J_e = \|r_e\|_W^2 + \lambda \|z_e\|_1 + P_{\mathrm{EC/TDS}} + P_{\mathrm{iso}}$$
where $\|r\|_W^2 = r^\top W r$ and $\lambda \ge 0$.

## Transport Models

**Evaporation (conservative scaling):**
$$x' = \gamma x,\quad \gamma \ge 1$$
Weighted least squares estimate:
$$\gamma^\star = \frac{x_u^\top W x_v}{x_u^\top W x_u}, \quad \gamma = \max(1, \gamma^\star)$$

**Single-endmember mixing:**
$$x' = (1 - f) x_u + f x_{\mathrm{end}},\quad f \in [0, 1]$$
Let $d = x_{\mathrm{end}} - x_u$. Then:
$$f^\star = \frac{d^\top W (x_v - x_u)}{d^\top W d}, \quad f = \min(1, \max(0, f^\star))$$

Affine form used in code:
$$A(\gamma) = \gamma I,\quad b(\gamma) = 0$$
$$A(f) = (1 - f) I,\quad b(f) = f x_{\mathrm{end}}$$

## Reaction Fit (Sparse, Constrained)

Given transport residual $r = x_v - (A x_u + b)$, solve:
$$\min_{z} \ \|r - S z\|_W^2 + \lambda \|z\|_1 \quad \text{s.t. } \ell \le z \le u$$

Coordinate-descent update (as implemented):
$$\rho_j = s_j^\top W r - \sum_{k \ne j} (s_j^\top W s_k) z_k$$
$$z_j \leftarrow \mathcal{S}\!\left(\frac{\rho_j}{s_j^\top W s_j}, \frac{\lambda}{2 s_j^\top W s_j}\right)$$
with soft-thresholding
$$\mathcal{S}(\alpha, \tau) = \mathrm{sign}(\alpha)\max(|\alpha| - \tau, 0)$$
followed by constraint clipping to $[\ell_j, u_j]$. By default, non-signed reactions are
restricted to $z_j \ge 0$ (dissolution-only).

## PHREEQC Constraint Logic

For mineral $j$ with saturation index $\mathrm{SI}_v^j$ at node $v$ and threshold $\tau$:

- If $\mathrm{SI}_u^j < -\tau$ and $\mathrm{SI}_v^j < -\tau$: dissolution only ($z_j \ge 0$).
- If $\mathrm{SI}_v^j > \tau$: precipitation only ($z_j \le 0$).
- Otherwise: free ($z_j \in \mathbb{R}$).

If PHREEQC is unavailable and `constraints_hard` is true, mineral extents revert to
dissolution-only for unsigned reactions.

## EC/TDS Penalty

Let $s = \mathbf{1}^\top x$ be the ionic sum. Linear predictors:
$$\widehat{\mathrm{EC}} = a_{\mathrm{EC}} s + b_{\mathrm{EC}},\quad \widehat{\mathrm{TDS}} = a_{\mathrm{TDS}} s + b_{\mathrm{TDS}}$$
Penalty term:
$$P_{\mathrm{EC/TDS}} = \eta_{\mathrm{EC}}(\mathrm{EC}_{\mathrm{obs}} - \widehat{\mathrm{EC}})^2
 + \eta_{\mathrm{TDS}}(\mathrm{TDS}_{\mathrm{obs}} - \widehat{\mathrm{TDS}})^2$$

## Isotope Penalty (LMWL + d-excess)
Local meteoric water line (LMWL):
$$\delta^2\mathrm{H} = a + b\,\delta^{18}\mathrm{O}$$
Default (Gibrilla et al., 2022): $\delta^2\mathrm{H} = 8.66 + 7.22\,\delta^{18}\mathrm{O}$.

Define:
$$E = \delta^2\mathrm{H} - (a + b\,\delta^{18}\mathrm{O}) \quad \text{(LMWL residual)}$$
$$d = \delta^2\mathrm{H} - 8\,\delta^{18}\mathrm{O} \quad \text{(d-excess)}$$

Penalty per edge (using upstream $u$ and downstream $v$):
$$P_{\mathrm{iso}} = w_{\mathrm{iso}}\left(P_E + w_d P_d\right)$$
Evaporation hypothesis:
$$P_E = \max(0, |E_u| - |E_v|),\quad P_d = \max(0, d_v - d_u)$$
Mixing hypothesis:
$$P_E = \max(0, |E_v| - |E_u|),\quad P_d = 0$$

## Transport Model Selection and Probabilities
For candidate scores $\{J_c\}$, define:
$$w_c = \exp(-(J_c - J_{\min}))$$
Transport probability for model $m$:
$$p(m) = \frac{\sum_{c \in m} w_c}{\sum_{c} w_c}$$

## Probabilistic Edge Inference (when head is missing)
For wells $i, j$ with estimated heads $\hat h_i, \hat h_j$ and uncertainties $\sigma_i, \sigma_j$:
$$\Delta h_{ij} = \hat h_i - \hat h_j,\quad \sigma_{ij} = \sqrt{\sigma_i^2 + \sigma_j^2}$$
$$p(i \to j) = \Phi\!\left(\frac{\Delta h_{ij}}{\sigma_{ij}}\right)$$

Gradient guard (distance $d_{ij}$ in km):
$$g_{ij} = \frac{|\Delta h_{ij}|}{1000 d_{ij}}$$
If $g_{ij} < g_{\min}$, the probability is pulled toward $0.5$:
$$p \leftarrow 0.5 + (p - 0.5)\frac{g_{ij}}{g_{\min}}$$

Edges are retained when $p \ge p_{\min}$ and then pruned to the $k$ strongest neighbors.

## CoDA & Nitrate Source Discrimination

### Compositional Data Analysis (CoDA)
To provide robust context independent of total concentration (TDS), we transform the 7-ion sub-composition $S^7 = \{Ca, Mg, Na, K, HCO_3, Cl, SO_4\}$ into Euclidean space using Isometric Log-Ratios (ilr).

**Sequential Binary Partition (SBP)**:
$$
ilr_i = \sqrt{\frac{r s}{r+s}} \ln\left( \frac{g(x_+)}{g(x_-)} \right)
$$
where $g(\cdot)$ is the geometric mean, $x_+$ is the group of 'numerator' ions (size $r$), and $x_-$ is the group of 'denominator' ions (size $s$).

### Bayesian Evidence Accumulation
We distinguish Manure ($M$) vs. Fertilizer ($F$) using a log-odds update. Given prior $P(M)$, the initial logit is:
$$
L_0 = \ln\left(\frac{P(M)}{1-P(M)}\right)
$$

For each evidence feature $\phi_k$ (e.g., $z$-scored $\ln(NO_3/Cl)$):
$$
L_{final} = L_0 + \sum_{k} w_k \cdot \text{sign}_k \cdot \phi_k
$$

The posterior probability is the sigmoid of the final logit:
$$
P(M | \Phi) = \frac{1}{1 + e^{-L_{final}}}
$$
