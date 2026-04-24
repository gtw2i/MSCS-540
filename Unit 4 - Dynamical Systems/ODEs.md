# Ordinary Differential Equations: An Introduction

---

## 1. What Is a Differential Equation?

A **differential equation** is an equation that relates a function to its own derivatives. When the function depends on a single independent variable, we call it an **ordinary differential equation (ODE)**. When it depends on multiple variables (and involves partial derivatives), we get a **partial differential equation** — those are for another day.

The general form of a first-order ODE is:

$$\frac{dy}{dt} = f(t, y)$$

where $y = y(t)$ is the unknown function and $f$ is some expression involving $t$ and $y$.

The **order** of an ODE is the order of the highest derivative that appears. The **degree** is the power to which that highest derivative is raised (assuming it's polynomial). We'll focus almost entirely on first-order equations; higher-order systems can be reduced to first-order via a standard trick (more on that later).

---

## 2. Vocabulary and Classification

| Term | Meaning |
|---|---|
| **Autonomous** | $f$ does not depend explicitly on $t$: $\dot{y} = f(y)$ |
| **Linear** | $f$ is linear in $y$: $\dot{y} = a(t)y + b(t)$ |
| **Homogeneous (linear)** | $b(t) = 0$ |
| **Separable** | $f(t,y) = g(t) \cdot h(y)$ — can separate variables |
| **IVP** | Initial Value Problem — ODE + an initial condition $y(t_0) = y_0$ |

An **initial condition** pins down a specific solution from the infinite family of solutions that a general ODE admits. Without it, the general solution contains free constants.

---

## 3. Example 1 — Exponential Growth / Decay

Perhaps the simplest and most ubiquitous ODE:

$$\frac{dy}{dt} = r \, y, \qquad y(0) = y_0$$

This says: *the rate of change of $y$ is proportional to $y$ itself.*

**Separation of variables:**

$$\frac{dy}{y} = r \, dt \implies \ln |y| = rt + C \implies y(t) = y_0 e^{rt}$$

- If $r > 0$: exponential growth (population, compound interest, early-stage epidemic)
- If $r < 0$: exponential decay (radioactive decay, cooling)

The **doubling time** when $r > 0$ is $T_2 = \frac{\ln 2}{r}$, and the **half-life** when $r < 0$ is $T_{1/2} = \frac{\ln 2}{|r|}$.

---

## 4. Example 2 — Logistic Growth

Exponential growth is unrealistic long-term — resources are finite. The **logistic equation** introduces a carrying capacity $K$:

$$\frac{dy}{dt} = r \, y \left(1 - \frac{y}{K}\right), \qquad y(0) = y_0$$

**Qualitative behavior:**
- When $y \ll K$: the factor $(1 - y/K) \approx 1$, so growth is approximately exponential.
- When $y \to K$: the factor vanishes and $\dot{y} \to 0$ — the population stabilizes.
- **Fixed points (equilibria):** $y^* = 0$ (unstable) and $y^* = K$ (stable).

**Exact solution** (via partial fractions):

$$y(t) = \frac{K \, y_0}{y_0 + (K - y_0) e^{-rt}}$$

The curve is the familiar **S-shaped (sigmoidal)** curve. The inflection point — where growth is fastest — occurs at $y = K/2$.

---

## 5. Fixed Points and Stability

For an autonomous ODE $\dot{y} = f(y)$, a **fixed point** (or equilibrium) is a value $y^*$ such that $f(y^*) = 0$.

**Stability** is determined by the sign of $f'(y^*)$:

$$\text{Stable (attracting): } f'(y^*) < 0 \qquad \text{Unstable (repelling): } f'(y^*) > 0$$

Intuitively: if $f$ is decreasing through zero, nearby trajectories are pushed back toward $y^*$; if increasing, they flee.

**Phase line analysis:** Draw the real line for $y$, mark fixed points, and indicate the sign of $f(y)$ (hence the direction of $\dot{y}$) in each interval. This gives a complete qualitative picture without solving the ODE.

---

## 6. Systems of ODEs and Higher-Order Equations

A second-order ODE like Newton's law $m\ddot{x} = F(t, x, \dot{x})$ can be rewritten as a **first-order system** by introducing a new variable $v = \dot{x}$:

$$\dot{x} = v, \qquad \dot{v} = \frac{F(t, x, v)}{m}$$

More generally, any $n$-th order ODE can be reduced to a system of $n$ first-order equations. In vector form:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}), \qquad \mathbf{y} \in \mathbb{R}^n$$

This is the standard form that numerical solvers operate on. **Everything reduces to this.**

---

## 7. Numerical Methods

Most ODEs that arise in practice — especially those coupled to network structure or real data — do not have closed-form solutions. We approximate them numerically.

The setup: given $\dot{y} = f(t, y)$ and $y(t_0) = y_0$, we want to compute $y$ at times $t_0 < t_1 < t_2 < \cdots$, where the step size is $h = t_{k+1} - t_k$.

---

### 7.1 Euler's Method

The simplest possible scheme: use the current slope to extrapolate forward.

**Derivation** from Taylor expansion:

$$y(t + h) = y(t) + h \, \dot{y}(t) + \mathcal{O}(h^2) = y(t) + h \, f(t, y(t)) + \mathcal{O}(h^2)$$

Dropping the error term gives the **update rule**:

$$\boxed{y_{k+1} = y_k + h \cdot f(t_k, y_k)}$$

**Local truncation error:** $\mathcal{O}(h^2)$ per step.  
**Global error:** $\mathcal{O}(h)$ — Euler's method is first-order accurate.

**Pros:** Simple, cheap, easy to implement.  
**Cons:** Accumulates error quickly; can be unstable for stiff problems or large $h$.

---

### 7.2 Midpoint Method (Modified Euler / RK2)

Instead of using the slope at the start of the interval, take a half-step first to estimate the slope at the midpoint:

$$k_1 = f(t_k, y_k)$$
$$k_2 = f\!\left(t_k + \tfrac{h}{2},\; y_k + \tfrac{h}{2} k_1\right)$$
$$\boxed{y_{k+1} = y_k + h \cdot k_2}$$

**Global error:** $\mathcal{O}(h^2)$ — second-order accurate. Substantially better than Euler for the same step size.

---

### 7.3 Runge-Kutta 4 (RK4)

The workhorse of numerical ODE solving. Uses four slope estimates per step:

$$k_1 = f(t_k,\; y_k)$$
$$k_2 = f\!\left(t_k + \tfrac{h}{2},\; y_k + \tfrac{h}{2} k_1\right)$$
$$k_3 = f\!\left(t_k + \tfrac{h}{2},\; y_k + \tfrac{h}{2} k_2\right)$$
$$k_4 = f(t_k + h,\; y_k + h \, k_3)$$

$$\boxed{y_{k+1} = y_k + \frac{h}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right)}$$

The update is a **weighted average** of the four slopes: the midpoint estimates ($k_2$, $k_3$) receive twice the weight of the endpoint ones.

**Global error:** $\mathcal{O}(h^4)$ — fourth-order accurate. Halving $h$ reduces error by a factor of $\sim 16$.

---

### 7.4 Comparison Table

| Method | Stages (evals of $f$) | Global Error | Typical Use |
|---|---|---|---|
| Euler | 1 | $\mathcal{O}(h)$ | Teaching; quick prototypes |
| Midpoint (RK2) | 2 | $\mathcal{O}(h^2)$ | Better accuracy, low overhead |
| RK4 | 4 | $\mathcal{O}(h^4)$ | Standard for most problems |
| `scipy.solve_ivp` (RK45) | adaptive | $\mathcal{O}(h^5)$ | Production; handles stiffness |

---

### 7.5 Stiffness and Adaptive Step Size

A system is **stiff** when it contains dynamics operating on wildly different timescales. Euler and fixed-step RK4 may require extremely small $h$ to remain stable, making them impractical. Modern solvers like `scipy.solve_ivp` with `method='Radau'` or `'BDF'` handle stiff systems efficiently.

**Adaptive step size control:** the solver estimates the local error at each step and shrinks or enlarges $h$ to keep it within a tolerance. You specify `rtol` (relative tolerance) and `atol` (absolute tolerance); the solver does the rest.

```python
from scipy.integrate import solve_ivp

def f(t, y):
    return -0.5 * y  # exponential decay

sol = solve_ivp(f, t_span=(0, 10), y0=[1.0], dense_output=True)
```

---

## 8. Existence and Uniqueness (Brief)

The **Picard-Lindelöf theorem** guarantees that if $f(t, y)$ is continuous and **Lipschitz continuous** in $y$ near $(t_0, y_0)$, then the IVP $\dot{y} = f(t,y)$, $y(t_0) = y_0$ has a **unique solution** in some interval around $t_0$.

Lipschitz in $y$ means: there exists $L$ such that $|f(t, y_1) - f(t, y_2)| \leq L |y_1 - y_2|$ for all $y_1, y_2$ in a neighborhood. If $\partial f / \partial y$ is bounded, this is satisfied.

**Why does this matter practically?** It tells you when you can trust a numerical solution to be tracking the unique true solution, rather than an artifact of the algorithm.

---

## 9. ODEs on Networks (Preview)

The real payoff for this course: when we put ODEs *on a graph*, the graph structure becomes part of the dynamics. The prototypical example is **diffusion on a graph** (Laplacian flow):

$$\frac{d\mathbf{y}}{dt} = -L \, \mathbf{y}$$

where $\mathbf{y} \in \mathbb{R}^n$ is a value at each node and $L$ is the **graph Laplacian**. The solution is:

$$\mathbf{y}(t) = e^{-Lt} \mathbf{y}(0)$$

The spectral decomposition of $L$ — its eigenvalues and eigenvectors — completely governs how information spreads across the network and at what rates. This is the bridge between the ODE module and the core of the course.

---

## Summary

| Concept | Key Idea |
|---|---|
| ODE | Relates $y(t)$ to its derivatives |
| IVP | ODE + initial condition $\Rightarrow$ unique solution |
| Separable | Split $f(t,y) = g(t)h(y)$; integrate both sides |
| Fixed points | $f(y^*) = 0$; stable if $f'(y^*) < 0$ |
| Higher-order $\to$ system | Introduce velocity variables; reduce to $\dot{\mathbf{y}} = \mathbf{f}$ |
| Euler | $y_{k+1} = y_k + h f(t_k, y_k)$; first-order |
| RK4 | Four-stage weighted average; fourth-order |
| `solve_ivp` | Adaptive, production-grade solver |
| ODEs on graphs | $\dot{\mathbf{y}} = -L\mathbf{y}$; connects to Laplacian and spectral theory |