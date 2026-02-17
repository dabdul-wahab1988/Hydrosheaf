"""
Gaussian-Levenberg-Marquardt (GLM) optimization engine.
Wraps scipy.optimize.least_squares to mimic PEST++ functionality.
Supports parallel Jacobian calculation.
"""

from typing import Dict, List, Callable, Any, Optional, cast
import numpy as np
from scipy.optimize import least_squares
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

from .definitions import AdjustableParameter, Observation
from .problem import CalibrationProblem
from ..log import get_logger

logger = get_logger("calibration.glm")


class PESTGLM:
    """
    PEST++ GLM clone using scipy.optimize.
    """

    def __init__(
        self,
        parameters: List[AdjustableParameter],
        observations: List[Observation],
        model_runner: Callable[[Dict[str, float]], Dict[str, float]],
        n_workers: int = 1,  # Number of parallel workers
        worker_type: str = "thread",  # "thread" or "process"
    ):
        self.parameters = parameters
        self.observations = observations
        self.model_runner = model_runner
        self.n_workers = n_workers
        self.worker_type = worker_type

        # Internal state
        self.iteration = 0
        self.phi_history: List[float] = []
        self.param_history: List[Dict[str, float]] = []

        # Prepare arrays for scipy
        self.x0 = np.array([p.to_internal(p.value) for p in parameters])
        self.bounds = (
            [p.to_internal(p.lower_bound) for p in parameters],
            [p.to_internal(p.upper_bound) for p in parameters],
        )
        self.obs_values = np.array([o.value for o in observations])
        self.weights = np.array([o.weight for o in observations])
        self.obs_names = [o.name for o in observations]
        self.par_names = [p.name for p in parameters]

    @classmethod
    def from_problem(
        cls,
        problem: CalibrationProblem,
        n_workers: int = 1,
        worker_type: str = "thread",
    ) -> "PESTGLM":
        """
        Create PESTGLM instance from a CalibrationProblem.
        """
        return cls(
            parameters=problem.get_parameters(),
            observations=problem.get_observations(),
            model_runner=problem.run_model,
            n_workers=n_workers,
            worker_type=worker_type,
        )

    def _get_residuals(
        self, x: np.ndarray, model_result: Optional[Dict[str, float]] = None
    ) -> np.ndarray:
        """
        Calculate weighted residuals given internal parameters x.
        If model_result is provided, skip running the model.
        """
        # 1. Map internal x to model parameters
        param_dict = {}
        for i, p in enumerate(self.parameters):
            real_val = p.from_internal(x[i])
            param_dict[p.name] = real_val

        # 2. Run Model if needed
        if model_result is None:
            model_result = self.model_runner(param_dict)

        # 3. Measurement Residuals
        obs_residuals = []
        for i, o in enumerate(self.observations):
            sim_val = model_result.get(o.name, -9999.0)
            res = o.weight * (sim_val - o.value)
            obs_residuals.append(res)

        # 4. Tikhonov Regularization Residuals
        # PEST++ adds equations: w_reg * (P_current - P_prior) = 0
        reg_residuals = []
        for i, p in enumerate(self.parameters):
            if (
                p.prior_mean is not None
                and p.prior_sigma is not None
                and p.prior_sigma > 0
            ):
                # Weight = 1 / sigma
                w_reg = 1.0 / p.prior_sigma

                # Calculation in internal (transformed) space for stability
                current_val_internal = x[i]
                prior_val_internal = p.to_internal(p.prior_mean)

                res = w_reg * (current_val_internal - prior_val_internal)
                reg_residuals.append(res)

        # Combine
        return np.array(obs_residuals + reg_residuals)

    def _objective_function(self, x: np.ndarray) -> np.ndarray:
        """
        Function minimized by least_squares.
        """
        residuals = self._get_residuals(x)

        # Logging
        # residuals contains [obs_1, ..., obs_N, reg_1, ..., reg_M]
        n_obs = len(self.observations)
        obs_part = residuals[:n_obs]
        reg_part = residuals[n_obs:]

        phi_meas = np.sum(obs_part**2)
        phi_reg = np.sum(reg_part**2)
        phi_total = phi_meas + phi_reg

        # Track phi history for convergence plotting
        self.phi_history.append(phi_total)

        if self.iteration % 1 == 0:
            logger.info(
                f"Iter {self.iteration}: Phi_total = {phi_total:.4e} (Meas: {phi_meas:.4e}, Reg: {phi_reg:.4e})"
            )
        self.iteration += 1

        return residuals

    def _calculate_jacobian(
        self, x: np.ndarray, f0: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Parallel Jacobian Calculation using finite differences.
        J_ij = d(res_i) / d(x_j)
        """
        n_params = len(x)

        # 1. Base run (f0)
        # scipy sometimes passes f0, sometimes not. If we are providing 'jac', scipy expects
        # us to compute everything.
        # We need residuals at x.
        if f0 is None:
            f0 = self._get_residuals(x)

        n_residuals = len(f0)
        J = np.zeros((n_residuals, n_params))

        # 2. Perturbation steps
        # Use simple forward difference with step size relative to value
        # eps = 1e-6 (approx sqrt(machine_eps))
        eps = 1e-6
        steps = np.abs(x) * eps
        steps[steps < 1e-8] = 1e-8  # Min step

        # 3. Prepare Tasks
        tasks = []
        for j in range(n_params):
            x_perturbed = x.copy()
            x_perturbed[j] += steps[j]
            tasks.append(x_perturbed)

        # 4. Helper for parallel execution
        # We need to map tasks -> model_results
        # Then we compute residuals locally (fast)

        # Since self.model_runner might not be picklable, we have to be careful with ProcessPool.
        # But if we use ThreadPool, it's shared memory.

        def run_task(x_p):
            # Map x_p to param_dict
            p_dict = {}
            for i, p in enumerate(self.parameters):
                p_dict[p.name] = p.from_internal(x_p[i])
            return self.model_runner(p_dict)

        model_results = []

        if self.n_workers > 1:
            Executor = (
                ThreadPoolExecutor
                if self.worker_type == "thread"
                else ProcessPoolExecutor
            )
            try:
                with Executor(max_workers=self.n_workers) as executor:
                    model_results = list(executor.map(run_task, tasks))
            except Exception as e:
                logger.error(f"Parallel execution failed: {e}. Falling back to serial.")
                model_results = [run_task(t) for t in tasks]
        else:
            model_results = [run_task(t) for t in tasks]

        # 5. Assemble Jacobian
        for j in range(n_params):
            f_perturbed = self._get_residuals(tasks[j], model_result=model_results[j])

            # Forward difference
            # J[:, j] = (f(x+h) - f(x)) / h
            col = (f_perturbed - f0) / steps[j]
            J[:, j] = col

        return J

    def calibrate(self, max_nfev: int = 50) -> Dict[str, Any]:
        """
        Run the calibration.
        """
        logger.info(
            f"Starting Calibration with {len(self.parameters)} parameters and {len(self.observations)} observations..."
        )
        if self.n_workers > 1:
            logger.info(
                f"Parallel Mode Enabled: {self.n_workers} workers ({self.worker_type})"
            )

        # We use a lambda to pass 'self' to the jacobian method cleanly if needed,
        # or just pass the bound method.
        # NOTE: least_squares(jac=callable) requires callable(x).
        # But we computed f0 inside objective_function?
        # Ideally least_squares calls fun(x), then jac(x).
        # If we want to reuse f0 in jac(x), we'd need caching.
        # For simplicity, we let jac(x) re-compute f0 (1 extra run),
        # or we rely on the fact that for many steps, the evaluation of J dominates.
        # (N+1 runs vs 1 run).

        # To strictly use parallel jacobian, we must provide it to least_squares.

        jac_arg: Any = "2-point"  # Default serial
        if self.n_workers > 1:
            jac_arg = self._calculate_jacobian

        result = least_squares(
            self._objective_function,
            self.x0,
            jac=cast(Any, jac_arg),
            bounds=self.bounds,
            method="trf",
            loss="linear",
            max_nfev=max_nfev,
            verbose=1,
            ftol=1e-4,
            xtol=1e-4,
            gtol=1e-4,
        )

        # Post-processing
        optimal_params = {}
        for i, p in enumerate(self.parameters):
            optimal_params[p.name] = p.from_internal(result.x[i])

        # Jacobian at solution
        if self.n_workers > 1:
            # result.jac is just the last J returned by our function
            jac_weighted = result.jac
        else:
            jac_weighted = result.jac

        # Approximate Covariance Matrix
        n_obs = len(self.observations)
        n_par = len(self.parameters)
        dof = max(1, n_obs - n_par)
        phi_final = result.cost * 2.0
        sigma2 = phi_final / dof

        try:
            cov_x = np.linalg.inv(jac_weighted.T @ jac_weighted) * sigma2
            uncertainties = np.sqrt(np.diag(cov_x))
        except np.linalg.LinAlgError:
            cov_x = None
            uncertainties = np.zeros(n_par)

        # Map uncertainties back to real space
        param_uncertainties = {}
        for i, p in enumerate(self.parameters):
            unc_internal = uncertainties[i]
            val_internal = result.x[i]
            if p.log_transform:
                deriv = np.log(10) * (10**val_internal)
                real_unc = unc_internal * deriv
            else:
                real_unc = unc_internal
            param_uncertainties[p.name] = real_unc

        return {
            "optimal_parameters": optimal_params,
            "parameter_uncertainties_95pc": {
                k: v * 1.96 for k, v in param_uncertainties.items()
            },
            "phi": phi_final,
            "n_iterations": result.nfev,
            "success": result.success,
            "message": result.message,
        }


def sample_initial_parameters(
    parameters: List[AdjustableParameter],
    rng: np.random.Generator
) -> Dict[str, float]:
    params: Dict[str, float] = {}
    for p in parameters:
        low = p.lower_bound
        high = p.upper_bound
        if p.log_transform:
            low = max(low, 1e-12)
            high = max(high, low * 10.0)
            log_low = np.log10(low)
            log_high = np.log10(high)
            value = 10 ** rng.uniform(log_low, log_high)
        else:
            value = rng.uniform(low, high)
        params[p.name] = value
    return params


def run_pestglm_multistart(
    problem: CalibrationProblem,
    n_starts: int = 6,
    max_nfev: int = 80,
    seed: int = 42,
    n_workers: int = 1,
    worker_type: str = "thread",
    score_fn: Optional[Callable[[Dict[str, float]], float]] = None,
) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    base_params = {p.name: p.value for p in problem.get_parameters()}

    best_score: Optional[float] = None
    best_params: Optional[Dict[str, float]] = None
    best_result: Optional[Dict[str, Any]] = None

    start_results: List[Dict[str, Any]] = []

    for start_idx in range(n_starts):
        if start_idx == 0:
            init_params = base_params
        else:
            init_params = sample_initial_parameters(problem.get_parameters(), rng)

        for p in problem.get_parameters():
            p.value = init_params.get(p.name, p.value)

        pest = PESTGLM.from_problem(problem, n_workers=n_workers, worker_type=worker_type)
        try:
            result = pest.calibrate(max_nfev=max_nfev)
        except Exception as exc:
            start_results.append({
                "start": start_idx + 1,
                "success": False,
                "error": str(exc),
            })
            continue

        candidate = result.get("optimal_parameters", init_params)
        if score_fn is None:
            phi = float(result.get("phi", np.inf))
            score = -phi
        else:
            score = score_fn(candidate)

        if score is None or not np.isfinite(score):
            start_results.append({
                "start": start_idx + 1,
                "success": result.get("success", False),
                "score": score,
                "phi": result.get("phi"),
                "error": "non-finite score",
            })
            continue

        start_results.append({
            "start": start_idx + 1,
            "success": result.get("success", False),
            "score": score,
            "phi": result.get("phi"),
        })

        if best_score is None or score > best_score:
            best_score = score
            best_params = candidate
            best_result = result

    return {
        "best_parameters": best_params or base_params,
        "best_result": best_result,
        "best_score": best_score,
        "start_results": start_results,
    }
