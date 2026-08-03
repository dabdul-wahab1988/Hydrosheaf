import unittest
from datetime import datetime
from unittest import mock

from hydrosheaf.api import (
    PipelineStageError,
    auto_disable_missing_modules,
    fit_network_pipeline,
    fit_network_with_priors,
)
from hydrosheaf.config import Config


class AutoDisableMissingModulesTests(unittest.TestCase):
    def test_disables_all_when_data_missing(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A"}]

        updated = auto_disable_missing_modules(samples, config)

        self.assertFalse(updated.phreeqc_enabled)
        self.assertFalse(updated.isotope_enabled)
        self.assertFalse(updated.nitrate_source_enabled)

    def test_keeps_isotopes_when_present(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A", "18O": -5.0, "2H": -35.0}]

        updated = auto_disable_missing_modules(samples, config)

        self.assertFalse(updated.phreeqc_enabled)
        self.assertTrue(updated.isotope_enabled)
        self.assertFalse(updated.nitrate_source_enabled)

    def test_keeps_phreeqc_and_nitrate_when_present(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A", "pH": 7.2, "NO3": 1.5}]

        updated = auto_disable_missing_modules(samples, config)

        self.assertTrue(updated.phreeqc_enabled)
        self.assertFalse(updated.isotope_enabled)
        self.assertTrue(updated.nitrate_source_enabled)

    def test_fit_network_with_priors_auto_disable(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A"}, {"site_id": "B"}]
        edges = [("A", "B")]
        captured = {}

        def _fake_fit_network(samples_arg, edges_arg, config_arg, **kwargs):
            captured["config"] = config_arg
            return []

        with mock.patch("hydrosheaf.api.fit_network", side_effect=_fake_fit_network):
            fit_network_with_priors(samples, edges, config)

        self.assertIn("config", captured)
        self.assertFalse(captured["config"].phreeqc_enabled)
        self.assertFalse(captured["config"].isotope_enabled)
        self.assertFalse(captured["config"].nitrate_source_enabled)

    def test_fit_network_pipeline_auto_disable(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A"}, {"site_id": "B"}]
        edges = [("A", "B")]
        captured = {}

        def _fake_fit_network(samples_arg, edges_arg, config_arg, **kwargs):
            captured["config"] = config_arg
            return []

        with mock.patch("hydrosheaf.api.fit_network", side_effect=_fake_fit_network):
            results, extras = fit_network_pipeline(samples, edges, config)

        self.assertEqual(results, [])
        self.assertIn("config", captured)
        self.assertFalse(captured["config"].phreeqc_enabled)
        self.assertFalse(captured["config"].isotope_enabled)
        self.assertFalse(captured["config"].nitrate_source_enabled)
        self.assertIn("edges", extras)
        self.assertEqual(
            extras["stage_status"]["network_fit"]["status"], "completed"
        )

    def test_pipeline_does_not_mutate_caller_samples_for_latent_nodes(self):
        config = Config(
            phreeqc_enabled=False,
            isotope_enabled=False,
            nitrate_source_enabled=False,
            latent_endmembers_enabled=True,
        )
        samples = [{"site_id": "A"}, {"site_id": "B"}]
        original = [dict(row) for row in samples]

        with mock.patch(
            "hydrosheaf.api.identify_latent_endmembers",
            return_value=[{"site_id": "latent_1"}],
        ):
            with mock.patch("hydrosheaf.api.fit_network", return_value=[]):
                _, extras = fit_network_pipeline(samples, [("A", "B")], config)

        self.assertEqual(samples, original)
        self.assertEqual(len(extras["virtual_nodes"]), 1)
        self.assertEqual(
            extras["stage_status"]["latent_endmembers"]["status"],
            "completed",
        )

    def test_strict_pipeline_fails_when_requested_nuclear_stage_is_skipped(self):
        config = Config(
            phreeqc_enabled=False,
            isotope_enabled=False,
            nitrate_source_enabled=False,
            sheaf_age_enabled=True,
            residence_time_tracer="3H",
        )
        samples = [{"site_id": "A"}, {"site_id": "B"}]

        with mock.patch("hydrosheaf.api.fit_network") as mocked_fit:
            with self.assertRaisesRegex(PipelineStageError, "nuclear_age"):
                fit_network_pipeline(
                    samples,
                    [("A", "B")],
                    config,
                    strict_stage_completion=True,
                )

        mocked_fit.assert_not_called()

    def test_strict_pipeline_fails_when_nuclear_solver_returns_no_posteriors(self):
        config = Config(
            phreeqc_enabled=False,
            isotope_enabled=False,
            nitrate_source_enabled=False,
            sheaf_age_enabled=True,
            residence_time_tracer="3H",
        )
        samples = [
            {"site_id": "A", "3H": 5.0, "sample_date": "2020-01-01"},
            {"site_id": "B", "3H": 2.0, "sample_date": "2020-01-01"},
        ]

        with mock.patch(
            "hydrosheaf.api.infer_network_ages_bayesian", return_value={}
        ):
            with mock.patch("hydrosheaf.api.fit_network") as mocked_fit:
                with self.assertRaisesRegex(PipelineStageError, "nuclear_age"):
                    fit_network_pipeline(
                        samples,
                        [("A", "B")],
                        config,
                        strict_stage_completion=True,
                    )

        mocked_fit.assert_not_called()

    def test_strict_pipeline_rejects_required_stage_not_requested(self):
        config = Config(
            phreeqc_enabled=False,
            isotope_enabled=False,
            nitrate_source_enabled=False,
            sheaf_age_enabled=False,
        )

        with self.assertRaisesRegex(ValueError, "not requested"):
            fit_network_pipeline(
                [{"site_id": "A"}, {"site_id": "B"}],
                [("A", "B")],
                config,
                strict_stage_completion=True,
                required_stages=["sheaf_refinement"],
            )

    def test_strict_pipeline_runs_nuclear_before_sheaf_refinement(self):
        config = Config(
            phreeqc_enabled=False,
            isotope_enabled=False,
            nitrate_source_enabled=False,
            sheaf_age_enabled=True,
            residence_time_tracer="3H",
        )
        samples = [
            {"site_id": "A", "3H": 5.0, "sample_date": "2020-01-01"},
            {"site_id": "B", "3H": 2.0, "sample_date": "2020-01-01"},
        ]
        posterior = {
            "A": {
                "mean_age_years": 3.0,
                "std_age_years": 0.5,
                "tracer_identifiable": True,
            },
            "B": {
                "mean_age_years": 7.0,
                "std_age_years": 0.8,
                "tracer_identifiable": True,
            },
        }

        def _refine(rows, edges, _config):
            by_id = {row["site_id"]: row for row in rows}
            self.assertEqual(by_id["A"]["mean_age_years"], 3.0)
            self.assertEqual(by_id["B"]["mean_age_std_years"], 0.8)
            return list(edges)

        with mock.patch(
            "hydrosheaf.api.infer_network_ages_bayesian", return_value=posterior
        ) as mocked_infer:
            with mock.patch(
                "hydrosheaf.api.refine_edges_with_sheaf", side_effect=_refine
            ):
                with mock.patch("hydrosheaf.api.fit_network", return_value=[]):
                    _, extras = fit_network_pipeline(
                        samples,
                        [("A", "B")],
                        config,
                        sheaf_refinement_enabled=True,
                        nuclear_inference_options={"sampler": "grid", "n_samples": 50},
                        strict_stage_completion=True,
                    )

        self.assertEqual(samples[0].get("mean_age_years"), None)
        infer_args, infer_kwargs = mocked_infer.call_args
        self.assertEqual(infer_args[0].number_of_edges(), 0)
        self.assertEqual(infer_kwargs["sampler"], "grid")
        self.assertEqual(infer_kwargs["n_samples"], 50)
        for stage in ("nuclear_age", "sheaf_refinement", "network_fit"):
            self.assertEqual(extras["stage_status"][stage]["status"], "completed")

    def test_fit_network_pipeline_nuclear_uses_median_sample_date(self):
        config = Config(
            phreeqc_enabled=False,
            isotope_enabled=False,
            nitrate_source_enabled=False,
            sheaf_age_enabled=True,
            residence_time_tracer="3H",
        )
        samples = [
            {"site_id": "A", "3H": 5.0, "sample_date": "2020-01-01"},
            {"site_id": "B", "3H": 2.0, "sample_date": "2022-01-01"},
        ]
        edges = [("A", "B")]

        with mock.patch("hydrosheaf.api.fit_network", return_value=[]):
            with mock.patch(
                "hydrosheaf.api.infer_network_ages_bayesian", return_value={}
            ) as mocked_infer:
                fit_network_pipeline(samples, edges, config)

        self.assertTrue(mocked_infer.called)
        args, _ = mocked_infer.call_args
        self.assertAlmostEqual(float(args[3]), 2021.0, places=6)

    def test_fit_network_pipeline_nuclear_falls_back_to_current_year(self):
        config = Config(
            phreeqc_enabled=False,
            isotope_enabled=False,
            nitrate_source_enabled=False,
            sheaf_age_enabled=True,
            residence_time_tracer="3H",
        )
        samples = [
            {"site_id": "A", "3H": 5.0},
            {"site_id": "B", "3H": 2.0},
        ]
        edges = [("A", "B")]

        with mock.patch("hydrosheaf.api.fit_network", return_value=[]):
            with mock.patch(
                "hydrosheaf.api.infer_network_ages_bayesian", return_value={}
            ) as mocked_infer:
                fit_network_pipeline(samples, edges, config)

        self.assertTrue(mocked_infer.called)
        args, _ = mocked_infer.call_args
        self.assertEqual(float(args[3]), float(datetime.now().year))

    def test_fit_network_with_priors_strict_requires_phreeqc_major_ions(self):
        config = Config(phreeqc_enabled=True, strict_input_validation=True)
        samples = [{"site_id": "A", "pH": 7.1}, {"site_id": "B", "pH": 7.3}]
        edges = [("A", "B")]

        with mock.patch("hydrosheaf.api.fit_network") as mocked_fit:
            with self.assertRaises(ValueError):
                fit_network_with_priors(samples, edges, config)

        mocked_fit.assert_not_called()


if __name__ == "__main__":
    unittest.main()
