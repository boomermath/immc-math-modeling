import csv
import tempfile
import unittest
from pathlib import Path

from flytrap_ipm import burn_history, integrate_burn_events, read_rows


class BurnIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.rows = read_rows(Path(__file__).parent / "Datasets" / "demoVFT.csv")

    def event_file(self, data):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "events.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["site", "burn_date"])
            writer.writeheader()
            writer.writerows(data)
        return path

    def test_existing_burn_ages_reset_after_fire(self):
        history = burn_history(self.rows)
        self.assertEqual([history["coastal.1", year] for year in (2016, 2017, 2018)],
                         [1, 2, 1])

    def test_matching_event_is_accepted(self):
        path = self.event_file([{"site": "coastal.1", "burn_date": "2018-04-10"}])
        self.assertEqual(integrate_burn_events(self.rows, path)["coastal.1", 2018], 1)

    def test_conflicting_event_is_rejected(self):
        path = self.event_file([{"site": "coastal.1", "burn_date": "2017-04-10"}])
        with self.assertRaisesRegex(ValueError, "Burn-history conflict"):
            integrate_burn_events(self.rows, path)

    def test_other_study_sites_are_rejected(self):
        path = self.event_file([{"site": "S1", "burn_date": "2017-03-29"}])
        with self.assertRaisesRegex(ValueError, "Unknown burn-event site"):
            integrate_burn_events(self.rows, path)


if __name__ == "__main__":
    unittest.main()
