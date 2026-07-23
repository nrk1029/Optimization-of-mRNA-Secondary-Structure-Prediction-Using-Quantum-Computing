import unittest

from src.evaluation.structure_metrics import compare_structures
from src.rna_utils import normalize_sequence, pairs_to_dot_bracket, validate_structure


class RNAUtilitiesTest(unittest.TestCase):
    def test_normalizes_dna_and_whitespace(self):
        self.assertEqual(normalize_sequence(" aug tc\n"), "AUGUC")

    def test_rejects_invalid_nucleotide(self):
        with self.assertRaises(ValueError):
            normalize_sequence("AUGX")

    def test_validates_and_round_trips_structure(self):
        sequence = "GGGAAACCC"
        structure = "(((...)))"
        pairs = validate_structure(sequence, structure)
        self.assertEqual(pairs_to_dot_bracket(len(sequence), pairs), structure)

    def test_metrics(self):
        metrics = compare_structures("((....))", "(......)")
        self.assertEqual(metrics.true_positives, 1)
        self.assertEqual(metrics.false_negatives, 1)
        self.assertAlmostEqual(metrics.precision, 1.0)
        self.assertAlmostEqual(metrics.recall, 0.5)
        self.assertAlmostEqual(metrics.f1, 2 / 3)

    def test_crossing_pairs_are_rejected(self):
        with self.assertRaises(ValueError):
            pairs_to_dot_bracket(6, {(0, 3), (1, 5)})


if __name__ == "__main__":
    unittest.main()
