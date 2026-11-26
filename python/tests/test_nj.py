import math

import pytest
from nj_py import average_distance, distance_matrix, nj


def msa(*pairs):
    return [{"identifier": name, "sequence": seq} for name, seq in pairs]


class TestBasic:
    def test_returns_dict_with_newick(self):
        result = nj(msa(("A", "ACGT"), ("B", "ACGA")))
        assert isinstance(result, dict)
        assert "newick" in result
        assert isinstance(result["newick"], str)

    def test_newick_format(self):
        result = nj(msa(("A", "ACGT"), ("B", "ACGA")))
        assert result["newick"].startswith("(")
        assert result["newick"].endswith(";")

    def test_two_identical_taxa_zero_distance(self):
        result = nj(msa(("A", "ACGT"), ("B", "ACGT")))
        assert result["newick"] == "(A:0.000,B:0.000);"

    def test_three_taxa_contains_all_leaves(self):
        result = nj(msa(("Sp1", "ACGT"), ("Sp2", "ACGA"), ("Sp3", "AGGT")))
        assert "Sp1" in result["newick"]
        assert "Sp2" in result["newick"]
        assert "Sp3" in result["newick"]

    def test_deterministic(self):
        sequences = msa(("A", "ACGTCG"), ("B", "ACG-GC"), ("C", "ACGCGT"))
        assert nj(sequences)["newick"] == nj(sequences)["newick"]


class TestSubstitutionModels:
    def test_jukes_cantor_dna(self):
        result = nj(
            msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")),
            substitution_model="JukesCantor",
        )
        assert result["newick"].endswith(";")

    def test_kimura2p_dna(self):
        result = nj(
            msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")),
            substitution_model="Kimura2P",
        )
        assert result["newick"].endswith(";")

    def test_pdiff_protein(self):
        result = nj(
            msa(("A", "ACDEFGH"), ("B", "ACDEFGK"), ("C", "ACDLFGH")),
            substitution_model="PDiff",
        )
        assert result["newick"].endswith(";")

    def test_poisson_protein(self):
        result = nj(
            msa(("A", "ACDEFGH"), ("B", "ACDEFGK"), ("C", "ACDLFGH")),
            substitution_model="Poisson",
        )
        assert result["newick"].endswith(";")


class TestErrors:
    def test_empty_msa_raises(self):
        with pytest.raises(ValueError, match="[Ee]mpty"):
            nj([])

    def test_poisson_on_dna_raises(self):
        with pytest.raises(ValueError):
            nj(msa(("A", "ACGT"), ("B", "ACGA")), substitution_model="Poisson")

    def test_jukes_cantor_on_protein_raises(self):
        with pytest.raises(ValueError):
            nj(msa(("A", "ACDEFGH"), ("B", "ACDEFGK")), substitution_model="JukesCantor")

    def test_kimura2p_on_protein_raises(self):
        with pytest.raises(ValueError):
            nj(msa(("A", "ACDEFGH"), ("B", "ACDEFGK")), substitution_model="Kimura2P")


class TestMixedCase:
    def test_lowercase_dna_detected_as_dna(self):
        result = nj(msa(("A", "acgt"), ("B", "acga")))
        assert result["newick"].endswith(";")

    def test_mixedcase_dna(self):
        result = nj(msa(("A", "AcGt"), ("B", "aCgA")))
        assert result["newick"].endswith(";")

    def test_lowercase_matches_uppercase(self):
        upper = nj(msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")))
        lower = nj(msa(("A", "acgt"), ("B", "acga"), ("C", "aggt")))
        assert upper["newick"] == lower["newick"]


class TestProteinCase:
    def test_lowercase_protein_matches_uppercase(self):
        upper = nj(
            msa(("A", "ACDEFGH"), ("B", "ACDEFGK"), ("C", "ACDQFGH")),
            substitution_model="Poisson",
        )
        lower = nj(
            msa(("A", "acdefgh"), ("B", "acdefgk"), ("C", "acdqfgh")),
            substitution_model="Poisson",
        )
        assert upper["newick"] == lower["newick"]

    def test_lowercase_protein_non_degenerate(self):
        # Distinct lowercase protein sequences must not collapse to distance 0.
        result = distance_matrix(
            msa(("A", "acdefgh"), ("B", "acdefgk")), substitution_model="Poisson"
        )
        assert result["matrix"][0][1] > 0.0


class TestSpecialNames:
    def test_sequence_name_with_spaces(self):
        result = nj(msa(("Seq A", "ACGT"), ("Seq B", "ACGA")))
        # Names with whitespace are single-quoted per the Newick grammar.
        assert "'Seq A'" in result["newick"]
        assert "'Seq B'" in result["newick"]

    def test_sequence_name_with_special_chars(self):
        # The pipe character is not reserved in Newick, so it stays unquoted.
        result = nj(msa(("Sp|1.1", "ACGT"), ("Sp|2.1", "ACGA")))
        assert "Sp|1.1" in result["newick"]
        assert "Sp|2.1" in result["newick"]

    def test_duplicate_identifiers_raise(self):
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            nj(msa(("A", "ACGT"), ("A", "ACGA")))


class TestDistanceMatrix:
    def test_returns_names_and_matrix_keys(self):
        result = distance_matrix(msa(("A", "ACGT"), ("B", "ACGA")))
        assert "names" in result
        assert "matrix" in result

    def test_names_match_input(self):
        result = distance_matrix(msa(("human", "ACGT"), ("chimp", "ACGA")))
        assert result["names"] == ["human", "chimp"]

    def test_matrix_is_n_by_n(self):
        result = distance_matrix(msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")))
        assert len(result["matrix"]) == 3
        assert all(len(row) == 3 for row in result["matrix"])

    def test_diagonal_is_zero(self):
        result = distance_matrix(msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")))
        for i in range(3):
            assert result["matrix"][i][i] == 0.0

    def test_matrix_is_symmetric(self):
        result = distance_matrix(msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")))
        for i in range(3):
            for j in range(3):
                assert result["matrix"][i][j] == result["matrix"][j][i]

    def test_pdiff_known_value(self):
        # one difference (T vs A at position 3) out of 4 → 0.25
        result = distance_matrix(msa(("A", "ACGT"), ("B", "ACGA")))
        assert abs(result["matrix"][0][1] - 0.25) < 1e-12

    def test_identical_sequences_zero_distance(self):
        result = distance_matrix(msa(("A", "ACGT"), ("B", "ACGT")))
        assert result["matrix"][0][1] == 0.0

    def test_jukes_cantor_known_value(self):
        # p=0.25 → d = -0.75 * ln(1 - 4/3 * 0.25)
        result = distance_matrix(msa(("A", "ACGT"), ("B", "ACGA")), substitution_model="JukesCantor")
        expected = -0.75 * math.log(1.0 - (4.0 / 3.0) * 0.25)
        assert abs(result["matrix"][0][1] - expected) < 1e-10

    def test_kimura2p_dna(self):
        result = distance_matrix(
            msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")),
            substitution_model="Kimura2P",
        )
        assert result["names"] == ["A", "B", "C"]
        assert result["matrix"][0][0] == 0.0

    def test_poisson_protein_known_value(self):
        # 1 diff (H vs K) out of 7 → p=1/7, d=-ln(1-1/7)
        result = distance_matrix(msa(("A", "ACDEFGH"), ("B", "ACDEFGK")), substitution_model="Poisson")
        expected = -math.log(1.0 - 1.0 / 7.0)
        assert abs(result["matrix"][0][1] - expected) < 1e-10

    def test_pdiff_protein(self):
        result = distance_matrix(msa(("A", "ACDEFGH"), ("B", "ACDEFGK")), substitution_model="PDiff")
        assert abs(result["matrix"][0][1] - 1.0 / 7.0) < 1e-12

    def test_empty_msa_raises(self):
        with pytest.raises(ValueError):
            distance_matrix([])

    def test_incompatible_model_raises(self):
        with pytest.raises(ValueError):
            distance_matrix(msa(("A", "ACGT"), ("B", "ACGA")), substitution_model="Poisson")

    def test_jukes_cantor_on_protein_raises(self):
        with pytest.raises(ValueError):
            distance_matrix(msa(("A", "ACDEFGH"), ("B", "ACDEFGK")), substitution_model="JukesCantor")


class TestAverageDistance:
    def test_returns_float(self):
        result = average_distance(msa(("A", "ACGT"), ("B", "ACGA")))
        assert isinstance(result, float)

    def test_identical_sequences_zero(self):
        result = average_distance(msa(("A", "ACGT"), ("B", "ACGT")))
        assert result == 0.0

    def test_two_taxa_equals_pairwise(self):
        # one difference out of 4 → 0.25
        result = average_distance(msa(("A", "ACGT"), ("B", "ACGA")))
        assert abs(result - 0.25) < 1e-12

    def test_three_taxa_known_value(self):
        # A↔B: 1/4=0.25 (T→A), A↔C: 1/4=0.25 (C→G), B↔C: 2/4=0.5 → avg = 1/3
        result = average_distance(msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")))
        assert abs(result - 1.0 / 3.0) < 1e-12

    def test_jukes_cantor_dna(self):
        result = average_distance(msa(("A", "ACGT"), ("B", "ACGA")), substitution_model="JukesCantor")
        expected = -0.75 * math.log(1.0 - (4.0 / 3.0) * 0.25)
        assert abs(result - expected) < 1e-10

    def test_empty_msa_raises(self):
        with pytest.raises(ValueError):
            average_distance([])

    def test_incompatible_model_raises(self):
        with pytest.raises(ValueError):
            average_distance(msa(("A", "ACGT"), ("B", "ACGA")), substitution_model="Poisson")


class TestBootstrap:
    def test_bootstrap_zero_samples(self):
        result = nj(
            msa(("A", "ACGTACGT"), ("B", "ACGAACGA"), ("C", "AGGTCGGT"), ("D", "TGGTCGGT")),
            n_bootstrap_samples=0,
        )
        assert result["newick"].endswith(";")

    def test_bootstrap_nonzero_returns_valid_newick(self):
        result = nj(
            msa(("A", "ACGTACGT"), ("B", "ACGAACGA"), ("C", "AGGTCGGT"), ("D", "TGGTCGGT")),
            n_bootstrap_samples=10,
        )
        assert isinstance(result["newick"], str)
        assert result["newick"].endswith(";")

    def test_on_event_called_for_bootstrap_progress(self):
        progress_events = []
        nj(
            msa(("A", "ACGTACGT"), ("B", "ACGAACGA"), ("C", "AGGTCGGT"), ("D", "TGGTCGGT")),
            n_bootstrap_samples=5,
            on_event=lambda event: progress_events.append(event)
            if event["type"] == "BootstrapProgress"
            else None,
        )
        assert len(progress_events) == 5
        assert progress_events[-1] == {"type": "BootstrapProgress", "completed": 5, "total": 5}


class TestReturnOptions:
    def test_return_distance_matrix(self):
        result = nj(
            msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")),
            return_distance_matrix=True,
        )
        assert "newick" in result
        assert "distance_matrix" in result
        assert result["distance_matrix"]["names"] == ["A", "B", "C"]

    def test_return_average_distance(self):
        result = nj(
            msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")),
            return_average_distance=True,
        )
        assert "newick" in result
        assert "average_distance" in result
        assert isinstance(result["average_distance"], float)

    def test_no_extra_keys_by_default(self):
        result = nj(msa(("A", "ACGT"), ("B", "ACGA")))
        assert "distance_matrix" not in result
        assert "average_distance" not in result


class TestRateVariation:
    SEQS = (
        ("A", "ACGTACGTACGTACGTACGTACGT"),
        ("B", "ACGTACGTACGAACGTAGGTACGT"),
        ("C", "ACGTTCGTACGTACATACGTACCT"),
        ("D", "ACGTACGTTCGTACGTACGTACGA"),
    )

    def test_gamma_changes_tree(self):
        plain = nj(msa(*self.SEQS), substitution_model="JukesCantor")
        gamma = nj(msa(*self.SEQS), substitution_model="JukesCantor", gamma_shape=0.5)
        assert plain["newick"] != gamma["newick"]

    def test_gamma_inflates_average_distance(self):
        plain = average_distance(msa(*self.SEQS), substitution_model="JukesCantor")
        gamma = average_distance(
            msa(*self.SEQS), substitution_model="JukesCantor", gamma_shape=0.5
        )
        assert gamma > plain

    def test_p_invar_accepted(self):
        result = nj(msa(*self.SEQS), substitution_model="JukesCantor", p_invar=0.3)
        assert result["newick"].endswith(";")

    def test_distance_matrix_accepts_rate_params(self):
        result = distance_matrix(
            msa(*self.SEQS), substitution_model="JukesCantor", gamma_shape=0.5
        )
        assert len(result["names"]) == 4

    def test_invalid_gamma_shape_raises(self):
        with pytest.raises(ValueError, match="[Gg]amma"):
            nj(msa(*self.SEQS), substitution_model="JukesCantor", gamma_shape=0.0)

    def test_invalid_p_invar_raises(self):
        with pytest.raises(ValueError, match="invariant"):
            nj(msa(*self.SEQS), substitution_model="JukesCantor", p_invar=1.0)


class TestSaturatedPairs:
    """Two groups of near-identical sequences that differ at every site between
    groups: under Jukes-Cantor the between-group pairs have no finite distance."""

    SEQS = (
        ("A", "ACGTACGTACGTACGTACGT"),
        ("B", "CATGCATGCATGCATGCATG"),
        ("C", "ACGTACGTACGAACGTACGT"),
        ("D", "ACGTACGTTCGTACGTACGT"),
        ("E", "CATGCATGCATCCATGCATG"),
    )

    def run(self, **kwargs):
        events = []
        result = nj(
            msa(*self.SEQS),
            substitution_model="JukesCantor",
            return_distance_matrix=True,
            return_average_distance=True,
            on_event=events.append,
            **kwargs,
        )
        warnings = [
            e["message"]
            for e in events
            if e["type"] == "Log" and e["level"] == "Warning"
        ]
        return result, warnings

    def test_newick_is_finite_and_warns(self):
        result, warnings = self.run()
        assert "NaN" not in result["newick"] and "inf" not in result["newick"]
        assert len(warnings) == 1 and "6 sequence pair(s)" in warnings[0]

    def test_default_cap_is_ten(self):
        result, _ = self.run()
        matrix = result["distance_matrix"]["matrix"]
        assert matrix[0][1] == 10.0
        assert all(math.isfinite(d) for row in matrix for d in row)
        assert math.isfinite(result["average_distance"])

    def test_custom_cap(self):
        result, warnings = self.run(max_distance=3.0)
        assert result["distance_matrix"]["matrix"][0][1] == 3.0
        assert "max_distance = 3" in warnings[0]

    def test_bootstrap_completes(self):
        result, warnings = self.run(n_bootstrap_samples=5)
        assert result["newick"].endswith(";")
        assert any("5 of 5 bootstrap replicates" in w for w in warnings)

    @pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
    def test_invalid_max_distance_raises(self, bad):
        with pytest.raises(ValueError, match="max_distance"):
            nj(msa(*self.SEQS), substitution_model="JukesCantor", max_distance=bad)

    def test_standalone_distance_matrix_is_raw(self):
        result = distance_matrix(msa(*self.SEQS), substitution_model="JukesCantor")
        assert math.isinf(result["matrix"][0][1])


# 190 equal exchangeabilities and 20 equal frequencies: the 20-state equal-rates model.
EQUAL_RATES_MATRIX = (
    "\n".join(" ".join(["1"] * i) for i in range(1, 20)) + "\n\n" + " ".join(["0.05"] * 20) + "\n"
)

PROTEINS = [
    {"identifier": "a", "sequence": "ARNDCQEGHILKMFPSTWYVARNDCQEGHILKMFPSTWYV"},
    {"identifier": "b", "sequence": "ARNDCQEGHILKMFPSTWYVARNDCQEGHILKMFPSTWYA"},
    {"identifier": "c", "sequence": "ARNDCQEGHILKMFPSTWYAARNDCQEGHILKMFPSTAYA"},
]


class TestCustomMatrix:
    def test_equal_rates_matrix_matches_the_closed_form(self):
        # Two sequences of 40 columns that differ in 6 of them.
        a = "ARNDCQEGHILKMFPSTWYV" * 2
        b = "WWWWWW" + a[6:]
        result = distance_matrix(
            [{"identifier": "a", "sequence": a}, {"identifier": "b", "sequence": b}],
            substitution_model="Custom",
            matrix=EQUAL_RATES_MATRIX,
        )
        p = 6 / 40
        want = -(19 / 20) * math.log(1 - p / (19 / 20))
        assert result["matrix"][0][1] == pytest.approx(want, rel=1e-8)

    def test_nj_and_average_distance_run(self):
        tree = nj(PROTEINS, substitution_model="Custom", matrix=EQUAL_RATES_MATRIX, n_bootstrap_samples=10)
        assert tree["newick"].endswith(";")
        for name in "abc":
            assert f"{name}:" in tree["newick"]
        mean = average_distance(PROTEINS, substitution_model="Custom", matrix=EQUAL_RATES_MATRIX)
        assert mean > 0

    def test_missing_matrix_raises(self):
        with pytest.raises(ValueError, match="InvalidMatrix"):
            nj(PROTEINS, substitution_model="Custom")

    def test_matrix_with_another_model_raises(self):
        with pytest.raises(ValueError, match="InvalidMatrix"):
            nj(PROTEINS, substitution_model="Poisson", matrix=EQUAL_RATES_MATRIX)

    def test_malformed_matrix_raises(self):
        with pytest.raises(ValueError, match="Invalid substitution matrix"):
            nj(PROTEINS, substitution_model="Custom", matrix="1 2 3")

    def test_custom_on_dna_raises(self):
        dna = [{"identifier": "a", "sequence": "ACGTACGT"}, {"identifier": "b", "sequence": "ACGTACGA"}]
        with pytest.raises(ValueError):
            nj(dna, substitution_model="Custom", matrix=EQUAL_RATES_MATRIX)


class TestBuiltinMatrices:
    @pytest.mark.parametrize("model", ["Lg", "Wag", "Jtt", "JttDcmut", "Dayhoff", "DayhoffDcmut"])
    def test_runs_and_is_finite(self, model):
        result = distance_matrix(PROTEINS, substitution_model=model)
        assert result["matrix"][0][1] > 0
        assert math.isfinite(result["matrix"][0][2])
        tree = nj(PROTEINS, substitution_model=model, n_bootstrap_samples=5)
        assert tree["newick"].endswith(";")

    def test_builtin_rejects_a_matrix(self):
        with pytest.raises(ValueError, match="InvalidMatrix"):
            nj(PROTEINS, substitution_model="Lg", matrix=EQUAL_RATES_MATRIX)

    def test_builtin_on_dna_raises(self):
        dna = [{"identifier": "a", "sequence": "ACGTACGT"}, {"identifier": "b", "sequence": "ACGTACGA"}]
        with pytest.raises(ValueError):
            nj(dna, substitution_model="Wag")


HKY = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]
DNA = [
    {"identifier": "a", "sequence": "ACGTACGTACGTAAGTACGTACGTACCTACGTACGT"},
    {"identifier": "b", "sequence": "ACGTACGTACGAAAGTACGTACGTACCTACGTACGA"},
    {"identifier": "c", "sequence": "ACGAACGTACGAAAGTACGTACGTACCTACGTTCGA"},
    {"identifier": "d", "sequence": "TCGAACGAACGAAAGTACTTACGTACCTACGTTCGA"},
]


class TestGtr:
    def test_equal_rates_and_frequencies_match_jukes_cantor(self):
        gtr = distance_matrix(DNA, substitution_model="Gtr", gtr_rates=[1] * 6, base_frequencies=[0.25] * 4)
        jc = distance_matrix(DNA, substitution_model="JukesCantor")
        for i in range(4):
            for j in range(4):
                assert gtr["matrix"][i][j] == pytest.approx(jc["matrix"][i][j], abs=1e-8)

    def test_pooled_frequencies_are_the_default(self):
        a = distance_matrix(DNA, substitution_model="Gtr", gtr_rates=HKY)
        b = distance_matrix(DNA, substitution_model="Gtr", gtr_rates=HKY, base_frequencies=[0.25] * 4)
        assert a["matrix"] != b["matrix"]
        assert all(math.isfinite(x) for row in a["matrix"] for x in row)

    def test_nj_with_bootstrap_and_rate_variation(self):
        tree = nj(DNA, substitution_model="Gtr", gtr_rates=HKY, gamma_shape=0.5, p_invar=0.1, n_bootstrap_samples=10)
        assert tree["newick"].endswith(";")
        assert average_distance(DNA, substitution_model="Gtr", gtr_rates=HKY) > 0

    def test_missing_rates_raise(self):
        with pytest.raises(ValueError, match="gtr_rates"):
            nj(DNA, substitution_model="Gtr")

    def test_wrong_number_of_rates_raises(self):
        with pytest.raises(ValueError, match="6 exchangeabilities"):
            nj(DNA, substitution_model="Gtr", gtr_rates=[1.0, 2.0])

    def test_rates_with_another_model_raise(self):
        with pytest.raises(ValueError, match="only used with the 'Gtr' model"):
            nj(DNA, substitution_model="Kimura2P", gtr_rates=HKY)

    def test_gtr_on_protein_raises(self):
        with pytest.raises(ValueError):
            nj(PROTEINS, substitution_model="Gtr", gtr_rates=HKY)
