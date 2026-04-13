import math

import pytest
from nj_py import average_distance, distance_matrix, nj


def msa(*pairs):
    return [{"identifier": name, "sequence": seq} for name, seq in pairs]


class TestBasic:
    def test_returns_string(self):
        result = nj(msa(("A", "ACGT"), ("B", "ACGA")))
        assert isinstance(result, str)

    def test_newick_format(self):
        result = nj(msa(("A", "ACGT"), ("B", "ACGA")))
        assert result.startswith("(")
        assert result.endswith(";")

    def test_two_identical_taxa_zero_distance(self):
        result = nj(msa(("A", "ACGT"), ("B", "ACGT")))
        assert result == "(A:0.000,B:0.000);"

    def test_three_taxa_contains_all_leaves(self):
        result = nj(msa(("Sp1", "ACGT"), ("Sp2", "ACGA"), ("Sp3", "AGGT")))
        assert "Sp1" in result
        assert "Sp2" in result
        assert "Sp3" in result

    def test_deterministic(self):
        sequences = msa(("A", "ACGTCG"), ("B", "ACG-GC"), ("C", "ACGCGT"))
        assert nj(sequences) == nj(sequences)


class TestSubstitutionModels:
    def test_jukes_cantor_dna(self):
        result = nj(
            msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")),
            substitution_model="JukesCantor",
        )
        assert result.endswith(";")

    def test_kimura2p_dna(self):
        result = nj(
            msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")),
            substitution_model="Kimura2P",
        )
        assert result.endswith(";")

    def test_pdiff_protein(self):
        result = nj(
            msa(("A", "ACDEFGH"), ("B", "ACDEFGK"), ("C", "ACDLFGH")),
            substitution_model="PDiff",
        )
        assert result.endswith(";")

    def test_poisson_protein(self):
        result = nj(
            msa(("A", "ACDEFGH"), ("B", "ACDEFGK"), ("C", "ACDLFGH")),
            substitution_model="Poisson",
        )
        assert result.endswith(";")


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
        assert result.endswith(";")

    def test_mixedcase_dna(self):
        result = nj(msa(("A", "AcGt"), ("B", "aCgA")))
        assert result.endswith(";")

    def test_lowercase_matches_uppercase(self):
        upper = nj(msa(("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")))
        lower = nj(msa(("A", "acgt"), ("B", "acga"), ("C", "aggt")))
        assert upper == lower


class TestSpecialNames:
    def test_sequence_name_with_spaces(self):
        result = nj(msa(("Seq A", "ACGT"), ("Seq B", "ACGA")))
        assert "Seq A" in result
        assert "Seq B" in result

    def test_sequence_name_with_special_chars(self):
        result = nj(msa(("Sp|1.1", "ACGT"), ("Sp|2.1", "ACGA")))
        assert "Sp|1.1" in result
        assert "Sp|2.1" in result


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
        assert result.endswith(";")

    def test_bootstrap_nonzero_returns_valid_newick(self):
        result = nj(
            msa(("A", "ACGTACGT"), ("B", "ACGAACGA"), ("C", "AGGTCGGT"), ("D", "TGGTCGGT")),
            n_bootstrap_samples=10,
        )
        assert isinstance(result, str)
        assert result.endswith(";")

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
