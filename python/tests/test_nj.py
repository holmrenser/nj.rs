import pytest
from nj_py import nj


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

    def test_on_progress_called_correct_number_of_times(self):
        calls = []
        nj(
            msa(("A", "ACGTACGT"), ("B", "ACGAACGA"), ("C", "AGGTCGGT"), ("D", "TGGTCGGT")),
            n_bootstrap_samples=5,
            on_progress=lambda current, total: calls.append((current, total)),
        )
        assert len(calls) == 5
        assert calls[-1] == (5, 5)
