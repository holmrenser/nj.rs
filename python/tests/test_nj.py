import pytest
import nj_py


def make_config(msa, substitution_model="PDiff", n_bootstrap_samples=0):
    return {
        "msa": [{"identifier": name, "sequence": seq} for name, seq in msa],
        "substitution_model": substitution_model,
        "n_bootstrap_samples": n_bootstrap_samples,
    }


class TestBasic:
    def test_returns_string(self):
        config = make_config([("A", "ACGT"), ("B", "ACGA")])
        result = nj_py.nj(config)
        assert isinstance(result, str)

    def test_newick_format(self):
        config = make_config([("A", "ACGT"), ("B", "ACGA")])
        result = nj_py.nj(config)
        assert result.startswith("(")
        assert result.endswith(";")

    def test_two_identical_taxa_zero_distance(self):
        config = make_config([("A", "ACGT"), ("B", "ACGT")])
        result = nj_py.nj(config)
        assert result == "(A:0.000,B:0.000);"

    def test_three_taxa_contains_all_leaves(self):
        config = make_config([("Sp1", "ACGT"), ("Sp2", "ACGA"), ("Sp3", "AGGT")])
        result = nj_py.nj(config)
        assert "Sp1" in result
        assert "Sp2" in result
        assert "Sp3" in result

    def test_deterministic(self):
        config = make_config([("A", "ACGTCG"), ("B", "ACG-GC"), ("C", "ACGCGT")])
        assert nj_py.nj(config) == nj_py.nj(config)


class TestSubstitutionModels:
    def test_jukes_cantor_dna(self):
        config = make_config(
            [("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")],
            substitution_model="JukesCantor",
        )
        result = nj_py.nj(config)
        assert result.endswith(";")

    def test_kimura2p_dna(self):
        config = make_config(
            [("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")],
            substitution_model="Kimura2P",
        )
        result = nj_py.nj(config)
        assert result.endswith(";")

    def test_pdiff_protein(self):
        config = make_config(
            [("A", "ACDEFGH"), ("B", "ACDEFGK"), ("C", "ACDLFGH")],
            substitution_model="PDiff",
        )
        result = nj_py.nj(config)
        assert result.endswith(";")

    def test_poisson_protein(self):
        config = make_config(
            [("A", "ACDEFGH"), ("B", "ACDEFGK"), ("C", "ACDLFGH")],
            substitution_model="Poisson",
        )
        result = nj_py.nj(config)
        assert result.endswith(";")


class TestErrors:
    def test_empty_msa_raises(self):
        config = {"msa": [], "substitution_model": "PDiff", "n_bootstrap_samples": 0}
        with pytest.raises(ValueError, match="[Ee]mpty"):
            nj_py.nj(config)

    def test_poisson_on_dna_raises(self):
        config = make_config([("A", "ACGT"), ("B", "ACGA")], substitution_model="Poisson")
        with pytest.raises(ValueError):
            nj_py.nj(config)

    def test_jukes_cantor_on_protein_raises(self):
        config = make_config(
            [("A", "ACDEFGH"), ("B", "ACDEFGK")],
            substitution_model="JukesCantor",
        )
        with pytest.raises(ValueError):
            nj_py.nj(config)

    def test_kimura2p_on_protein_raises(self):
        config = make_config(
            [("A", "ACDEFGH"), ("B", "ACDEFGK")],
            substitution_model="Kimura2P",
        )
        with pytest.raises(ValueError):
            nj_py.nj(config)


class TestBootstrap:
    def test_bootstrap_zero_samples_no_support_values(self):
        config = make_config(
            [("A", "ACGTACGT"), ("B", "ACGAACGA"), ("C", "AGGTCGGT"), ("D", "TGGTCGGT")],
            n_bootstrap_samples=0,
        )
        result = nj_py.nj(config)
        assert result.endswith(";")

    def test_bootstrap_nonzero_returns_valid_newick(self):
        config = make_config(
            [("A", "ACGTACGT"), ("B", "ACGAACGA"), ("C", "AGGTCGGT"), ("D", "TGGTCGGT")],
            n_bootstrap_samples=10,
        )
        result = nj_py.nj(config)
        assert isinstance(result, str)
        assert result.endswith(";")