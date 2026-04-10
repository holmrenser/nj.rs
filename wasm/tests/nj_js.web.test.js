import { expect } from "@esm-bundle/chai";

const loadWebBundle = async () => {
  const mod = await import("../pkg-web/nj_js.js");
  return mod.default ?? mod;
};

const makeConfig = (
  msa,
  substitutionModel = "PDiff",
  nBootstrapSamples = 0,
) => ({
  msa: msa.map(([identifier, sequence]) => ({ identifier, sequence })),
  substitution_model: substitutionModel,
  n_bootstrap_samples: nBootstrapSamples,
});

const makeDistConfig = (msa, substitutionModel = "PDiff") => ({
  msa: msa.map(([identifier, sequence]) => ({ identifier, sequence })),
  substitution_model: substitutionModel,
});

describe("nj_js web bundle — module shape", () => {
  it("exposes nj export", async () => {
    const wasm = await loadWebBundle();
    expect(wasm).to.have.property("nj").that.is.a("function");
  });

  it("exposes distance_matrix export", async () => {
    const wasm = await loadWebBundle();
    expect(wasm).to.have.property("distance_matrix").that.is.a("function");
  });

  it("exposes average_distance export", async () => {
    const wasm = await loadWebBundle();
    expect(wasm).to.have.property("average_distance").that.is.a("function");
  });
});

describe("nj_js web bundle — basic functionality", () => {
  it("two identical taxa produce zero-distance newick", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig([
      ["A", "ACGT"],
      ["B", "ACGT"],
    ]);
    const newick = wasm.nj(config);
    expect(newick).to.equal("(A:0.000,B:0.000);");
  });

  it("three taxa newick contains all leaf names", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig([
      ["Sp1", "ACGT"],
      ["Sp2", "ACGA"],
      ["Sp3", "AGGT"],
    ]);
    const newick = wasm.nj(config);
    expect(newick).to.include("Sp1");
    expect(newick).to.include("Sp2");
    expect(newick).to.include("Sp3");
    expect(newick).to.match(/^\(/);
    expect(newick).to.match(/;$/);
  });

  it("output is deterministic across calls", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig([
      ["A", "ACGTCG"],
      ["B", "ACG-GC"],
      ["C", "ACGCGT"],
    ]);
    expect(wasm.nj(config)).to.equal(wasm.nj(config));
  });
});

describe("nj_js web bundle — substitution models", () => {
  it("JukesCantor model runs on DNA", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig(
      [
        ["A", "ACGT"],
        ["B", "ACGA"],
        ["C", "AGGT"],
      ],
      "JukesCantor",
    );
    const newick = wasm.nj(config);
    expect(newick).to.be.a("string").and.match(/;$/);
  });

  it("Kimura2P model runs on DNA", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig(
      [
        ["A", "ACGT"],
        ["B", "ACGA"],
        ["C", "AGGT"],
      ],
      "Kimura2P",
    );
    const newick = wasm.nj(config);
    expect(newick).to.be.a("string").and.match(/;$/);
  });

  it("Poisson model runs on protein sequences", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig(
      [
        ["A", "ACDEFGH"],
        ["B", "ACDEFGK"],
        ["C", "ACDLFGH"],
      ],
      "Poisson",
    );
    const newick = wasm.nj(config);
    expect(newick).to.be.a("string").and.match(/;$/);
  });

  it("PDiff model runs on protein sequences", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig(
      [
        ["A", "ACDEFGH"],
        ["B", "ACDEFGK"],
        ["C", "ACDLFGH"],
      ],
      "PDiff",
    );
    const newick = wasm.nj(config);
    expect(newick).to.be.a("string").and.match(/;$/);
  });
});

describe("nj_js web bundle — error handling", () => {
  it("empty MSA throws", async () => {
    const wasm = await loadWebBundle();
    const config = {
      msa: [],
      substitution_model: "PDiff",
      n_bootstrap_samples: 0,
    };
    expect(() => wasm.nj(config)).to.throw();
  });

  it("Poisson model with DNA throws", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig(
      [
        ["A", "ACGT"],
        ["B", "ACGA"],
      ],
      "Poisson",
    );
    expect(() => wasm.nj(config)).to.throw();
  });

  it("JukesCantor model with protein throws", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig(
      [
        ["A", "ACDEFGH"],
        ["B", "ACDEFGK"],
      ],
      "JukesCantor",
    );
    expect(() => wasm.nj(config)).to.throw();
  });
});

describe("nj_js web bundle — distance_matrix", () => {
  it("returns object with names and matrix", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
    expect(result).to.have.property("names").that.is.an("array");
    expect(result).to.have.property("matrix").that.is.an("array");
  });

  it("names match input identifiers", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.distance_matrix(makeDistConfig([["human", "ACGT"], ["chimp", "ACGA"]]));
    expect(result.names).to.deep.equal(["human", "chimp"]);
  });

  it("produces n×n matrix", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]));
    expect(result.matrix).to.have.length(3);
    result.matrix.forEach((row) => expect(row).to.have.length(3));
  });

  it("diagonal is zero", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]));
    for (let i = 0; i < 3; i++) {
      expect(result.matrix[i][i]).to.equal(0.0);
    }
  });

  it("matrix is symmetric", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]));
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        expect(result.matrix[i][j]).to.equal(result.matrix[j][i]);
      }
    }
  });

  it("PDiff known value", async () => {
    const wasm = await loadWebBundle();
    // one difference (T vs A at pos 3) out of 4 → 0.25
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
    expect(Math.abs(result.matrix[0][1] - 0.25)).to.be.lessThan(1e-12);
  });

  it("identical sequences give zero off-diagonal", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGT"]]));
    expect(result.matrix[0][1]).to.equal(0.0);
  });

  it("JukesCantor known value", async () => {
    const wasm = await loadWebBundle();
    const expected = -0.75 * Math.log(1.0 - (4.0 / 3.0) * 0.25);
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "JukesCantor"));
    expect(Math.abs(result.matrix[0][1] - expected)).to.be.lessThan(1e-10);
  });

  it("Poisson protein known value", async () => {
    const wasm = await loadWebBundle();
    const expected = -Math.log(1.0 - 1.0 / 7.0);
    const result = wasm.distance_matrix(makeDistConfig([["A", "ACDEFGH"], ["B", "ACDEFGK"]], "Poisson"));
    expect(Math.abs(result.matrix[0][1] - expected)).to.be.lessThan(1e-10);
  });

  it("empty MSA throws", async () => {
    const wasm = await loadWebBundle();
    expect(() => wasm.distance_matrix({ msa: [], substitution_model: "PDiff" })).to.throw();
  });

  it("incompatible model throws", async () => {
    const wasm = await loadWebBundle();
    expect(() =>
      wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "Poisson")),
    ).to.throw();
  });
});

describe("nj_js web bundle — average_distance", () => {
  it("returns a number", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
    expect(result).to.be.a("number");
  });

  it("identical sequences returns zero", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGT"]]));
    expect(result).to.equal(0.0);
  });

  it("two taxa equals pairwise PDiff", async () => {
    const wasm = await loadWebBundle();
    const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
    expect(Math.abs(result - 0.25)).to.be.lessThan(1e-12);
  });

  it("three taxa known value", async () => {
    const wasm = await loadWebBundle();
    // A↔B: 0.25 (T→A), A↔C: 0.25 (C→G), B↔C: 0.5 → avg = 1/3
    const result = wasm.average_distance(
      makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]),
    );
    expect(Math.abs(result - 1.0 / 3.0)).to.be.lessThan(1e-12);
  });

  it("JukesCantor known value", async () => {
    const wasm = await loadWebBundle();
    const expected = -0.75 * Math.log(1.0 - (4.0 / 3.0) * 0.25);
    const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "JukesCantor"));
    expect(Math.abs(result - expected)).to.be.lessThan(1e-10);
  });

  it("empty MSA throws", async () => {
    const wasm = await loadWebBundle();
    expect(() => wasm.average_distance({ msa: [], substitution_model: "PDiff" })).to.throw();
  });

  it("incompatible model throws", async () => {
    const wasm = await loadWebBundle();
    expect(() =>
      wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "Poisson")),
    ).to.throw();
  });
});

describe("nj_js web bundle — bootstrap", () => {
  it("bootstrap with n>0 returns valid newick", async () => {
    const wasm = await loadWebBundle();
    const config = makeConfig(
      [
        ["A", "ACGTACGT"],
        ["B", "ACGAACGA"],
        ["C", "AGGTCGGT"],
        ["D", "TGGTCGGT"],
      ],
      "PDiff",
      10,
    );
    const newick = wasm.nj(config);
    expect(newick).to.be.a("string").and.match(/;$/);
  });
});
