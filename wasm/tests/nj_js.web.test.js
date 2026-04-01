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

describe("nj_js web bundle — module shape", () => {
  it("exposes nj export", async () => {
    const wasm = await loadWebBundle();
    expect(wasm).to.have.property("nj").that.is.a("function");
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
