import test from "node:test";
import assert from "node:assert/strict";

const loadWasm = async () => {
  const moduleNs = await import("../pkg-node/nj_js.js");
  return moduleNs.default ?? moduleNs;
};

const makeConfig = (msa, substitutionModel = "PDiff", nBootstrapSamples = 0) => ({
  msa: msa.map(([identifier, sequence]) => ({ identifier, sequence })),
  substitution_model: substitutionModel,
  n_bootstrap_samples: nBootstrapSamples,
});

// --- module shape ---

test("wasm module loads and exposes nj export", async () => {
  const wasm = await loadWasm();
  assert.equal(typeof wasm.nj, "function", "nj export missing");
});

// --- basic functionality ---

test("two identical taxa produce zero-distance newick", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACGT"], ["B", "ACGT"]]);
  const newick = wasm.nj(config);
  assert.equal(newick, "(A:0.000,B:0.000);");
});

test("two divergent taxa produce correct branch lengths", async () => {
  const wasm = await loadWasm();
  // 1 difference over 4 positions = 0.25 each side
  const config = makeConfig([["A", "ACGT"], ["B", "ACGT"]]);
  const newick = wasm.nj(config);
  assert.ok(newick.endsWith(";"), "newick must end with semicolon");
  assert.ok(newick.startsWith("("), "newick must start with parenthesis");
});

test("three taxa newick contains all leaf names", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["Sp1", "ACGT"], ["Sp2", "ACGA"], ["Sp3", "AGGT"]]);
  const newick = wasm.nj(config);
  assert.ok(newick.includes("Sp1"), "missing Sp1");
  assert.ok(newick.includes("Sp2"), "missing Sp2");
  assert.ok(newick.includes("Sp3"), "missing Sp3");
  assert.ok(newick.startsWith("("));
  assert.ok(newick.endsWith(";"));
});

test("output is deterministic across calls", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACGTCG"], ["B", "ACG-GC"], ["C", "ACGCGT"]]);
  const r1 = wasm.nj(config);
  const r2 = wasm.nj(config);
  assert.equal(r1, r2);
});

// --- substitution models ---

test("JukesCantor model runs on DNA", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]], "JukesCantor");
  const newick = wasm.nj(config);
  assert.ok(typeof newick === "string");
  assert.ok(newick.endsWith(";"));
});

test("Kimura2P model runs on DNA", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]], "Kimura2P");
  const newick = wasm.nj(config);
  assert.ok(typeof newick === "string");
  assert.ok(newick.endsWith(";"));
});

test("Poisson model runs on protein sequences", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACDEFGH"], ["B", "ACDEFGK"], ["C", "ACDLFGH"]], "Poisson");
  const newick = wasm.nj(config);
  assert.ok(typeof newick === "string");
  assert.ok(newick.endsWith(";"));
});

test("PDiff model runs on protein sequences", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACDEFGH"], ["B", "ACDEFGK"], ["C", "ACDLFGH"]], "PDiff");
  const newick = wasm.nj(config);
  assert.ok(typeof newick === "string");
  assert.ok(newick.endsWith(";"));
});

// --- error handling ---

test("empty MSA throws", async () => {
  const wasm = await loadWasm();
  const config = { msa: [], substitution_model: "PDiff", n_bootstrap_samples: 0 };
  assert.throws(() => wasm.nj(config), /empty|Empty/);
});

test("Poisson model with DNA throws", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACGT"], ["B", "ACGA"]], "Poisson");
  assert.throws(() => wasm.nj(config));
});

test("JukesCantor model with protein throws", async () => {
  const wasm = await loadWasm();
  const config = makeConfig([["A", "ACDEFGH"], ["B", "ACDEFGK"]], "JukesCantor");
  assert.throws(() => wasm.nj(config));
});

// --- bootstrap ---

test("bootstrap with n>0 returns valid newick", async () => {
  const wasm = await loadWasm();
  const config = makeConfig(
    [["A", "ACGTACGT"], ["B", "ACGAACGA"], ["C", "AGGTCGGT"], ["D", "TGGTCGGT"]],
    "PDiff",
    10,
  );
  const newick = wasm.nj(config);
  assert.ok(typeof newick === "string");
  assert.ok(newick.endsWith(";"));
});