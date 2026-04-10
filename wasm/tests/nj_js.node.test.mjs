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

// --- distance_matrix ---

const makeDistConfig = (msa, substitutionModel = "PDiff") => ({
  msa: msa.map(([identifier, sequence]) => ({ identifier, sequence })),
  substitution_model: substitutionModel,
});

test("wasm module exposes distance_matrix export", async () => {
  const wasm = await loadWasm();
  assert.equal(typeof wasm.distance_matrix, "function");
});

test("distance_matrix returns object with names and matrix", async () => {
  const wasm = await loadWasm();
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
  assert.ok(Array.isArray(result.names), "names should be an array");
  assert.ok(Array.isArray(result.matrix), "matrix should be an array");
});

test("distance_matrix names match input identifiers", async () => {
  const wasm = await loadWasm();
  const result = wasm.distance_matrix(makeDistConfig([["human", "ACGT"], ["chimp", "ACGA"]]));
  assert.deepEqual(result.names, ["human", "chimp"]);
});

test("distance_matrix produces n×n matrix", async () => {
  const wasm = await loadWasm();
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]));
  assert.equal(result.matrix.length, 3);
  assert.ok(result.matrix.every((row) => row.length === 3));
});

test("distance_matrix diagonal is zero", async () => {
  const wasm = await loadWasm();
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]));
  for (let i = 0; i < 3; i++) {
    assert.equal(result.matrix[i][i], 0.0);
  }
});

test("distance_matrix is symmetric", async () => {
  const wasm = await loadWasm();
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]));
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      assert.equal(result.matrix[i][j], result.matrix[j][i]);
    }
  }
});

test("distance_matrix PDiff known value", async () => {
  const wasm = await loadWasm();
  // one difference (T vs A at pos 3) out of 4 → 0.25
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
  assert.ok(Math.abs(result.matrix[0][1] - 0.25) < 1e-12);
});

test("distance_matrix identical sequences give zero off-diagonal", async () => {
  const wasm = await loadWasm();
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGT"]]));
  assert.equal(result.matrix[0][1], 0.0);
});

test("distance_matrix JukesCantor known value", async () => {
  const wasm = await loadWasm();
  // p=0.25 → d = -0.75 * ln(1 - 4/3 * 0.25)
  const expected = -0.75 * Math.log(1.0 - (4.0 / 3.0) * 0.25);
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "JukesCantor"));
  assert.ok(Math.abs(result.matrix[0][1] - expected) < 1e-10);
});

test("distance_matrix Poisson protein known value", async () => {
  const wasm = await loadWasm();
  // 1 diff (H vs K) out of 7 → p=1/7, d=-ln(1-1/7)
  const expected = -Math.log(1.0 - 1.0 / 7.0);
  const result = wasm.distance_matrix(makeDistConfig([["A", "ACDEFGH"], ["B", "ACDEFGK"]], "Poisson"));
  assert.ok(Math.abs(result.matrix[0][1] - expected) < 1e-10);
});

test("distance_matrix empty MSA throws", async () => {
  const wasm = await loadWasm();
  assert.throws(() => wasm.distance_matrix({ msa: [], substitution_model: "PDiff" }));
});

test("distance_matrix incompatible model throws", async () => {
  const wasm = await loadWasm();
  assert.throws(() => wasm.distance_matrix(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "Poisson")));
});

// --- average_distance ---

test("wasm module exposes average_distance export", async () => {
  const wasm = await loadWasm();
  assert.equal(typeof wasm.average_distance, "function");
});

test("average_distance returns a number", async () => {
  const wasm = await loadWasm();
  const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
  assert.equal(typeof result, "number");
});

test("average_distance identical sequences returns zero", async () => {
  const wasm = await loadWasm();
  const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGT"]]));
  assert.equal(result, 0.0);
});

test("average_distance two taxa equals pairwise PDiff", async () => {
  const wasm = await loadWasm();
  // one difference out of 4 → 0.25
  const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]]));
  assert.ok(Math.abs(result - 0.25) < 1e-12);
});

test("average_distance three taxa known value", async () => {
  const wasm = await loadWasm();
  // A↔B: 0.25 (T→A), A↔C: 0.25 (C→G), B↔C: 0.5 → avg = 1/3
  const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"], ["C", "AGGT"]]));
  assert.ok(Math.abs(result - 1.0 / 3.0) < 1e-12);
});

test("average_distance JukesCantor known value", async () => {
  const wasm = await loadWasm();
  const expected = -0.75 * Math.log(1.0 - (4.0 / 3.0) * 0.25);
  const result = wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "JukesCantor"));
  assert.ok(Math.abs(result - expected) < 1e-10);
});

test("average_distance empty MSA throws", async () => {
  const wasm = await loadWasm();
  assert.throws(() => wasm.average_distance({ msa: [], substitution_model: "PDiff" }));
});

test("average_distance incompatible model throws", async () => {
  const wasm = await loadWasm();
  assert.throws(() => wasm.average_distance(makeDistConfig([["A", "ACGT"], ["B", "ACGA"]], "Poisson")));
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