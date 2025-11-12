import test from "node:test";
import assert from "node:assert/strict";

test("wasm module loads and exposes nj export", async () => {
  const moduleNs = await import("../pkg-node/nj_js.js");
  const wasm = moduleNs.default ?? moduleNs;

  assert.equal(typeof wasm.nj, "function", "nj export missing");
});

test("nj function works", async () => {
  const moduleNs = await import("../pkg-node/nj_js.js");
  const wasm = moduleNs.default ?? moduleNs;

  const msa = [
    { header: "seq1", sequence: "ACGT" },
    { header: "seq2", sequence: "ACGA" },
    { header: "seq3", sequence: "AGGT" },
  ];

  const config = {
    msa,
    hide_internal: true,
  };

  const newick = await wasm.nj(config);
  assert.equal(
    newick,
    "((seq2:0.250,seq1:0.000):0.125,seq3:0.125);",
    "Unexpected Newick output"
  );
});
