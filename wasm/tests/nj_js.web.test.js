import { expect } from "@esm-bundle/chai";

const loadWebBundle = async () => {
  const mod = await import("../pkg-web/nj_js.js");
  return mod.default ?? mod;
};

describe("nj_js web bundle", () => {
  it("exposes nj export", async () => {
    const wasm = await loadWebBundle();
    expect(wasm).to.have.property("nj").that.is.a("function");
  });

  it("computes Newick tree from MSA FASTA", async () => {
    const wasm = await loadWebBundle();

    const msa = [
      { header: "seq1", sequence: "ACGT" },
      { header: "seq2", sequence: "ACGA" },
      { header: "seq3", sequence: "AGGT" },
    ];

    const config = {
      msa,
      hide_internal: true,
    };

    const newick = wasm.nj(config);
    expect(newick).to.be.a("string");
    expect(newick).to.equal("((seq2:0.250,seq1:0.000):0.125,seq3:0.125);");
  });
});
