const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const root = path.resolve(__dirname, "..");
const api = require(path.join(root, "app.js"));
const data = JSON.parse(fs.readFileSync(path.join(root, "peptide-cancer-data.json"), "utf8"));
const html = fs.readFileSync(path.join(root, "index.html"), "utf8");
const source = fs.readFileSync(path.join(root, "app.js"), "utf8");

assert.equal(data.study.samples.length, 11);
assert.equal(new Set(data.candidates.map(candidate => candidate.id)).size, data.candidates.length);
assert.ok(data.candidates.length > 0);
let combinations = 0;
for (const candidate of data.candidates) {
  const groups = api.sampleRows(candidate);
  assert.deepEqual(groups.map(group => group.n), [4, 4, 3], candidate.gene);
  for (const sample of candidate.samples) {
    assert.ok(Number.isFinite(sample.rawCount) && sample.rawCount >= 0);
    assert.ok(Number.isFinite(sample.normalizedCount) && sample.normalizedCount >= 0);
    const metadata = data.study.samples.find(item => item.id === sample.sampleId);
    assert.ok(metadata);
    assert.ok(Math.abs(sample.rawCount / metadata.sizeFactor - sample.normalizedCount) < 1e-6);
  }
  for (const peptide of candidate.peptides) {
    const evidence = api.peptideEvidence(peptide);
    if (evidence.comparable) {
      assert.equal(evidence.identicalPositions + evidence.mismatchPositions.length, evidence.human.length);
      assert.equal(evidence.sequenceConserved, Boolean(evidence.rawHuman &&
        evidence.rawHuman === evidence.rawMouse && !evidence.human.includes("-") && !evidence.mouse.includes("-")));
      assert.ok(evidence.mismatchPositions.every(position => position >= 1 && position <= evidence.human.length));
    }
    assert.equal(evidence.sequenceConserved, peptide.sequenceConserved, candidate.gene);
  }
  const defaultDecision = api.candidateDecision(candidate,
    { minEffect: .5, stableOnly: true, sequenceRule: "identical", requireCPM: false });
  assert.equal(defaultDecision.identical, candidate.evidence.conservedPeptideCount, candidate.gene);
  assert.equal(defaultDecision.conserved, candidate.evidence.conservedCleavagePeptideCount, candidate.gene);
  for (const minEffect of [0, .5, 1]) for (const stableOnly of [false, true])
    for (const sequenceRule of ["identical", "cleavage"]) for (const requireCPM of [false, true]) {
      const settings = { minEffect, stableOnly, sequenceRule, requireCPM };
      const result = api.candidateDecision(candidate, settings);
      assert.equal(result.eligible, result.measurable && result.available &&
        Math.abs(candidate.effect.log2FoldChange) >= minEffect &&
        (!stableOnly || result.stable) &&
        (!requireCPM || candidate.effect.normalizationDirectionAgrees === true) &&
        (sequenceRule === "identical" ? result.identical > 0 : result.conserved > 0));
      const proposal = api.followup(candidate, settings);
      assert.ok(proposal.title.includes(candidate.gene));
      assert.ok(proposal.text.length > 80);
      combinations++;
    }
}
for (const minEffect of [0, .5, 1]) for (const stableOnly of [false, true])
  for (const sequenceRule of ["identical", "cleavage"]) for (const requireCPM of [false, true]) {
    const settings = { minEffect, stableOnly, sequenceRule, requireCPM };
    const ranked = api.rankCandidates(data.candidates, settings);
    assert.ok(ranked.every(result => result.eligible));
    assert.equal(new Set(ranked.map(result => result.candidate.id)).size, ranked.length);
    const exportText = api.shortlistCSV(data, settings);
    assert.equal(exportText.trim().split("\r\n").length, ranked.length + 1);
    assert.ok(!exportText.includes("undefined"));
    if (ranked.length) assert.ok(exportText.includes("not measured peptide abundance"));
  }
assert.equal(api.mean([]), null);
assert.equal(api.mean([1, 2, 3]), 2);
assert.equal(api.fmt(null), "not available");
assert.equal(api.peptideEvidence(null).comparable, false);
assert.equal(api.peptideEvidence({ humanSequence: "ABC", mouseSequence: "ADC" }).identicalPositions, 2);
assert.deepEqual(api.peptideEvidence({ humanSequence: "ABC", mouseSequence: "ADC" }).mismatchPositions, [2]);
assert.equal(api.peptideEvidence({ humanSequence: "ABC", mouseSequence: "AB" }).comparable, false);
assert.equal(api.followup(null, {}).title, "No candidate selected");
assert.ok(api.niceAxis([0, 0], true).high > 0);
for (const values of [[-3, -.5, 1.5], [0, 10], []]) {
  const axis = api.niceAxis(values);
  assert.ok(Number.isFinite(axis.low) && Number.isFinite(axis.high) && axis.high > axis.low);
  assert.ok(values.every(value => value >= axis.low && value <= axis.high));
}

// Minimal DOM smoke test: run the real UI code with the real derived data,
// exercise every candidate and selection rule, and inspect SVG coordinates.
// Full browser layout/accessibility checks are performed separately.
class Element {
  constructor(tag) {
    this.tagName = tag; this.children = []; this.attrs = {}; this.listeners = {};
    this.dataset = {}; this.disabled = false; this.hidden = false; this.checked = false;
    this.textContent = ""; this._value = undefined;
  }
  append(...items) { this.children.push(...items); }
  replaceChildren(...items) { this.children = items; this._value = undefined; }
  setAttribute(key, value) { this.attrs[key] = String(value); }
  getAttribute(key) { return this.attrs[key]; }
  addEventListener(type, handler) { (this.listeners[type] ||= []).push(handler); }
  dispatch(type) { for (const handler of this.listeners[type] || []) handler({ target: this }); }
  get value() { return this._value ?? this.children[0]?.value ?? ""; }
  set value(value) { this._value = String(value); }
  scrollIntoView() {}
  remove() {}
  click() { this.dispatch("click"); }
}
const ids = new Map([...html.matchAll(/\bid="([^"]+)"/g)].map(match => [match[1], new Element("div")]));
ids.get("effect-select").value = ".5";
ids.get("sequence-rule").value = "identical";
ids.get("stable-only").checked = true;
const errors = [], downloads = [];
const document = {
  readyState: "complete", body: new Element("body"),
  getElementById: id => ids.get(id) || null,
  createElement: tag => new Element(tag),
  createElementNS: (_, tag) => new Element(tag),
  querySelectorAll: selector => selector === "#suggestions button" ?
    ids.get("suggestions").children.filter(element => element.tagName === "button") : []
};
const context = {
  document, window: { UNIVERSITY_LAB_DEMO_DATA: { candidateFile: "peptide-cancer-data.json" } },
  fetch: async () => ({ ok: true, json: async () => structuredClone(data) }),
  console: { error: error => errors.push(String(error)) },
  URL: { createObjectURL: blob => { downloads.push(blob); return "blob:test"; }, revokeObjectURL() {} },
  Blob: class { constructor(parts, options) { this.parts = parts; this.options = options; } },
  setTimeout: handler => { handler(); return 1; }
};
vm.runInNewContext(source, context);
function assertFiniteSVG(element) {
  for (const [key, value] of Object.entries(element.attrs)) {
    if (["x", "y", "x1", "y1", "x2", "y2", "cx", "cy", "r", "width", "height"].includes(key))
      assert.ok(Number.isFinite(Number(value)), `${key}=${value}`);
  }
  for (const child of element.children) if (child instanceof Element) assertFiniteSVG(child);
}
function descendants(element) {
  return [element, ...element.children.flatMap(child => child instanceof Element ? descendants(child) : [])];
}
(async () => {
  await new Promise(resolve => setImmediate(resolve));
  assert.deepEqual(errors, []);
  assert.match(ids.get("load-status").textContent, /loaded/);
  for (const candidate of data.candidates) {
    ids.get("candidate-select").value = candidate.id;
    ids.get("candidate-select").dispatch("change");
    assert.ok(ids.get("readout").textContent.includes(candidate.gene));
    assertFiniteSVG(ids.get("hero-viz"));
    assert.ok(ids.get("sequence-viz").children.length > 0);
    for (let peptide = 0; peptide < candidate.peptides.length; peptide++) {
      ids.get("peptide-select").value = peptide;
      ids.get("peptide-select").dispatch("change");
      assert.ok(!ids.get("sequence-caption").textContent.includes("undefined"));
    }
  }
  for (const minEffect of ["0", ".5", "1"]) for (const stableOnly of [false, true])
    for (const sequenceRule of ["identical", "cleavage"]) for (const requireCPM of [false, true]) {
      ids.get("effect-select").value = minEffect;
      ids.get("stable-only").checked = stableOnly;
      ids.get("sequence-rule").value = sequenceRule;
      ids.get("cpm-agreement").checked = requireCPM;
      ids.get("sequence-rule").dispatch("change");
      const ranked = api.rankCandidates(data.candidates, { minEffect: Number(minEffect), stableOnly, sequenceRule, requireCPM });
      assert.equal(Number(ids.get("metric-three").textContent), ranked.length);
      if (ranked.length) {
        ids.get("export-shortlist").click();
        assert.ok(downloads.at(-1).parts[0].includes("GSE266899"));
        const button = descendants(ids.get("shortlist")).find(element =>
          element.getAttribute("aria-label")?.startsWith("Compare "));
        button?.click();
        assertFiniteSVG(ids.get("comparison-viz"));
        ids.get("clear-comparison").click();
        assert.equal(ids.get("compare-region").hidden, true);
      }
      assertFiniteSVG(ids.get("hero-viz"));
    }
  assert.deepEqual(errors, []);
  console.log(`PASS: ${data.candidates.length} candidate views, ${combinations} candidate/rule combinations, sequence identities, shortlist CSV and real-data DOM/SVG smoke checks.`);
})().catch(error => { console.error(error); process.exitCode = 1; });
