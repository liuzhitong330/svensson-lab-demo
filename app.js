(function () {
  "use strict";
  const finite = value => typeof value === "number" && Number.isFinite(value);
  const fmt = (value, digits = 2) => finite(value) ? value.toFixed(digits) : "not available";
  const signed = value => finite(value) ? `${value > 0 ? "+" : ""}${fmt(value)}` : "not available";
  const mean = values => values.length && values.every(finite) ?
    values.reduce((sum, value) => sum + value, 0) / values.length : null;
  const asText = value => typeof value === "string" ? value : "";
  const csv = rows => rows.map(row => row.map(value =>
    `"${String(value ?? "").replace(/"/g, '""')}"`).join(",")).join("\r\n") + "\r\n";
  const GROUPS = [
    { id: "LFD", label: "Low-fat diet" },
    { id: "HFD", label: "High-fat diet" },
    { id: "SEMA", label: "HFD + semaglutide" }
  ];

  function niceAxis(values, startAtZero = false) {
    const numbers = values.filter(finite);
    const minimum = numbers.length ? Math.min(0, ...numbers) : 0;
    const maximum = numbers.length ? Math.max(0, ...numbers) : 1;
    const power = 10 ** Math.floor(Math.log10((maximum - minimum || 1) / 4));
    const step = [1, 2, 2.5, 5, 10].find(value => (maximum - minimum || 1) / 4 / power <= value) * power;
    const low = startAtZero ? 0 : Math.floor(minimum / step) * step;
    const high = Math.max(low + step, Math.ceil(maximum / step) * step);
    return { low, high, ticks: Array.from({ length: Math.round((high - low) / step) + 1 },
      (_, index) => Number((low + step * index).toPrecision(12))) };
  }

  function peptideEvidence(peptide) {
    const rawHuman = asText(peptide?.humanSequence), rawMouse = asText(peptide?.mouseSequence);
    const human = asText(peptide?.alignedHuman) || rawHuman;
    const mouse = asText(peptide?.alignedMouse) || rawMouse;
    // Use the pipeline's projected pairwise alignment, retaining insertions/deletions.
    const comparable = Boolean(human && mouse && human.length === mouse.length);
    const mismatchPositions = comparable ? Array.from(human).flatMap((aa, index) =>
      aa === mouse[index] ? [] : [index + 1]) : [];
    return {
      human, mouse, rawHuman, rawMouse, comparable, mismatchPositions,
      identicalPositions: comparable ? Array.from(human).filter((aa, index) => aa === mouse[index] && aa !== "-").length : null,
      cleavageConserved: peptide?.cleavageConserved === true,
      sequenceConserved: Boolean(rawHuman && rawHuman === rawMouse && !human.includes("-") && !mouse.includes("-"))
    };
  }

  function candidateDecision(candidate, settings) {
    const effect = candidate?.effect || {};
    const measurable = candidate?.evidence?.measurable === true;
    const conserved = (candidate?.peptides || []).filter(peptide =>
      peptide.cleavageConserved === true && asText(peptide.mouseSequence)).length;
    const identical = (candidate?.peptides || []).filter(peptide =>
      peptide.cleavageConserved === true && peptideEvidence(peptide).sequenceConserved).length;
    const sequenceEligible = settings.sequenceRule === "cleavage" ? conserved > 0 : identical > 0;
    const available = finite(effect.log2FoldChange);
    const passesMagnitude = available && Math.abs(effect.log2FoldChange) >= settings.minEffect;
    const stable = effect.stableDirection === true;
    const reasons = [];
    if (!measurable) reasons.push("insufficient precursor RNA coverage");
    if (!conserved) reasons.push("no mapped peptide with predicted cleavage in both species");
    else if (!sequenceEligible) reasons.push("no identical peptide with both-species cleavage support");
    if (!available) reasons.push("RNA contrast unavailable");
    else if (!passesMagnitude) reasons.push("below the selected contrast floor");
    if (settings.stableOnly && !stable) reasons.push("direction depends on an individual animal");
    if (settings.requireCPM && effect.normalizationDirectionAgrees !== true)
      reasons.push("direction changes under counts-per-million normalization");
    return { candidate, measurable, conserved, identical, stable, available, passesMagnitude,
      eligible: reasons.length === 0, reasons };
  }

  function rankCandidates(candidates, settings) {
    return candidates.map(candidate => candidateDecision(candidate, settings))
      .filter(result => result.eligible)
      .sort((a, b) => Number(b.stable) - Number(a.stable) ||
        Math.abs(b.candidate.effect.log2FoldChange) - Math.abs(a.candidate.effect.log2FoldChange) ||
        a.candidate.id.localeCompare(b.candidate.id));
  }

  function sampleRows(candidate) {
    return GROUPS.map(group => {
      const samples = (candidate?.samples || []).filter(sample => sample.group === group.id);
      return { ...group, samples, n: samples.length,
        mean: mean(samples.map(sample => sample.normalizedCount)) };
    });
  }

  function followup(candidate, settings) {
    if (!candidate) return { title: "No candidate selected", text: "Load verified evidence before choosing an experiment." };
    const decision = candidateDecision(candidate, settings);
    const name = candidate.gene;
    if (!decision.measurable) return {
      title: `${name}: establish the RNA source first`,
      text: `This precursor does not meet the dataset’s expression-coverage rule. Check its expression and likely producing cells in independent material before committing to peptide detection; a large contrast at low coverage is not reliable prioritization evidence.`
    };
    if (!decision.conserved) return {
      title: `${name}: resolve the mouse assay sequence first`,
      text: `The current mapping does not provide a mouse peptide with the required predicted cleavage context. Verify precursor identity and processing before choosing a synthetic sequence; human precursor orthology alone is not a mouse peptide assay design.`
    };
    if (!decision.stable) return {
      title: `${name}: replicate the precursor signal before expanding`,
      text: `Omitting one animal changes the direction of the RNA contrast. Preserve the sequence candidate, but prioritize independent precursor measurements and cellular-source checks before a larger bioactivity screen.`
    };
    if (settings.sequenceRule !== "cleavage" && !decision.identical) return {
      title: `${name}: test the species-specific sequence explicitly`,
      text: `A mapped mouse sequence has predicted cleavage support, but it is not identical to the human candidate. It is excluded by the current sequence rule. A species-matched detection assay, followed by a direct human-versus-mouse comparison if warranted, would resolve the sequence choice.`
    };
    if (settings.requireCPM && candidate.effect.normalizationDirectionAgrees !== true) return {
      title: `${name}: resolve normalization sensitivity first`,
      text: `The precursor’s contrast changes direction under the counts-per-million check. Inspect library composition and independent precursor measurements before using this contrast to prioritize a peptide experiment. The alternative normalization is a diagnostic, not a second significance test.`
    };
    if (!decision.passesMagnitude) return {
      title: `${name}: consider a comparison candidate, not a selected lead`,
      text: `This precursor is below the chosen RNA-contrast floor. It may be useful as an expression comparison alongside a shortlisted precursor, but is not a validated inactive peptide control.`
    };
    return {
      title: `${name}: look for the endogenous peptide before testing function`,
      text: `Prioritize the mapped mouse sequence for targeted peptide detection in independent, carefully annotated material. Confirm identity against an appropriate standard and document the producing cells. A detected peptide would justify controlled cellular follow-up; RNA and cleavage prediction alone do not.`
    };
  }

  function shortlistCSV(data, settings, selectedIds = []) {
    const results = rankCandidates(data.candidates || [], settings);
    const selected = new Set(selectedIds);
    return csv([
      ["Rank", "Mouse gene", "Mouse gene ID", "Human gene", "Human UniProt",
        "Eligible human peptide sequences @ 1-based starts", "Mapped mouse peptide sequences @ 1-based starts",
        "Descriptive log2 contrast SEMA vs HFD", "Leave-one-out minimum", "Leave-one-out maximum",
        "Stable direction", "Mapped peptides with predicted cleavage in both species", "RNA coverage",
        "CPM log2 contrast", "Normalization directions agree", "Require CPM agreement",
        "Minimum absolute contrast", "Require stable direction", "Sequence rule", "Selected for comparison",
        "Next measurement", "Source", "Interpretation limit"],
      ...results.map((result, index) => {
        const candidate = result.candidate;
        const peptides = candidate.peptides.filter(peptide => peptide.cleavageConserved === true &&
          (settings.sequenceRule === "cleavage" || peptideEvidence(peptide).sequenceConserved));
        return [index + 1, candidate.gene, candidate.id, candidate.humanGene, candidate.uniprot,
          peptides.map(peptide => `${peptide.humanSequence}@${peptide.humanStart}`).join(" | "),
          peptides.map(peptide => `${peptide.mouseSequence}@${peptide.mouseStart}`).join(" | "),
          candidate.effect.log2FoldChange, candidate.effect.looMin, candidate.effect.looMax,
          result.stable, result.conserved, candidate.evidence.detectedN,
          candidate.effect.cpmLog2FoldChange, candidate.effect.normalizationDirectionAgrees, Boolean(settings.requireCPM),
          settings.minEffect, settings.stableOnly, settings.sequenceRule || "identical", selected.has(candidate.id),
          followup(candidate, settings).text, data.study?.accession || "GSE266899",
          "Exploratory precursor RNA and sequence prioritization; not measured peptide abundance, secretion or anticancer activity."];
      })
    ]);
  }

  const API = { finite, fmt, signed, mean, csv, niceAxis, peptideEvidence, candidateDecision,
    rankCandidates, sampleRows, followup, shortlistCSV };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  if (typeof document === "undefined") return;

  const $ = id => document.getElementById(id);
  const make = (tag, text, className) => {
    const element = document.createElement(tag);
    if (text !== undefined) element.textContent = text;
    if (className) element.className = className;
    return element;
  };
  const put = (id, text) => { if ($(id)) $(id).textContent = text; };
  function svgElement(tag, attrs, text) {
    const element = document.createElementNS("http://www.w3.org/2000/svg", tag);
    Object.entries(attrs || {}).forEach(([key, value]) => element.setAttribute(key, String(value)));
    if (text !== undefined) element.textContent = text;
    return element;
  }
  function label(svg, x, y, text, attrs = {}) {
    svg.append(svgElement("text", { x, y, fill: "#555", "font-size": 12,
      "font-family": "system-ui,sans-serif", ...attrs }, text));
  }
  function table(headers, rows) {
    const element = make("table"), head = make("thead"), heading = make("tr"), body = make("tbody");
    headers.forEach(text => { const cell = make("th", text); cell.scope = "col"; heading.append(cell); });
    head.append(heading);
    rows.forEach(row => {
      const tr = make("tr");
      row.forEach(text => tr.append(make("td", String(text))));
      body.append(tr);
    });
    element.append(head, body); return element;
  }
  function download(filename, content) {
    const url = URL.createObjectURL(new Blob([content], { type: "text/csv;charset=utf-8;" }));
    const anchor = make("a"); anchor.href = url; anchor.download = filename;
    document.body.append(anchor); anchor.click(); anchor.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  }

  let data = null, selectedId = null;
  const comparisonIds = new Set();
  const settings = () => ({ minEffect: Number($("effect-select").value),
    stableOnly: $("stable-only").checked, sequenceRule: $("sequence-rule").value,
    requireCPM: $("cpm-agreement").checked });
  const selectedCandidate = () => data?.candidates.find(candidate => candidate.id === selectedId);

  function drawSamples(candidate) {
    const svg = $("hero-viz"); svg.replaceChildren();
    const width = Math.max(280, Math.min(600, Math.round(svg.getBoundingClientRect?.().width || 600)));
    svg.setAttribute("viewBox", `0 0 ${width} 330`);
    const groups = sampleRows(candidate);
    const values = groups.flatMap(group => group.samples.map(sample => sample.normalizedCount))
      .filter(value => finite(value) && value >= 0);
    if (!values.length) {
      label(svg, width / 2, 150, "No usable sample values.", { "text-anchor": "middle" });
      return;
    }
    const transformed = values.map(value => Math.log2(value + 1));
    const { low, high, ticks } = niceAxis(transformed, true);
    const top = 27, bottom = 253, left = 49, right = width - 12;
    const y = value => bottom - (value - low) / (high - low) * (bottom - top);
    svg.append(svgElement("title", {}, `${candidate.gene}: individual late-time-point pancreatic-lesion RNA values`));
    svg.append(svgElement("desc", {}, "Each labeled dot is one mouse. Values are plotted as log2 of normalized count plus one. Horizontal bars mark arithmetic group means on the normalized-count scale. These animals are not paired across groups."));
    ticks.forEach(value => {
      svg.append(svgElement("line", { x1: left, x2: right, y1: y(value), y2: y(value), stroke: "#e5e5e5" }));
      label(svg, left - 9, y(value) + 4, fmt(value, 0), { "text-anchor": "end" });
    });
    label(svg, 16, 145, "log₂(normalized count + 1)", {
      transform: "rotate(-90 16 145)", "text-anchor": "middle", "font-size": 11
    });
    groups.forEach((group, groupIndex) => {
      const center = left + (groupIndex + .5) * (right - left) / 3;
      const color = group.id === "SEMA" ? "#1f7a8c" : group.id === "HFD" ? "#444" : "#999";
      group.samples.forEach((sample, index) => {
        if (!finite(sample.normalizedCount) || sample.normalizedCount < 0) return;
        const x = center + (index - (group.samples.length - 1) / 2) * 15;
        const yy = y(Math.log2(sample.normalizedCount + 1));
        const description = `${sample.sampleId}: ${fmt(sample.normalizedCount)} normalized counts; raw count ${sample.rawCount}.`;
        const dot = svgElement("g", { tabindex: 0, "aria-label": description });
        dot.append(svgElement("circle", { cx: x, cy: yy, r: 10, fill: "transparent" }));
        dot.append(svgElement("circle", { cx: x, cy: yy, r: 5, fill: color, stroke: "#fff", "stroke-width": 1 }));
        dot.append(svgElement("title", {}, description));
        label(dot, x, yy - 10, sample.sampleId, { "text-anchor": "middle", "font-size": 10, fill: color });
        svg.append(dot);
      });
      if (finite(group.mean)) svg.append(svgElement("line", {
        x1: center - 39, x2: center + 39,
        y1: y(Math.log2(group.mean + 1)), y2: y(Math.log2(group.mean + 1)),
        stroke: color, "stroke-width": 2
      }));
      const groupLabel = width < 470 ? ["Low-fat", "High-fat", "HFD + sema"][groupIndex] : group.label;
      label(svg, center, 280, groupLabel, { "text-anchor": "middle", "font-size": 11, fill: color });
      label(svg, center, 301, `n = ${group.n}`, { "text-anchor": "middle", "font-size": 11 });
    });
    put("sample-caption", "One dot per mouse; horizontal bars show arithmetic group means before log display. Counts use median-ratio library-size normalization, not peptide concentrations. Groups are not paired.");
  }

  function renderSequence() {
    const candidate = selectedCandidate();
    const peptide = candidate?.peptides?.[Number($("peptide-select").value)];
    const container = $("sequence-viz"); container.replaceChildren();
    if (!peptide) {
      container.append(make("p", "No mapped peptide sequence is available for this precursor.", "sequence-empty"));
      put("sequence-caption", "Do not substitute precursor orthology for peptide identity or cleavage evidence.");
      return;
    }
    const evidence = peptideEvidence(peptide);
    const scroll = make("div", undefined, "sequence-scroll");
    const grid = make("div", undefined, "sequence-alignment");
    const rows = [
      { label: "Human", sequence: evidence.human, left: peptide.leftFlankHuman, right: peptide.rightFlankHuman },
      { label: "Mouse", sequence: evidence.mouse, left: peptide.leftFlankMouse, right: peptide.rightFlankMouse }
    ];
    rows.forEach(row => {
      grid.append(make("span", row.label, "sequence-label"));
      const line = make("div", undefined, "sequence-row");
      line.append(make("span", asText(row.left) || "?", "flank"), make("span", "│", "cut"));
      Array.from(row.sequence).forEach((aa, index) => {
        const mismatch = evidence.comparable && evidence.mismatchPositions.includes(index + 1);
        const span = make("span", aa, `residue ${evidence.comparable ? mismatch ? "mismatch" : "conserved" : ""}`);
        span.title = `${row.label} alignment column ${index + 1}: ${aa}${mismatch ? " (differs between species)" : ""}`;
        line.append(span);
      });
      if (!row.sequence) line.append(make("span", "not mapped", "sequence-empty"));
      line.append(make("span", "│", "cut"), make("span", asText(row.right) || "?", "flank"));
      grid.append(line);
    });
    scroll.append(grid); container.append(scroll);
    const positions = evidence.mismatchPositions.length ?
      `Differences at alignment columns ${evidence.mismatchPositions.join(", ")}.` :
      evidence.comparable ? "The mapped peptide sequences are identical." :
      "No residue-by-residue identity is asserted for an absent or unequal-length mapping.";
    container.append(make("p", `Teal residues: identical between species. Underlined gray residues: differences; a dash denotes an alignment gap. Vertical rules separate the peptide from its flanking sequence. ${positions}`, "sequence-legend"));
    const coordinates = [
      finite(peptide.humanStart) ? `human start ${peptide.humanStart}` : "",
      finite(peptide.mouseStart) ? `mouse start ${peptide.mouseStart}` : ""
    ].filter(Boolean).join("; ");
    put("sequence-caption",
      `${evidence.comparable ? `${evidence.identicalPositions}/${evidence.human.length} mapped residues identical. ` : ""}` +
      `Cleavage sites predicted in both species: ${evidence.cleavageConserved ? "yes" : "not confirmed"}. ` +
      `${coordinates ? `${coordinates} (1-based). ` : ""}` +
      "Predicted processing is not demonstrated endogenous peptide release.");
  }

  function renderCandidate() {
    const candidate = selectedCandidate(); if (!candidate) return;
    $("candidate-select").value = selectedId;
    put("candidate-description", `${candidate.gene} · human ${candidate.humanGene || "mapping unavailable"}${candidate.description ? ` · ${candidate.description}` : ""}` +
      (candidate.humanGene === "BRINP2" ? " · BRP precursor reference; not a lesion-data lead." : ""));
    drawSamples(candidate);
    const decision = candidateDecision(candidate, settings());
    const effect = candidate.effect;
    put("readout", `${candidate.gene}: HFD + semaglutide versus HFD has a descriptive log₂ contrast of ${signed(effect.log2FoldChange)}. ` +
      `Leave-one-animal-out range: ${signed(effect.looMin)} to ${signed(effect.looMax)}. ` +
      (decision.eligible ? "Meets the current follow-up rules; peptide presence remains unmeasured." :
        `Not on the current shortlist: ${decision.reasons.join("; ")}.`) +
      (effect.normalizationDirectionAgrees !== true ?
        " Its direction changes with CPM normalization: treat this as normalization-sensitive." : ""));
    const select = $("peptide-select"); select.replaceChildren();
    (candidate.peptides || []).forEach((peptide, index) => {
      const e = peptideEvidence(peptide);
      const option = make("option", `Peptide ${index + 1} · ${e.rawHuman.length} human aa · ` +
        (e.comparable ? `${e.identicalPositions}/${e.human.length} identical residues` : "mapping incomplete"));
      option.value = index; select.append(option);
    });
    select.disabled = !candidate.peptides?.length;
    if (!candidate.peptides?.length) select.append(make("option", "No mapped peptide"));
    else {
      const preferred = candidate.peptides.findIndex(peptide => peptide.cleavageConserved &&
        (settings().sequenceRule === "cleavage" || peptideEvidence(peptide).sequenceConserved));
      select.value = String(Math.max(0, preferred));
    }
    renderSequence();
    const metadata = new Map((data.study.samples || []).map(sample => [sample.id, sample]));
    $("candidate-detail").replaceChildren(table(
      ["Sample", "Group", "Raw count", "Normalized count", "GEO ID"],
      (candidate.samples || []).map(sample => [sample.sampleId, sample.group,
        finite(sample.rawCount) ? sample.rawCount : "not available", fmt(sample.normalizedCount),
        metadata.get(sample.sampleId)?.geoAccession || "see source manifest"])));
    const detail = $("evidence-detail"); detail.replaceChildren();
    detail.append(make("p", `Mouse gene: ${candidate.id}. Human UniProt: ${candidate.uniprot || "not available"}.`, "small-note"));
    detail.append(make("p", `Main-comparison mice with ≥10 raw counts: ${candidate.evidence?.detectedN ?? "not available"}/7. ` +
      `Mapped peptides with predicted cleavage in both species: ${decision.conserved}. ` +
      "Orthology, local sequence identity and cleavage prediction are different evidence layers.", "small-note"));
    detail.append(make("p", `CPM-normalized contrast: ${signed(effect.cpmLog2FoldChange)}. ` +
      `Direction agrees with the median-ratio result: ${effect.normalizationDirectionAgrees ? "yes" : "no"}. ` +
      "Both descriptive ratios use a +1 pseudocount in their respective units.", "small-note"));
    (candidate.limitations || []).forEach(text => detail.append(make("p", text, "small-note")));
    const next = followup(candidate, settings());
    put("next-title", next.title); put("next-step", next.text);
    document.querySelectorAll("#suggestions button").forEach(button =>
      button.setAttribute("aria-pressed", String(button.dataset.candidate === candidate.id)));
  }

  function chooseCandidate(id, scroll = false) {
    selectedId = id; renderCandidate(); renderShortlist();
    if (scroll) $("candidate-select").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function renderShortlist() {
    const ranked = rankCandidates(data.candidates, settings());
    put("metric-three", String(ranked.length));
    const eligible = new Set(ranked.map(result => result.candidate.id));
    for (const id of comparisonIds) if (!eligible.has(id)) comparisonIds.delete(id);
    const container = $("shortlist"); container.replaceChildren();
    if (!ranked.length) {
      container.append(make("p", "No precursor meets all current rules. Relax the contrast floor or inspect direction-sensitive candidates; a zero-length shortlist is a valid result.", "methods-note"));
    } else {
      const list = make("table", undefined, "candidate-table");
      const head = make("thead"), tr = make("tr");
      ["Precursor", "RNA log₂ contrast", "Omit-one range", "Compare"].forEach(text => {
        const cell = make("th", text); cell.scope = "col"; tr.append(cell);
      });
      head.append(tr); list.append(head);
      const body = make("tbody");
      ranked.slice(0, 12).forEach((result, index) => {
        const candidate = result.candidate;
        const row = make("tr", undefined, candidate.id === selectedId ? "selected-row" : "");
        const name = make("td"), button = make("button", `${index + 1}. ${candidate.gene}`, "gene-button");
        button.type = "button"; button.addEventListener("click", () => chooseCandidate(candidate.id, true));
        name.append(button, make("span", `${result.conserved} cleavage-mapped peptide${result.conserved === 1 ? "" : "s"}`, "gene-meta"));
        const effectCell = make("td", signed(candidate.effect.log2FoldChange));
        if (!result.stable) effectCell.append(make("span", "direction-sensitive", "gene-meta"));
        if (!candidate.effect.normalizationDirectionAgrees)
          effectCell.append(make("span", "normalization-sensitive", "gene-meta"));
        const range = make("td", `${signed(candidate.effect.looMin)} to ${signed(candidate.effect.looMax)}`);
        const action = make("td"), compare = make("button", comparisonIds.has(candidate.id) ? "Selected" : "Compare", "plain-button compare-button");
        compare.type = "button"; compare.setAttribute("aria-pressed", String(comparisonIds.has(candidate.id)));
        compare.setAttribute("aria-label", `Compare ${candidate.gene}`);
        compare.disabled = !comparisonIds.has(candidate.id) && comparisonIds.size >= 3;
        compare.addEventListener("click", () => {
          comparisonIds.has(candidate.id) ? comparisonIds.delete(candidate.id) : comparisonIds.add(candidate.id);
          renderShortlist();
        });
        action.append(compare); row.append(name, effectCell, range, action); body.append(row);
      });
      list.append(body);
      const wrap = make("div", undefined, "table-wrap"); wrap.append(list); container.append(wrap);
    }
    put("shortlist-status", `${ranked.length} candidates meet these rules. ` +
      (ranked.length > 12 ? "The first 12 are displayed; the export contains all qualifying candidates. " : "") +
      "Order: stable direction first, then absolute RNA contrast. This is not an activity score.");
    $("export-shortlist").disabled = ranked.length === 0;
    $("clear-comparison").disabled = comparisonIds.size === 0;
    drawComparison();
  }

  function drawComparison() {
    const candidates = [...comparisonIds].map(id => data.candidates.find(candidate => candidate.id === id)).filter(Boolean);
    $("compare-region").hidden = candidates.length === 0;
    const svg = $("comparison-viz"); svg.replaceChildren();
    if (!candidates.length) return;
    const width = Math.max(280, Math.min(600, Math.round(svg.getBoundingClientRect?.().width || 600)));
    svg.setAttribute("viewBox", `0 0 ${width} 210`);
    const values = candidates.flatMap(candidate =>
      [candidate.effect.looMin, candidate.effect.looMax, candidate.effect.log2FoldChange]).filter(finite);
    const { low, high, ticks } = niceAxis(values);
    const left = 76, right = width - 15, bottom = 165;
    const x = value => left + (value - low) / (high - low) * (right - left);
    svg.append(svgElement("title", {}, "Shortlist robustness comparison"));
    ticks.forEach(value => {
      svg.append(svgElement("line", { x1: x(value), x2: x(value), y1: 15, y2: bottom,
        stroke: value === 0 ? "#777" : "#e8e8e8", "stroke-dasharray": value === 0 ? "4 3" : "" }));
      label(svg, x(value), bottom + 18, signed(value), { "text-anchor": "middle", "font-size": 11 });
    });
    candidates.forEach((candidate, index) => {
      const yy = 35 + index * 50;
      label(svg, left - 12, yy + 4, candidate.gene, { "text-anchor": "end" });
      const effect = candidate.effect;
      if (finite(effect.looMin) && finite(effect.looMax)) {
        svg.append(svgElement("line", { x1: x(effect.looMin), x2: x(effect.looMax), y1: yy, y2: yy,
          stroke: "#779498", "stroke-width": 2 }));
        [effect.looMin, effect.looMax].forEach(value => svg.append(svgElement("line", {
          x1: x(value), x2: x(value), y1: yy - 5, y2: yy + 5, stroke: "#779498"
        })));
      }
      if (finite(effect.log2FoldChange)) {
        const point = svgElement("circle", { cx: x(effect.log2FoldChange), cy: yy, r: 5, fill: "#1f7a8c", tabindex: 0,
          "aria-label": `${candidate.gene}: ${signed(effect.log2FoldChange)}; omit-one range ${signed(effect.looMin)} to ${signed(effect.looMax)}.` });
        svg.append(point);
      }
    });
    label(svg, (left + right) / 2, 205, width < 470 ? "log₂ contrast · sema / HFD" :
      "Descriptive log₂ RNA contrast · semaglutide versus HFD", { "text-anchor": "middle", "font-size": 11 });
    if (candidates.length > 1) {
      const largest = [...candidates].sort((a, b) => Math.abs(b.effect.log2FoldChange) - Math.abs(a.effect.log2FoldChange))[0];
      const narrowest = [...candidates].sort((a, b) =>
        (a.effect.looMax - a.effect.looMin) - (b.effect.looMax - b.effect.looMin))[0];
      put("compare-insight", largest.id !== narrowest.id ?
        `${largest.gene} has the larger RNA contrast, while ${narrowest.gene} has the narrower omission range. This separates effect magnitude from dependence on an individual animal; neither establishes peptide activity.` :
        `${largest.gene} has both the largest absolute contrast and narrowest omission range among these selections. Inspect the underlying samples before choosing a detection experiment.`);
    } else put("compare-insight", "Add another precursor to compare effect magnitude with omission sensitivity.");
  }

  function renderMethods() {
    const container = $("methods-detail"); container.replaceChildren();
    const normalization = data.study.normalization || data.normalization;
    if (typeof normalization === "string") container.append(make("p", normalization, "methods-note"));
    else if (normalization) Object.entries(normalization).forEach(([key, value]) => {
      if (typeof value === "string" || typeof value === "number")
        container.append(make("p", `${key}: ${value}`, "methods-note"));
    });
    const rules = data.selectionRules;
    if (typeof rules === "string") container.append(make("p", rules, "methods-note"));
    else if (rules) Object.entries(rules).forEach(([key, value]) => {
      if (typeof value === "string" || typeof value === "number")
        container.append(make("p", `${key}: ${value}`, "methods-note"));
    });
    if (data.summary) {
      const summary = data.summary;
      container.append(make("p", `Coverage: ${summary.publishedPrecursors} precursors in the published library; ` +
        `${summary.shortPrecursors} have published peptide sequences of 5–25 residues; ` +
        `${summary.sequenceVerifiedPrecursors} passed the unique orthology and pinned-sequence checks used here. ` +
        "This explorer is the sequence-audited subset, not a complete search of all secreted proteins.", "methods-note"));
    }
    if (data.study.depositionDiscrepancy)
      container.append(make("p", data.study.depositionDiscrepancy, "methods-note"));
  }

  async function load() {
    try {
      const config = window.UNIVERSITY_LAB_DEMO_DATA || {};
      const response = await fetch(`${config.candidateFile || "peptide-cancer-data.json"}?v=${encodeURIComponent(config.version || "peptide-cancer-1")}`);
      if (!response.ok) throw new Error(`Candidate data request returned ${response.status}.`);
      data = await response.json();
      if (!Array.isArray(data.candidates) || !data.candidates.length || !Array.isArray(data.study?.samples))
        throw new Error("Verified candidate data are incomplete.");
      const select = $("candidate-select"); select.replaceChildren();
      data.candidates.forEach(candidate => {
        const option = make("option", `${candidate.gene} · human ${candidate.humanGene || "unmapped"}` +
          (candidate.effect.normalizationDirectionAgrees !== true ? " · normalization-sensitive" : ""));
        option.value = candidate.id; select.append(option);
      });
      select.disabled = false;
      const ranked = rankCandidates(data.candidates, settings());
      const initial = ranked.find(result => result.candidate.gene === "Ntn1");
      selectedId = initial?.candidate.id || ranked[0]?.candidate.id || data.candidates[0].id;
      const suggestions = $("suggestions"); suggestions.replaceChildren(make("span", "Inspect:"));
      const examples = [...ranked.slice(0, 3).map(result => result.candidate)];
      const reference = data.candidates.find(candidate => candidate.humanGene === "BRINP2");
      const fragile = data.candidates.find(candidate => candidate.evidence?.measurable && !candidate.effect?.stableDirection);
      const additional = reference || fragile;
      if (additional && !examples.some(candidate => candidate.id === additional.id)) examples.push(additional);
      const normalizationSensitive = ranked.find(result => !result.candidate.effect.normalizationDirectionAgrees)?.candidate;
      if (normalizationSensitive && !examples.some(candidate => candidate.id === normalizationSensitive.id))
        examples.push(normalizationSensitive);
      examples.forEach(candidate => {
        const button = make("button", `${candidate.gene}${candidate.humanGene === "BRINP2" ? " (reference)" : ""}`);
        button.type = "button"; button.dataset.candidate = candidate.id;
        button.addEventListener("click", () => chooseCandidate(candidate.id));
        suggestions.append(button);
      });
      select.addEventListener("change", () => chooseCandidate(select.value));
      $("peptide-select").addEventListener("change", renderSequence);
      ["effect-select", "stable-only", "sequence-rule", "cpm-agreement"].forEach(id => $(id).addEventListener("change", () => {
        renderCandidate(); renderShortlist();
      }));
      $("clear-comparison").addEventListener("click", () => { comparisonIds.clear(); renderShortlist(); });
      $("export-shortlist").addEventListener("click", () => download("svensson-precursor-followup-shortlist.csv",
        shortlistCSV(data, settings(), [...comparisonIds])));
      put("metric-one", String(data.study.samples.length));
      put("metric-two", String(data.candidates.length));
      put("selection-rule", "Always require ≥10 raw counts in at least 3 of the 7 HFD/semaglutide mice and mapped cleavage evidence. Strict mode also requires an identical human–mouse peptide. These rules prioritize detection experiments, not activity.");
      const ntn1 = data.candidates.find(candidate => candidate.gene === "Ntn1");
      const coch = data.candidates.find(candidate => candidate.gene === "Coch");
      if (ntn1 && coch && reference) put("findings-note",
        `Ntn1 is the initial inspection example: its RNA contrast stays between ${signed(ntn1.effect.looMin)} and ${signed(ntn1.effect.looMax)} when one animal is omitted, ` +
        `compared with ${signed(coch.effect.looMin)} to ${signed(coch.effect.looMax)} for Coch. ` +
        `BRINP2 has ${reference.evidence.detectedN}/7 main-comparison samples reaching 10 raw counts. ` +
        "These differences favor an evidence-based detection shortlist over simply picking the biggest contrast or the best-known peptide.");
      ranked.slice(0, 2).forEach(result => comparisonIds.add(result.candidate.id));
      renderCandidate(); renderShortlist(); renderMethods();
      window.addEventListener?.("resize", () => {
        if (selectedCandidate()) drawSamples(selectedCandidate());
        drawComparison();
      });
      put("load-status", "Verified sample data and predicted sequence evidence loaded. All selections run locally.");
    } catch (error) {
      put("load-status", `The checked data could not be loaded: ${error.message} No analysis result is shown.`);
      console.error(error);
    }
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", load, { once: true });
  else load();
})();
