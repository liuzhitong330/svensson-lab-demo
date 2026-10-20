(function () {
  "use strict";
  const finite = value => typeof value === "number" && Number.isFinite(value);
  const fmt = (value, digits = 1) => finite(value) ? value.toFixed(digits) : "not available";
  const percent = value => finite(value) ? value * 100 : null;
  const mean = values => values.length && values.every(finite) ?
    values.reduce((sum, value) => sum + value, 0) / values.length : null;
  const quoteCSV = value => `"${String(value ?? "").replace(/"/g, '""')}"`;
  const csv = rows => rows.map(row => row.map(quoteCSV).join(",")).join("\r\n") + "\r\n";

  function niceAxis(values, nonnegative = false) {
    const numbers = values.filter(finite);
    const minimum = numbers.length ? Math.min(0, ...numbers) : 0;
    const maximum = numbers.length ? Math.max(0, ...numbers) : 1;
    const rawStep = (maximum - minimum || 1) / 4;
    const power = 10 ** Math.floor(Math.log10(rawStep));
    const fraction = rawStep / power;
    const factor = [1, 2, 2.5, 5, 10].find(value => fraction <= value);
    const step = factor * power;
    const low = nonnegative ? 0 : Math.floor(minimum / step) * step;
    let high = Math.ceil(maximum / step) * step;
    if (high <= low) high = low + step * 4;
    const count = Math.round((high - low) / step);
    const ticks = Array.from({ length: count + 1 }, (_, index) =>
      Number((low + index * step).toPrecision(12)));
    return { low, high, step, ticks };
  }

  function peptideResult(condition, thresholdPercent, timeHours) {
    const feeding = condition?.feeding?.find(row => row.time_h === Number(timeHours));
    const residuals = [
      percent(condition?.ns1?.retained_activity),
      percent(condition?.neuro2a?.retained_activity),
      percent(feeding?.retained_inhibition)
    ];
    const complete = Boolean(condition?.ns1 && condition?.neuro2a && feeding &&
      condition.ns1.n >= 3 && condition.neuro2a.n >= 3 && feeding.n >= 3 &&
      !condition.ns1.missing_count && !condition.neuro2a.missing_count && !feeding.missing_count);
    const available = residuals.every(finite) && finite(Number(thresholdPercent)) && Number(thresholdPercent) >= 0;
    const eligible = complete && available && residuals.every(value => Math.abs(value) <= Number(thresholdPercent));
    const status = !available ? "unavailable" : !complete ? "repeat missing measurement" :
      eligible ? "provisional control candidate" : "does not meet selected ceiling";
    return { condition, feeding, residuals, complete, available, eligible, status };
  }

  function peptideCSV(data, ceiling, timeHours) {
    return csv([
      ["Peptide", "Sequence", "Ceiling percent", "Feeding time hours", "NS-1 retained activity percent",
        "Neuro2a retained activity percent", "Feeding retained inhibition percent", "Status", "Notes"],
      ...data.conditions.slice(2).map(condition => {
        const result = peptideResult(condition, ceiling, timeHours);
        return [condition.id, condition.sequence, ceiling, timeHours, ...result.residuals,
          result.status, "Descriptive screen, not an equivalence test. Feeding n=3 experiment means; cell n=3 except NS-1 R3A n=2."];
      })
    ]);
  }

  function peptideIdentity(data, condition) {
    const source = data.repo_identity.results.find(row => row.gene === "BRINP2" && row.minimum_spacing_aa === 5);
    const mutationPosition = condition?.mutation_position;
    return {
      source,
      mutationPosition,
      precursorPosition: source && finite(mutationPosition) ?
        source.brp_exact_match_start_1based + mutationPosition - 1 : null
    };
  }

  function cancerResult(data, endpointKey, icbRequested, contrastKey) {
    const endpoint = data.endpoints[endpointKey];
    if (!endpoint) return null;
    const icb = endpointKey !== "mass_change_week22" && Boolean(icbRequested);
    const groups = ["nc", "wd", "sema"].map(key =>
      endpoint.groups.find(group => group.key === key + (icb ? "_icb" : "")));
    if (groups.some(group => !group || !finite(group.mean))) return null;
    const target = groups[contrastKey === "nc" ? 0 : 2], reference = groups[1];
    const difference = target.mean - reference.mean;
    const differenceUnit = endpoint.unit === "%" ? "percentage points" : endpoint.unit;
    return { endpointKey, endpoint, icb, groups, target, reference, difference, differenceUnit };
  }

  function cancerCSV(data, endpointKey, icb, contrast) {
    const result = cancerResult(data, endpointKey, icb, contrast);
    if (!result) return csv([["Error"], ["Endpoint or group unavailable"]]);
    const rows = [["Endpoint", "Unit", "Group", "ICB", "Source record index", "Value", "Source cell",
      "Biological unit", "Selected contrast", "Difference of group means", "Difference unit"]];
    result.groups.forEach(group => group.values.forEach((value, index) => {
      rows.push([result.endpoint.title, result.endpoint.unit, group.label, result.icb ? "yes" : "no",
        index + 1, value, group.sourceCells?.[index] || "", result.endpoint.biologicalUnit,
        `${result.target.label} minus ${result.reference.label}`, result.difference, result.differenceUnit]);
    }));
    return csv(rows);
  }

  function nextControl(result) {
    if (!result) return { title: "Check the missing measurement", text: "No comparison can be made until the source group is available." };
    const targetIsSema = result.target.key.startsWith("sema");
    const context = result.icb ? "with the same checkpoint-blockade context" : "without checkpoint blockade";
    if (result.endpointKey === "mass_change_week22") {
      return {
        title: targetIsSema ? "Separate body-mass change from tumor response" : "Separate diet composition from body mass",
        text: targetIsSema ?
          "Proposed next control: collect food intake, body composition and tumor response under diet-matched vehicle and intake-matched controls. If the tumor response follows matched intake or body composition, investigate host-metabolic explanations; if it differs despite comparable host measures, test receptor-dependent cellular activity. Neither pattern alone proves a mechanism." :
          "Proposed next control: compare nutrition-matched cohorts and prespecified body-composition measurements before attributing a tumor difference to body mass. NC is a normal-chow cohort, not a diet-switch or weight-matched control."
      };
    }
    if (result.endpointKey === "tumor_day16") {
      return {
        title: targetIsSema ? "Distinguish tumor-cell activity from host effects" : "Test the contribution of the host diet context",
        text: targetIsSema ?
          `Proposed next control: pair a receptor-validated tumor-cell assay and viability readout with host-metabolic and immune measurements ${context}. A receptor-dependent cell effect would prioritize cell-signaling follow-up. An absent cell effect alongside a host change would prioritize metabolic or immune hypotheses, but would not prove them. Preserve experiment and animal IDs.` :
          `Proposed next control: retain matched diet, immune and tumor-growth records ${context}, and add an intake-matched comparator. A normal-chow versus Western-diet contrast changes several host factors; it is not an isolated drug comparison.`
      };
    }
    return {
      title: targetIsSema ? "Check immune function as its own outcome" : "Test whether the immune phenotype follows the diet context",
      text: targetIsSema ?
        `Proposed next control: validate tumor-antigen-specific cell abundance and function ${context}, alongside tumor size in newly linked animals. A functional change would prioritize immune follow-up; a tumor change without it would prioritize other host or tumor-cell mechanisms. These are hypotheses. Match flow gates and sampling time; the existing panels cannot be linked per animal.` :
        `Proposed next control: retain comparable antigen-specific flow gates and sampling time ${context}; examine diet and immune function in an independently replicated cohort. Report denominators and batch IDs, especially when a group contains only three mouse summaries.`
    };
  }

  const API = { finite, fmt, mean, niceAxis, peptideResult, peptideCSV, peptideIdentity, cancerResult, cancerCSV, nextControl, csv };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  if (typeof document === "undefined") return;

  const $ = id => document.getElementById(id);
  const make = (tag, text, className) => {
    const element = document.createElement(tag);
    if (text !== undefined) element.textContent = text;
    if (className) element.className = className;
    return element;
  };
  function put(id, text) { if ($(id)) $(id).textContent = text; }
  function svgElement(tag, attrs, text) {
    const element = document.createElementNS("http://www.w3.org/2000/svg", tag);
    Object.entries(attrs || {}).forEach(([key, value]) => element.setAttribute(key, String(value)));
    if (text !== undefined) element.textContent = text;
    return element;
  }
  function table(headers, rows) {
    const element = make("table"), head = make("thead"), body = make("tbody"), heading = make("tr");
    headers.forEach(label => {
      const th = make("th", label); th.scope = "col"; heading.append(th);
    });
    head.append(heading);
    rows.forEach(row => {
      const tr = make("tr"); row.forEach(value => tr.append(make("td", String(value)))); body.append(tr);
    });
    element.append(head, body); return element;
  }
  function download(filename, content) {
    const url = URL.createObjectURL(new Blob([content], { type: "text/csv;charset=utf-8;" }));
    const anchor = make("a"); anchor.href = url; anchor.download = filename;
    document.body.append(anchor); anchor.click(); anchor.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  }
  function axisLabel(svg, x, y, text, extra = {}) {
    svg.append(svgElement("text", { x, y, fill: "#555", "font-size": 11,
      "font-family": "system-ui, sans-serif", ...extra }, text));
  }

  const config = window.UNIVERSITY_LAB_DEMO_DATA || {};
  let peptideData, cancerData;

  function drawPeptide(result, ceiling) {
    const svg = $("hero-viz"); if (!svg) return;
    svg.replaceChildren();
    svg.append(svgElement("title", {}, "Retained Fos activity across two cell lines"));
    svg.append(svgElement("desc", {}, "Each point is an alanine-substituted BRP peptide. Select a point by click or keyboard. The shaded box shows the selected maximum residual activity in both cell lines; feeding is checked separately."));
    const width = 600, height = 360, left = 58, right = 24, top = 26, bottom = 58;
    const xMax = 350, yMax = 250;
    const x = value => left + value / xMax * (width-left-right);
    const y = value => height-bottom-value / yMax * (height-top-bottom);
    const axis = { stroke: "#bbb", "stroke-width": 1 };
    svg.append(svgElement("rect", { x: x(0), y: y(ceiling), width: x(ceiling)-x(0),
      height: y(0)-y(ceiling), fill: "#e4efed", stroke: "#8ea7a1", "stroke-dasharray": "3 3" }));
    [0,50,100,150,200,250,300,350].forEach(value => {
      svg.append(svgElement("line", { x1:x(value), x2:x(value), y1:y(0), y2:y(250), stroke:"#eee" }));
      axisLabel(svg, x(value), y(0)+18, value, { "text-anchor":"middle" });
    });
    [0,50,100,150,200,250].forEach(value => {
      svg.append(svgElement("line", { x1:x(0), x2:x(350), y1:y(value), y2:y(value), stroke:"#eee" }));
      axisLabel(svg, left-9, y(value)+4, value, { "text-anchor":"end" });
    });
    svg.append(svgElement("line", { x1:x(0), x2:x(350), y1:y(0), y2:y(0), ...axis }));
    svg.append(svgElement("line", { x1:x(0), x2:x(0), y1:y(0), y2:y(250), ...axis }));
    axisLabel(svg, 322, 348, "NS-1 retained Fos activity (%)", { "text-anchor":"middle" });
    axisLabel(svg, 15, 164, "Neuro2a retained Fos activity (%)", { transform:"rotate(-90 15 164)", "text-anchor":"middle" });
    const labelOffset = { L8A:[5,-10], L5A:[6,-8], C12A:[5,-10], T1A:[6,-6], H2A:[6,-6],
      R3A:[6,-7], I4A:[6,13], R6A:[-6,-9], R7A:[6,-8], F9A:[-7,15], N10A:[6,-6], L11A:[6,-7] };
    peptideData.conditions.slice(2).forEach(condition => {
      const xx = percent(condition.ns1.retained_activity), yy = percent(condition.neuro2a.retained_activity);
      if (!finite(xx) || !finite(yy)) return;
      const chosen = condition.id === result.condition.id;
      const group = svgElement("g", { role:"button", tabindex:0,
        "aria-label":`${condition.id}: NS-1 ${fmt(xx)}%, Neuro2a ${fmt(yy)}%. Select peptide.`,
        "aria-pressed":String(chosen), "data-peptide":condition.id, style:"cursor:pointer" });
      group.append(svgElement("circle", { cx:x(xx), cy:y(yy), r:14, fill:"transparent" }));
      group.append(svgElement("circle", { cx:x(xx), cy:y(yy), r:chosen ? 6.5 : 4.5,
        fill:chosen ? "#1f7a8c" : "#8a918e", stroke:"#fff", "stroke-width":1.5 }));
      const offset = labelOffset[condition.id] || [6,-6];
      axisLabel(group, x(xx)+offset[0], y(yy)+offset[1], condition.id,
        { "text-anchor":offset[0]<0 ? "end":"start", fill:chosen ? "#1f7a8c":"#555",
          "font-weight":chosen ? "600":"400" });
      const activate = () => { $("peptide-select").value = condition.id; renderPeptide(); };
      group.addEventListener("click", activate);
      group.addEventListener("keydown", event => {
        if (event.key === "Enter" || event.key === " ") { event.preventDefault(); activate(); }
      });
      svg.append(group);
    });
  }

  function renderPeptide() {
    const condition = peptideData.conditions.find(row => row.id === $("peptide-select").value);
    if (!condition) return;
    const ceiling = Number($("residual-select").value), time = Number($("time-select").value);
    const result = peptideResult(condition, ceiling, time);
    drawPeptide(result, ceiling);
    put("readout", `${condition.id}: ${result.status}. At the ${fmt(ceiling,0)}% ceiling and ${time} h, ` +
      `retained responses are ${fmt(result.residuals[0])}% in NS-1, ${fmt(result.residuals[1])}% in Neuro2a, ` +
      `and ${fmt(result.residuals[2])}% for feeding inhibition. This is a mean-based triage rule, not proof of inactivity or equivalence.`);
    const details = $("peptide-detail");
    if (details) {
      const rows = [
        ["NS-1 Fos", condition.ns1, "relative expression"],
        ["Neuro2a Fos", condition.neuro2a, "relative expression"],
        [`Feeding at ${time} h`, result.feeding, "g"]
      ].map(([label, assay, unit]) => [label,
        assay ? assay.replicates.map(v=>fmt(v,3)).join(", ") : "not available",
        assay ? `${fmt(assay.mean,3)} ± ${fmt(assay.sd,3)} ${unit}` : "not available",
        assay ? `${assay.n}${assay.missing_count ? `; ${assay.missing_count} missing` : ""}` : "0"]);
      details.replaceChildren(table(["Measurement", "Reported values", "Mean ± sample SD", "n"], rows));
      details.append(make("p", "Fos: biological replicates. Feeding: three independent experiment means, each from three mice. No observations are paired across assays. Normalized vehicle values of 1 do not imply zero raw assay noise.", "small-note"));
    }
    const accepted = peptideData.conditions.slice(2).map(c=>peptideResult(c,ceiling,time)).filter(r=>r.eligible);
    put("metric-three", String(accepted.length));
    const shortlist = $("shortlist");
    if (shortlist) {
      shortlist.replaceChildren(make("p", accepted.length ?
        `Provisional shortlist: ${accepted.map(r=>r.condition.id).join(", ")} (${accepted.length} of 12 mutants).` :
        "No mutant meets all three selected response ceilings with complete observations."));
      shortlist.append(make("p", result.eligible ?
        `Next for ${condition.id}: repeat independent dose-response and viability controls before using it as an inactive comparator. Test the selected feeding window explicitly.` :
        !result.complete ? `Next for ${condition.id}: repeat its missing NS-1 measurement before making a control-selection decision.` :
        `Next for ${condition.id}: retain the discordant assay in the validation panel. Low activity in one cell line does not establish an inactive control in another cell type or in vivo.`));
    }
    const identity = peptideIdentity(peptideData, condition);
    put("repo-detail", `Selected sequence: ${condition.sequence}. ${condition.id} changes peptide residue ` +
      `${identity.mutationPosition}, corresponding to BRINP2 (Q9C0B6) residue ${identity.precursorPosition}. ` +
      `The lab’s pinned repository maps BRP to residues ${identity.source.brp_exact_match_start_1based}–` +
      `${identity.source.brp_exact_match_end_1based}, between KK384–385 and KR398–399. ` +
      `Its upstream motif component detects ${identity.source.cleavage_count} sites in BRINP2; ` +
      "spacing 3 versus 5 aa leaves this mapping unchanged. This identity check does not establish physiological processing or receptor specificity.");
    document.querySelectorAll("#suggestions button").forEach(button =>
      button.setAttribute("aria-pressed", String(button.dataset.peptide === condition.id)));
  }

  function drawCancer(result) {
    const svg = $("cancer-viz"); if (!svg) return;
    svg.replaceChildren();
    svg.append(svgElement("title", {}, result.endpoint.title));
    svg.append(svgElement("desc", {}, "Dots show individual mouse observations; horizontal lines show group means, not confidence intervals. Mouse identifiers are local to each panel. Tumor and immune panels pool two experiments without available batch identifiers."));
    const all = result.groups.flatMap(group=>group.values).filter(finite);
    const left = 68, right = 18, top = 28, bottom = 65, height = 340;
    const {low, high, ticks} = niceAxis(all, result.endpointKey !== "mass_change_week22");
    const y = value => height-bottom-(value-low)/(high-low)*(height-top-bottom);
    for(const value of ticks) {
      svg.append(svgElement("line", { x1:left,x2:600-right,y1:y(value),y2:y(value),stroke:"#e5e5e5" }));
      axisLabel(svg,left-8,y(value)+4,fmt(value,result.endpointKey==="tumor_day16"?0:1),{"text-anchor":"end"});
    }
    if (low<0 && high>0) svg.append(svgElement("line",{x1:left,x2:582,y1:y(0),y2:y(0),stroke:"#999","stroke-dasharray":"3 3"}));
    result.groups.forEach((group,index)=>{
      const center=150+index*170, chosen=group.key===result.target.key, color=chosen?"#1f7a8c":"#8b928f";
      group.values.forEach((value,i)=>{
        const jitter=((i*7)%13-6)*3.1;
        const point=svgElement("circle",{cx:center+jitter,cy:y(value),r:3.7,fill:color,opacity:.64});
        point.append(svgElement("title",{},`${group.label}, source record ${i+1}: ${fmt(value,3)} ${result.endpoint.unit}`));
        svg.append(point);
      });
      svg.append(svgElement("line",{x1:center-34,x2:center+34,y1:y(group.mean),y2:y(group.mean),stroke:chosen?"#164b43":"#333","stroke-width":2.4}));
      axisLabel(svg,center,294,group.label,{"text-anchor":"middle",fill:chosen?"#1f7a8c":"#444"});
      axisLabel(svg,center,312,`n = ${group.n}`,{"text-anchor":"middle"});
    });
    axisLabel(svg,15,160,result.endpoint.unit,{transform:"rotate(-90 15 160)","text-anchor":"middle"});
  }

  function renderCancer() {
    const endpoint = $("endpoint-select").value;
    const icbSelect = $("icb-select");
    icbSelect.disabled = endpoint === "mass_change_week22";
    if (icbSelect.disabled) icbSelect.value = "off";
    const result = cancerResult(cancerData, endpoint, icbSelect.value === "on", $("contrast-select").value);
    if (!result) { put("cancer-readout","Source group unavailable; no comparison computed."); return; }
    drawCancer(result);
    put("cancer-caption", `${result.endpoint.panel}. ${result.endpoint.biologicalUnit}. ` +
      (result.icb ? "ICB = anti-PD-1 plus anti-CTLA-4. " : "") +
      "Dots are source records; bars are means. No cross-panel pairing or inferential test.");
    put("cancer-readout", `${result.target.label} minus ${result.reference.label}: ` +
      `${result.difference>0?"+":""}${fmt(result.difference,2)} ${result.differenceUnit}. ` +
      `Means are ${fmt(result.target.mean,2)} versus ${fmt(result.reference.mean,2)} ${result.endpoint.unit}. ` +
      "This describes the available groups; it is not an adjusted or causal estimate.");
    const detail=$("cancer-detail");
    if(detail) {
      detail.replaceChildren(table(["Group","Source n",`Mean ± SD (${result.endpoint.unit})`],
        result.groups.map(g=>[g.label,g.n,`${fmt(g.mean,2)} ± ${fmt(g.sd,2)}`])));
      const details=make("details"), summary=make("summary","Inspect individual source values");
      details.append(summary);
      result.groups.forEach(g=>details.append(make("p",`${g.label}: ${g.values.map(v=>fmt(v,2)).join(", ")} ${result.endpoint.unit}.`)));
      detail.append(details);
    }
    const next=nextControl(result); put("next-title",next.title); put("next-control",next.text);
  }

  async function load() {
    try {
      const get = async path => {
        const response=await fetch(`${path}?v=${encodeURIComponent(config.version||"1")}`);
        if(!response.ok) throw new Error(`Data request failed (${response.status}).`);
        return response.json();
      };
      [peptideData,cancerData] = await Promise.all([
        get(config.peptideFile||"brp-data.json"), get(config.cancerFile||"cancer-data.json")
      ]);
      if (!Array.isArray(peptideData.conditions) || peptideData.conditions.length!==14) throw new Error("Unexpected peptide dataset.");
      const select=$("peptide-select"); select.replaceChildren();
      peptideData.conditions.slice(2).forEach(condition=>{
        const option=make("option",`${condition.id} · ${condition.sequence}`); option.value=condition.id; select.append(option);
      });
      select.value=config.defaultPeptide||"L8A";
      const suggestions=$("suggestions");
      if(suggestions) {
        suggestions.replaceChildren();
        ["L8A","L5A","C12A","R3A"].forEach(id=>{
          const button=make("button",id); button.type="button"; button.dataset.peptide=id;
          button.addEventListener("click",()=>{select.value=id;renderPeptide();}); suggestions.append(button);
        });
      }
      ["peptide-select","residual-select","time-select"].forEach(id=>$(id).addEventListener("change",renderPeptide));
      ["endpoint-select","icb-select","contrast-select"].forEach(id=>$(id).addEventListener("change",renderCancer));
      $("export-peptides")?.addEventListener("click",()=>download(
        `brp-control-shortlist-${$("residual-select").value}pct-${$("time-select").value}h.csv`,
        peptideCSV(peptideData,Number($("residual-select").value),Number($("time-select").value))));
      $("export-cancer")?.addEventListener("click",()=>download(
        `cancer-${$("endpoint-select").value}-${$("icb-select").value}.csv`,
        cancerCSV(cancerData,$("endpoint-select").value,$("icb-select").value==="on",$("contrast-select").value)));
      renderPeptide(); renderCancer();
      put("load-status","Published source data loaded. Controls recompute the descriptive comparisons locally.");
    } catch(error) {
      put("load-status",`Could not load the analysis: ${error.message} Please reload the public page. No result is available.`);
      console.error(error);
    }
  }
  if(document.readyState==="loading") document.addEventListener("DOMContentLoaded",load,{once:true});
  else load();
})();
