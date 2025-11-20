// ---------- Utilities ----------
const cleanSeq = (s) =>
  (s || "").toUpperCase().split("").filter(c => /[A-Z]/.test(c)).join("");

// Chou–Fasman parameters 
const P_ALPHA = {
  E:1.53, A:1.45, L:1.34, H:1.24, M:1.20, Q:1.17,
  W:1.14, V:1.14, F:1.12, K:1.07, I:1.00, D:0.98,
  T:0.82, S:0.79, R:0.79, C:0.77, N:0.73, Y:0.61,
  P:0.59, G:0.53
};
const P_BETA = {
  M:1.67, V:1.65, I:1.60, C:1.30, Y:1.29, F:1.28,
  Q:1.23, L:1.22, T:1.20, W:1.19, A:0.97, R:0.90,
  G:0.81, D:0.80, K:0.74, S:0.72, H:0.71, N:0.65,
  P:0.62, E:0.26
};

const pa = r => P_ALPHA[r] ?? 0.0;
const pb = r => P_BETA[r]  ?? 0.0;

const windowHasMinOver1 = (vals, minCount) =>
  vals.filter(v => v > 1.0).length >= minCount;

function extendSegment(start, end, prop, thresholdSum, n){
  let changed = true;
  while(changed){
    changed = false;
    // Left extension
    if(start > 0 && (end - start + 1) >= 3){
      const w = [prop[start-1], ...prop.slice(start, start+3)];
      if(w.reduce((a,b)=>a+b,0) >= thresholdSum){
        start--; changed = true;
      }
    }
    // Right extension
    if(end < n-1 && (end - start + 1) >= 3){
      const w = [...prop.slice(end-2, end+1), prop[end+1]];
      if(w.reduce((a,b)=>a+b,0) >= thresholdSum){
        end++; changed = true;
      }
    }
  }
  return [start, end];
}

function mergeSegments(segs){
  if(!segs.length) return [];
  segs.sort((a,b)=>a[0]-b[0]);
  const out=[segs[0]];
  for(let i=1;i<segs.length;i++){
    const [s,e]=segs[i]; const last=out[out.length-1];
    if(s<=last[1]+1){
      last[1]=Math.max(last[1],e);
    } else {
      out.push([s,e]);
    }
  }
  return out;
}

function formatSegs(seq, segs){
  return segs.map(([s,e])=>({
    start:s+1, end:e+1, length:e-s+1, seq:seq.slice(s,e+1)
  }));
}

// ---------- Core Chou–Fasman ----------
function predictChouFasman(seq){
  const SEQ = cleanSeq(seq);
  const n = SEQ.length;
  const Palpha = [...SEQ].map(pa);
  const Pbeta  = [...SEQ].map(pb);

  // Helix nucleation (6-mer, >=4 residues with Pa>1)
  let helixSegs = [];
  const wH = 6;
  for(let i=0;i<=n-wH;i++){
    const vals = Palpha.slice(i,i+wH);
    if(windowHasMinOver1(vals,4)){
      let s=i, e=i+wH-1;
      [s,e] = extendSegment(s,e,Palpha,4.0,n);
      helixSegs.push([s,e]);
    }
  }
  helixSegs = mergeSegments(helixSegs);

  // Strand nucleation (5-mer, >=3 residues with Pb>1)
  let strandSegs = [];
  const wS = 5;
  for(let i=0;i<=n-wS;i++){
    const vals = Pbeta.slice(i,i+wS);
    if(windowHasMinOver1(vals,3)){
      let s=i, e=i+wS-1;
      [s,e] = extendSegment(s,e,Pbeta,4.0000001,n);
      strandSegs.push([s,e]);
    }
  }
  strandSegs = mergeSegments(strandSegs);

  // Label layers
  const Lh = Array(n).fill('-');
  const Ls = Array(n).fill('-');
  helixSegs.forEach(([s,e])=>{ for(let k=s;k<=e;k++) Lh[k]='H'; });
  strandSegs.forEach(([s,e])=>{ for(let k=s;k<=e;k++) Ls[k]='S'; });

  // Conflicts
  const conflictMask = Lh.map((c,i)=> c==='H' && Ls[i]==='S');
  const conflicts = [];
  let i=0;
  while(i<n){
    if(conflictMask[i]){
      let j=i;
      while(j+1<n && conflictMask[j+1]) j++;
      conflicts.push([i,j]);
      i=j+1;
    } else i++;
  }

  // Final labels (with temporary X)
  const final = Array(n).fill('-');
  for(let k=0;k<n;k++){
    if(Lh[k]==='H') final[k]='H';
    if(Ls[k]==='S') final[k] = final[k]==='-' ? 'S' : 'X';
  }

  // Resolve conflicts by avg Pa vs avg Pb
  const resolved=[];
  conflicts.forEach(([s,e])=>{
    const len=e-s+1;
    const avgPa = Palpha.slice(s,e+1).reduce((a,b)=>a+b,0)/len;
    const avgPb = Pbeta .slice(s,e+1).reduce((a,b)=>a+b,0)/len;
    const chosen = avgPb>avgPa ? 'S' : 'H';
    for(let k=s;k<=e;k++) final[k]=chosen;
    resolved.push({
      start:s+1,
      end:e+1,
      length:len,
      avg_pa:+avgPa.toFixed(3),
      avg_pb:+avgPb.toFixed(3),
      chosen,
      seq:SEQ.slice(s,e+1),
    });
  });

  // Clean any leftover X
  for(let k=0;k<n;k++){
    if(final[k]==='X') final[k]='H';
  }

  const finalStr = final.join('');
  const blocks = [];
  for(let i=0;i<n;i+=60){
    blocks.push({
      range:`${i+1}-${Math.min(i+60,n)}`,
      seq: SEQ.slice(i,i+60),
      ann: finalStr.slice(i,i+60),
    });
  }

  return {
    clean_seq: SEQ,
    helix_segments: formatSegs(SEQ, helixSegs),
    strand_segments: formatSegs(SEQ, strandSegs),
    conflicts: resolved,
    final_string: finalStr,
    blocks,
    length:n,
  };
}

// ---------- GOR placeholder ----------
function predictGOR(seq){
  const SEQ = cleanSeq(seq);
  const dash = "-".repeat(SEQ.length);
  return {
    clean_seq: SEQ,
    helix_segments: [],
    strand_segments: [],
    conflicts: [],
    final_string: dash,
    blocks: [{
      range:`1-${SEQ.length}`,
      seq: SEQ,
      ann: dash,
    }],
    length: SEQ.length,
  };
}

// ---------- UI ----------
const els = {
  seq: document.getElementById('sequence'),
  method: document.getElementById('method'),
  run: document.getElementById('run'),
  results: document.getElementById('results'),
  summary: document.getElementById('summary'),
  blocks: document.getElementById('blocks'),
  helixList: document.getElementById('helix-list'),
  strandList: document.getElementById('strand-list'),
  conflicts: document.getElementById('conflicts'),
  copyBtn: document.getElementById('copy'),
  dlBtn: document.getElementById('download'),
};

function renderList(container, items){
  container.innerHTML = '';
  if(!items || !items.length){
    container.innerHTML = '<li>None</li>';
    return;
  }
  container.innerHTML = items.map(s =>
    `<li>${s.start}–${s.end} (len=${s.length}) — <span class="dim">${s.seq}</span></li>`
  ).join('');
}

function colorizeAnn(text){
  let html = '';
  for(const ch of text){
    if(ch === 'H') html += `<span class="ann-h">H</span>`;
    else if(ch === 'S') html += `<span class="ann-s">S</span>`;
    else html += `<span class="ann-gap">${ch}</span>`;
  }
  return html;
}

function renderBlocks(container, blocks){
  container.innerHTML = blocks.map(b => (
    `<div class="block">
       <div class="line meta">${b.range}</div>
       <div class="line seq">${b.seq}</div>
       <div class="line ann">${colorizeAnn(b.ann)}</div>
     </div>`
  )).join('');
}

function buildTextReport(seq, method, res){
  const lines = [];
  lines.push(`Method: ${method === 'chou_fasman' ? 'Chou–Fasman' : 'GOR (placeholder)'}`);
  lines.push(`Sequence length: ${seq.length}`);
  lines.push(``);
  lines.push(`FINAL ASSIGNMENT:`);
  lines.push(res.final_string);
  lines.push(``);
  lines.push(`HELICAL REGIONS (H):`);
  if(res.helix_segments.length){
    res.helix_segments.forEach(s =>
      lines.push(`  ${s.start}-${s.end} len=${s.length} ${s.seq}`)
    );
  } else {
    lines.push(`  None`);
  }
  lines.push(``);
  lines.push(`BETA-STRAND REGIONS (S):`);
  if(res.strand_segments.length){
    res.strand_segments.forEach(s =>
      lines.push(`  ${s.start}-${s.end} len=${s.length} ${s.seq}`)
    );
  } else {
    lines.push(`  None`);
  }
  lines.push(``);
  lines.push(`CONFLICTS & RESOLUTION:`);
  if(res.conflicts.length){
    res.conflicts.forEach(c =>
      lines.push(`  ${c.start}-${c.end} len=${c.length} avgPa=${c.avg_pa} avgPb=${c.avg_pb} -> ${c.chosen}  ${c.seq}`)
    );
  } else {
    lines.push(`  None`);
  }
  lines.push(``);
  lines.push(`Blocks (60-res):`);
  res.blocks.forEach(b => {
    lines.push(``);
    lines.push(b.range);
    lines.push(b.seq);
    lines.push(b.ann);
  });
  return lines.join('\n');
}

els.run.addEventListener('click', () => {
  const raw = els.seq.value;
  const method = els.method.value;
  const seq = cleanSeq(raw);

  if(!seq){
    alert("Please paste a valid amino-acid sequence (letters only).");
    return;
  }

  const res = method === 'chou_fasman'
    ? predictChouFasman(seq)
    : predictGOR(seq);

  els.summary.textContent =
    `Sequence length: ${res.length} • Method: ${
      method === 'chou_fasman' ? 'Chou–Fasman' : 'GOR (placeholder)'
    }`;

  renderBlocks(els.blocks, res.blocks);
  renderList(els.helixList, res.helix_segments);
  renderList(els.strandList, res.strand_segments);
  renderList(els.conflicts, res.conflicts);

  els.results.classList.remove('hidden');

  const report = buildTextReport(res.clean_seq, method, res);

  els.copyBtn.onclick = () => {
    navigator.clipboard.writeText(report).catch(() => {});
  };

  els.dlBtn.onclick = () => {
    const blob = new Blob([report], {type:'text/plain'});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `prediction_${method}.txt`;
    a.click();
    URL.revokeObjectURL(a.href);
  };
});