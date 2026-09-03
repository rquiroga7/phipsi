let viewer, model;
let currentPhi = 165, currentPsi = 165;
let animating = false;
let planeShapes = [];
let clashShapes = [];
let alanineShapes = [];
let showVDW = false;
let showPlanes = false;
let showClashes = false;
let trailClashes = false;
let clashesOrange = false; // original var
// For 3Dmol we keep one model; PDB has 2 models but we use only first
let atomsCache = []; // flat list of atom objects from model
const PDB_FALLBACK = `MODEL     1
ATOM      1  CA  LYS    15      -7.384  -4.928   3.141  1.00  0.00           C
ATOM      2  C   LYS    15      -8.891  -4.877   3.059  1.00  0.00           C
ATOM      3  O   LYS    15      -9.587  -5.900   3.158  1.00  0.00           O
ATOM      4  N   ALA    16      -9.429  -3.683   2.895  1.00  0.00           N
ATOM      5  CA  ALA    16     -10.875  -3.493   2.809  1.00  0.00           C
ATOM      6  C   ALA    16     -11.258  -2.040   2.980  1.00  0.00           C
ATOM      7  O   ALA    16     -10.470  -1.201   3.435  1.00  0.00           O
ATOM      8  CB  ALA    16     -11.526  -4.406   3.864  1.00  0.00           C
ATOM      9  H   ALA    16      -8.750  -2.871   2.796  1.00  0.00           H
ATOM     10  HA  ALA    16     -11.202  -3.798   1.797  1.00  0.00           H
ATOM     11 1HB  ALA    16     -11.267  -5.470   3.700  1.00  0.00           H
ATOM     12 2HB  ALA    16     -11.207  -4.150   4.892  1.00  0.00           H
ATOM     13 3HB  ALA    16     -12.629  -4.343   3.839  1.00  0.00           H
ATOM     14  N   ARG    17     -12.489  -1.722   2.616  1.00  0.00           N
ATOM     15  CA  ARG    17     -12.987  -0.352   2.720  1.00  0.00           C
ATOM     16  H   ARG    17     -13.068  -2.515   2.211  1.00  0.00           H
ENDMDL
`;

function initViewer(){
  const elem = document.getElementById('viewer');
  // original background #d0d0d0
  viewer = $3Dmol.createViewer(elem, {backgroundColor: '#d0d0d0'});
  // handle resize
  window.addEventListener('resize', ()=> viewer.resize());
}

function loadTripeptide(){
  // try local pdb first, fallback to 1CRN if fails
  fetch('phipsi-16atoms.pdb').then(r=>{
    if(!r.ok) throw new Error('no local');
    return r.text();
  }).then(txt=>{
    // PDB contains MODEL 1 and MODEL 2 – keep only MODEL 1 (first 18 lines)
    let pdbText = txt;
    if(txt.includes('MODEL')){
      // extract first MODEL ... ENDMDL
      const m = txt.match(/MODEL\s+1([\s\S]*?)ENDMDL/);
      if(m) pdbText = m[1];
    }
    viewer.clear();
    // clear shapes
    planeShapes=[]; clashShapes=[]; alanineShapes=[];
    model = viewer.addModel(pdbText,'pdb');
    // assign bonds if not present
    // Apply original styling
    applyOriginalStyle();
    viewer.zoomTo();
    // center as per original pp1.spt: center around -10.8 -3.48 3.12, rotate, zoom 150, set rotationRadius 5.78
    // 3Dmol zoomTo does similar; we can set view
    viewer.render();
    cacheAtoms();
    // compute initial phi/psi from coordinates (should be ~165)
    const v = getPhiPsi();
    if(!isNaN(v.phi)) currentPhi = v.phi;
    if(!isNaN(v.psi)) currentPsi = v.psi;
    updateDihedralsDisplay();
    drawRamaHeat();
    // highlight phi initially (green)
    highlightPhiPsi('phi');
    // setup initial rock? not needed
  }).catch(e=>{
    console.warn('Failed local PDB, using embedded fallback',e);
    let pdbText = PDB_FALLBACK;
    viewer.clear();
    planeShapes=[]; clashShapes=[]; alanineShapes=[];
    model = viewer.addModel(pdbText,'pdb');
    applyOriginalStyle();
    viewer.zoomTo();
    viewer.render();
    cacheAtoms();
    const v = getPhiPsi();
    if(!isNaN(v.phi)) currentPhi = v.phi;
    if(!isNaN(v.psi)) currentPsi = v.psi;
    updateDihedralsDisplay();
    drawRamaHeat();
    highlightPhiPsi('phi');
  });
}

function loadPDB(pdbId='1CRN'){
  // kept for backward compat with old page (resi-input)
  fetch(`https://files.rcsb.org/download/${pdbId}.pdb`).then(r=>r.text()).then(txt=>{
    viewer.clear();
    planeShapes=[]; clashShapes=[]; alanineShapes=[];
    model = viewer.addModel(txt,'pdb');
    viewer.setStyle({}, {cartoon:{color:'spectrum'}});
    viewer.zoomTo();
    viewer.render();
    cacheAtoms();
    updateDihedralsDisplay();
    drawRamaHeat();
  }).catch(e=>alert('Failed to load PDB:'+e));
}

function cacheAtoms(){
  // model.selectedAtoms returns array with elem, atom, resi, resn, chain, x,y,z, serial etc.
  atomsCache = model.selectedAtoms({});
  // 3Dmol atoms are references to internal objects; moving them via x,y,z will affect rendering after setStyle? Need to ensure we keep reference.
}

function applyOriginalStyle(){
  // clear previous styles
  // wireframe 0.07, spacefill 0.3 (hydrogen 0.2), colors as per pp1.spt
  // Use stick + sphere
  // First, set all atoms to ball-and-stick white bonds
  viewer.setStyle({}, {stick:{radius:0.15, color:'white'}, sphere:{scale:0.28}});
  // per-element sphere colors
  // CA special gray #707070, other C #c8c8c8
  // We need to override by selection
  // C atoms (elem C)
  viewer.setStyle({elem:'C'}, {sphere:{color:'#c8c8c8'}, stick:{color:'#c8c8c8'}});
  // need to re-color CA specifically – 3DMol selection by atom name CA
  viewer.setStyle({atom:'CA'}, {sphere:{color:'#707070'}});
  viewer.setStyle({elem:'N'}, {sphere:{color:'#6580ff'}, stick:{color:'#6580ff'}});
  viewer.setStyle({elem:'O'}, {sphere:{color:'#ff6060'}, stick:{color:'#ff6060'}});
  viewer.setStyle({elem:'H'}, {sphere:{scale:0.2, color:'white'}});
  // bonds white globally – but stick color white is needed for bond cylinders
  // In 3Dmol, stick color is derived from atom color? Actually stick color per atom determines bond color gradient.
  // To enforce white bonds as original: set stick color white for all
  // We'll keep sphere colors as above, but stick white
  viewer.setStyle({}, {stick:{radius:0.12, color:'white'}});
  // Re-apply sphere colors after ( spheres remain)
  viewer.setStyle({elem:'C'}, {sphere:{color:'#c8c8c8'}});
  viewer.setStyle({atom:'CA'}, {sphere:{color:'#707070'}});
  viewer.setStyle({elem:'N'}, {sphere:{color:'#6580ff'}});
  viewer.setStyle({elem:'O'}, {sphere:{color:'#ff6060'}});
  viewer.setStyle({elem:'H'}, {sphere:{color:'white', scale:0.2}});
  // highlight current bond (phi initially)
  // will be done via highlightPhiPsi
  // partial double bonds: not easily rendered, keep as sticks
}

function highlightPhiPsi(which){
  // color bond green #80ff80 for selected, white for other
  // We'll approximate by adding colored cylinders over the bond
  // First remove previous highlights
  if(viewer) {
    // remove previous highlight shapes if any (tagged)
    // we stored them as highlightShapes
    if(window._highlightShapes){
      window._highlightShapes.forEach(s=>viewer.removeShape(s));
    }
    window._highlightShapes=[];
    const getAtom = (resi, name)=> findAtom(resi, name);
    let a1,b1;
    if(which==='phi'){
      a1 = getAtom(16,'N'); b1 = getAtom(16,'CA');
    } else {
      a1 = getAtom(16,'CA'); b1 = getAtom(16,'C');
    }
    if(a1 && b1){
      const shape = viewer.addCylinder({
        start:{x:a1.x,y:a1.y,z:a1.z},
        end:{x:b1.x,y:b1.y,z:b1.z},
        radius:0.18,
        color:'#80ff80',
        fromCap:1,toCap:1
      });
      window._highlightShapes.push(shape);
      viewer.render();
    }
    // update label colors
    const phiAngleEl = document.getElementById('phiangle');
    const psiAngleEl = document.getElementById('psiangle');
    const phiLabel = document.getElementById('phi-label');
    const psiLabel = document.getElementById('psi-label');
    if(which==='phi'){
      if(phiAngleEl){phiAngleEl.style.color='#00c800';phiAngleEl.style.fontWeight='bold';}
      if(psiAngleEl){psiAngleEl.style.color='#a0a0a0';psiAngleEl.style.fontWeight='normal';}
      if(phiLabel){phiLabel.style.color='#00c800';}
      if(psiLabel){psiLabel.style.color='#a0a0a0';}
    } else {
      if(phiAngleEl){phiAngleEl.style.color='#a0a0a0';phiAngleEl.style.fontWeight='normal';}
      if(psiAngleEl){psiAngleEl.style.color='#00c800';psiAngleEl.style.fontWeight='bold';}
      if(phiLabel){phiLabel.style.color='#a0a0a0';}
      if(psiLabel){psiLabel.style.color='#00c800';}
    }
  }
}

// ---------- atom helpers ----------
function findAtom(resi, name){
  const atoms = model.selectedAtoms({});
  // chain may be '' – ignore chain, just match resi and atom
  return atoms.find(a=>a.resi===resi && a.atom===name) || null;
}
function findAtomByChain(chain, resi, name){ // for generic
  const atoms = model.selectedAtoms({});
  return atoms.find(a=>a.chain===chain && a.resi===resi && a.atom===name) || null;
}

function dihedral(p1,p2,p3,p4){
  const v1 = vecSub(p2,p1);
  const v2 = vecSub(p3,p2);
  const v3 = vecSub(p4,p3);
  const n1 = cross(v1,v2);
  const n2 = cross(v2,v3);
  const y = dot(cross(n1,n2), normalize(v2));
  const x = dot(n1,n2);
  return Math.atan2(y,x) * 180/Math.PI;
}
function vecSub(a,b){return {x:a.x-b.x,y:a.y-b.y,z:a.z-b.z};}
function cross(a,b){return {x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x};}
function dot(a,b){return a.x*b.x+a.y*b.y+a.z*b.z;}
function normalize(a){const m=Math.hypot(a.x,a.y,a.z)||1;return {x:a.x/m,y:a.y/m,z:a.z/m};}
function getPhiPsi(){
  // For tripeptide: phi = C(LYS15)-N(ALA16)-CA(ALA16)-C(ALA16)
  // psi = N(ALA16)-CA(ALA16)-C(ALA16)-N(ARG17)
  const Cprev = findAtom(15,'C');
  const N = findAtom(16,'N');
  const CA = findAtom(16,'CA');
  const C = findAtom(16,'C');
  const Nnext = findAtom(17,'N');
  if(!N||!CA||!C) return {phi:NaN,psi:NaN};
  let phi = Cprev ? dihedral(Cprev,N,CA,C) : NaN;
  let psi = Nnext ? dihedral(N,CA,C,Nnext) : NaN;
  // also support generic chain fallback if tripeptide not found (e.g., 1CRN)
  if(isNaN(phi) && isNaN(psi)){
    // try generic via residueOrder (old logic)
    const chain = getFirstChain();
    const resi = parseInt(document.getElementById('resi-input')?.value || 16,10);
    const N2 = findAtomByChain(chain,resi,'N');
    const CA2 = findAtomByChain(chain,resi,'CA');
    const C2 = findAtomByChain(chain,resi,'C');
    if(N2&&CA2&&C2){
      const Cprev2 = findAtomByChain(chain,resi-1,'C');
      const Nnext2 = findAtomByChain(chain,resi+1,'N');
      if(Cprev2) phi = dihedral(Cprev2,N2,CA2,C2);
      if(Nnext2) psi = dihedral(N2,CA2,C2,Nnext2);
    }
  }
  return {phi,psi};
}
function getFirstChain(){
  const atoms = model.selectedAtoms({});
  if(!atoms.length) return 'A';
  return atoms[0].chain || 'A';
}
function updateDihedralsDisplay(){
  const v = getPhiPsi();
  if(!isNaN(v.phi)) currentPhi = v.phi;
  if(!isNaN(v.psi)) currentPsi = v.psi;
  const phiAngleEl = document.getElementById('phiangle');
  const psiAngleEl = document.getElementById('psiangle');
  const phiValEl = document.getElementById('phi-val');
  const psiValEl = document.getElementById('psi-val');
  if(phiAngleEl) phiAngleEl.innerHTML = (isNaN(v.phi)?'n/a':v.phi.toFixed(1)+'&deg;');
  if(psiAngleEl) psiAngleEl.innerHTML = (isNaN(v.psi)?'n/a':v.psi.toFixed(1)+'&deg;');
  if(phiValEl) phiValEl.textContent = isNaN(v.phi)?'n/a':v.phi.toFixed(1);
  if(psiValEl) psiValEl.textContent = isNaN(v.psi)?'n/a':v.psi.toFixed(1);
  // also update legacy phiangle/psiangle colors handled elsewhere
  updatePlotMarker(v.phi, v.psi);
  // re-draw highlight position if needed (bond moved, need to update cylinder)
  // we keep highlight but positions are stale after rotation – remove and re-add on next highlight call?
  // For now, after each rotation we re-highlight same bond
}

function rotateAtomsAboutAxis(axisPoint1, axisPoint2, angleDeg, predicate){
  const angle = angleDeg*Math.PI/180;
  const u = {x:axisPoint2.x-axisPoint1.x,y:axisPoint2.y-axisPoint1.y,z:axisPoint2.z-axisPoint1.z};
  const un = normalize(u);
  const cosA = Math.cos(angle), sinA=Math.sin(angle);
  const atoms = model.selectedAtoms({});
  atoms.forEach(a=>{
    if(!predicate(a)) return;
    const p = {x:a.x-axisPoint1.x,y:a.y-axisPoint1.y,z:a.z-axisPoint1.z};
    const dotu = dot(un,p);
    const crossu = cross(un,p);
    const rotated = {
      x: un.x*dotu*(1-cosA) + p.x*cosA + crossu.x*sinA,
      y: un.y*dotu*(1-cosA) + p.y*cosA + crossu.y*sinA,
      z: un.z*dotu*(1-cosA) + p.z*cosA + crossu.z*sinA
    };
    a.x = rotated.x + axisPoint1.x;
    a.y = rotated.y + axisPoint1.y;
    a.z = rotated.z + axisPoint1.z;
  });
  // need to update model? In 3Dmol, atoms are references, but viewer needs to update?
  // We must tell model to update? The atoms array is live; after changing, we need to render.
}

function isMovingAtomForPhi(atom){
  // phi: rotation about N-CA of Ala16, moving atoms are C-terminal side:
  // includes: CA? (on axis), C, O, CB, HA, etc of residue 16 except N and H (amide), plus all of residue 17
  // For simplicity:
  if(atom.resi===15) return false; // Lys upstream
  if(atom.resi===16){
    // N and H are on N side, stay
    if(atom.atom==='N' || atom.atom==='H') return false;
    // everything else of 16 moves (CA is on axis but okay, C,O,CB, HA, HB etc)
    return true;
  }
  if(atom.resi===17) return true;
  // for generic fallback (1CRN with resi param)
  return false;
}
function isMovingAtomForPsi(atom){
  // psi: about CA-C of Ala16, moving: O of 16 + all of 17
  // CA side stays: N, CA, CB, HA of 16 stay
  if(atom.resi===15) return false;
  if(atom.resi===16){
    // only O (and maybe C? C is on axis) moves
    if(atom.atom==='O' || atom.atom==='C') {
      // C is on axis, movement zero but include
      return true;
    }
    return false;
  }
  if(atom.resi===17) return true;
  return false;
}
function genericMovingPredicate(chain, resi, type){
  // for fallback generic protein (1CRN) when tripeptide not found
  return (a)=>{
    if(a.chain!==chain) return false;
    if(type==='phi'){
      // phi: rotate residues >= resi, but exclude N/H as above
      if(a.resi < resi) return false;
      if(a.resi===resi && (a.atom==='N' || a.atom==='H' || a.atom==='H1' || a.atom==='H2' || a.atom==='H3')) return false;
      return true;
    } else {
      // psi: rotate residues > resi plus O of resi
      if(a.resi > resi) return true;
      if(a.resi===resi && (a.atom==='O' || a.atom==='OXT')) return true;
      return false;
    }
  };
}

function adjustDihedral(type, deltaDeg){
  if(animating) return;
  // determine which atoms define axis and predicate
  // try tripeptide atoms first
  const N = findAtom(16,'N');
  const CA = findAtom(16,'CA');
  const C = findAtom(16,'C');
  let axisA, axisB, predicate;
  let isTripeptide = !!(N && CA && C);
  if(isTripeptide){
    if(type==='phi'){
      if(!N||!CA) return;
      axisA=N; axisB=CA;
      predicate = isMovingAtomForPhi;
      highlightPhiPsi('phi');
    } else {
      if(!CA||!C) return;
      axisA=CA; axisB=C;
      predicate = isMovingAtomForPsi;
      highlightPhiPsi('psi');
    }
    // animate with steps of 2° (original) but for 15° => 8 steps of ~1.875 or 10 steps of 1.5
    // use 7 steps of variable? We'll use 8 steps for smoothness
    animating=true;
    animateRotation(axisA,axisB,deltaDeg,predicate,()=>{
      animating=false;
      // after rotation, need to reapply styles? Atoms moved, but styles remain
      // Update vdW clashes/planes if visible
      if(showPlanes) updatePlanes();
      if(showClashes) updateClashes();
      updateDihedralsDisplay();
      viewer.render();
      // re-highlight to update cylinder position
      highlightPhiPsi(type);
    });
  } else {
    // generic fallback
    const chain = getFirstChain();
    const resiEl = document.getElementById('resi-input');
    const resi = resiEl ? parseInt(resiEl.value,10) : 1;
    const N2 = findAtomByChain(chain,resi,'N');
    const CA2 = findAtomByChain(chain,resi,'CA');
    const C2 = findAtomByChain(chain,resi,'C');
    if(type==='phi'){
      if(!N2||!CA2){animating=false;return}
      axisA=N2; axisB=CA2;
      predicate = genericMovingPredicate(chain,resi,'phi');
      highlightPhiPsi('phi');
    } else {
      if(!CA2||!C2){animating=false;return}
      axisA=CA2; axisB=C2;
      predicate = genericMovingPredicate(chain,resi,'psi');
      highlightPhiPsi('psi');
    }
    animating=true;
    animateRotation(axisA,axisB,deltaDeg,predicate,()=>{
      animating=false;
      if(showPlanes) updatePlanes();
      if(showClashes) updateClashes();
      updateDihedralsDisplay();
      viewer.render();
      highlightPhiPsi(type);
    });
  }
}

function rotateDirect(type, deltaDeg, cb){
  const N = findAtom(16,'N');
  const CA = findAtom(16,'CA');
  const C = findAtom(16,'C');
  let isTripeptide = !!(N && CA && C);
  let axisA, axisB, predicate;
  if(isTripeptide){
    if(type==='phi'){
      axisA=N; axisB=CA; predicate=isMovingAtomForPhi;
    } else {
      axisA=CA; axisB=C; predicate=isMovingAtomForPsi;
    }
    animateRotation(axisA,axisB,deltaDeg,predicate,cb);
  } else {
    const chain = getFirstChain();
    const resi = parseInt(document.getElementById('resi-input')?.value||16,10);
    const N2 = findAtomByChain(chain,resi,'N');
    const CA2 = findAtomByChain(chain,resi,'CA');
    const C2 = findAtomByChain(chain,resi,'C');
    if(type==='phi'){
      axisA=N2; axisB=CA2; predicate=genericMovingPredicate(chain,resi,'phi');
    } else {
      axisA=CA2; axisB=C2; predicate=genericMovingPredicate(chain,resi,'psi');
    }
    if(!axisA||!axisB){ if(cb) cb(); return; }
    animateRotation(axisA,axisB,deltaDeg,predicate,cb);
  }
}

function animateRotation(axisA,axisB,totalDeg,predicate,cb){
  // replicate original: step 2° for smooth, here totalDeg may be 15 => 7-8 steps
  const steps = 10; // original had delta/step =2 => 10 steps for 20°
  // For 15°, steps=10 gives 1.5° per step – still smooth
  const step = totalDeg/steps;
  let i=0;
  // keep references to axis points at start; they move? Axis points are on static side, so not moving – use initial positions
  // However axis atoms themselves may be considered moving? For phi, CA is on axis but also considered moving? But axis stays.
  // So we keep axisPoints as objects (they won't move because predicate excludes them? N is not moving for phi, CA is moving? but CA is on axis so its position irrelevant)
  // To avoid drift, we copy axis coordinates
  const aA = {x:axisA.x, y:axisA.y, z:axisA.z};
  const aB = {x:axisB.x, y:axisB.y, z:axisB.z};
  function stepFn(){
    if(i<steps){
      rotateAtomsAboutAxis(aA,aB,step,predicate);
      // need to update shapes that depend on positions? planes/clashes handled after
      viewer.render();
      i++;
      setTimeout(stepFn,40); // ~40ms per step => 400ms total, similar to original delay
    } else {
      // re-cache? not needed
      if(cb) cb();
    }
  }
  stepFn();
}

// ---------- Controls: Alanine, Peptide Bonds, Planes, vdW, Clashes ----------
function doAlanine(checked){
  // original dots 60% black for Ala16
  // remove previous
  alanineShapes.forEach(s=>viewer.removeShape(s));
  alanineShapes=[];
  if(checked){
    // Add black dots (small spheres) around Ala16 atoms? Use surface dots
    // Approximate with translucent black spheres scaled 0.9 for each atom of Ala16?
    const alaAtoms = model.selectedAtoms({resi:16});
    alaAtoms.forEach(a=>{
      const s = viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:0.45, color:'black', opacity:0.6});
      alanineShapes.push(s);
    });
    viewer.render();
  } else {
    viewer.render();
  }
}
function doPeptideBonds(checked){
  // color peptide bonds magenta #ff80ff vs white
  // We'll add/remove magenta cylinders over peptide bonds
  if(window._pepShapes){
    window._pepShapes.forEach(s=>viewer.removeShape(s));
    window._pepShapes=[];
  } else window._pepShapes=[];
  if(checked){
    const pairs = [
      [findAtom(15,'C'), findAtom(16,'N')],
      [findAtom(16,'C'), findAtom(17,'N')],
      [findAtom(15,'C'), findAtom(15,'O')],
      [findAtom(16,'C'), findAtom(16,'O')]
    ];
    pairs.forEach(pair=>{
      const a=pair[0], b=pair[1];
      if(a&&b){
        const s=viewer.addCylinder({start:{x:a.x,y:a.y,z:a.z}, end:{x:b.x,y:b.y,z:b.z}, radius:0.12, color:'#ff80ff', fromCap:1,toCap:1});
        window._pepShapes.push(s);
      }
    });
    viewer.render();
  } else {
    viewer.render();
  }
}
function clearPlanes(){
  planeShapes.forEach(s=>viewer.removeShape(s));
  planeShapes=[];
  // also remove custom
  if(window._customPlanes){
    window._customPlanes.forEach(s=>viewer.removeShape(s));
    window._customPlanes=[];
  }
}
function updatePlanes(){
  if(!showPlanes) return;
  clearPlanes();
  // define plane corners: original hPlane {15.CA} {16.H} {16.CA} {15.O}
  // and sPlane {16.CA} {17.H} {17.CA} {16.O}
  const p1a = findAtom(15,'CA'), p1b=findAtom(16,'H'), p1c=findAtom(16,'CA'), p1d=findAtom(15,'O');
  const p2a = findAtom(16,'CA'), p2b=findAtom(17,'H'), p2c=findAtom(17,'CA'), p2d=findAtom(16,'O');
  const drawPlane = (a,b,c,d, color)=>{
    if(!a||!b||!c||!d) return;
    // Create two triangles: a-b-c and a-c-d
    const verts = [
      new $3Dmol.Vector3(a.x,a.y,a.z),
      new $3Dmol.Vector3(b.x,b.y,b.z),
      new $3Dmol.Vector3(c.x,c.y,c.z),
      new $3Dmol.Vector3(d.x,d.y,d.z)
    ];
    // normals: compute for each triangle
    const n1 = cross(vecSub(b,a), vecSub(c,a));
    const n2 = cross(vecSub(c,a), vecSub(d,a));
    const nn1 = normalize(n1), nn2=normalize(n2);
    const normals = [
      new $3Dmol.Vector3(nn1.x,nn1.y,nn1.z),
      new $3Dmol.Vector3(nn1.x,nn1.y,nn1.z),
      new $3Dmol.Vector3(nn1.x,nn1.y,nn1.z),
      new $3Dmol.Vector3(nn2.x,nn2.y,nn2.z)
    ];
    // we need 4 vertices, but faces reference 0,1,2 and 0,2,3
    const colors = [{r:0.6,g:0.6,b:1},{r:0.6,g:0.6,b:1},{r:0.6,g:0.6,b:1},{r:0.6,g:0.6,b:1}];
    // use addCustom
    try{
      const shape = viewer.addCustom({vertexArr:verts, normalArr:normals, faceArr:[0,1,2,0,2,3], color:colors});
      // set opacity via material? addCustom may not support opacity directly, we can set alpha in color?
      // For now use line border + custom
      window._customPlanes = window._customPlanes||[];
      window._customPlanes.push(shape);
    } catch(e){
      // fallback: draw border lines
      const col = color || '#a0a0ff';
      const l1 = viewer.addLine({start:{x:a.x,y:a.y,z:a.z}, end:{x:b.x,y:b.y,z:b.z}, color:col});
      const l2 = viewer.addLine({start:{x:b.x,y:b.y,z:b.z}, end:{x:c.x,y:c.y,z:c.z}, color:col});
      const l3 = viewer.addLine({start:{x:c.x,y:c.y,z:c.z}, end:{x:d.x,y:d.y,z:d.z}, color:col});
      const l4 = viewer.addLine({start:{x:d.x,y:d.y,z:d.z}, end:{x:a.x,y:a.y,z:a.z}, color:col});
      planeShapes.push(l1,l2,l3,l4);
    }
    // add border lines regardless for visibility
    const border = viewer.addLine({start:{x:a.x,y:a.y,z:a.z}, end:{x:b.x,y:b.y,z:b.z}, color:color});
    planeShapes.push(border);
  };
  // fallback drawing with simple lines if addCustom fails
  // Instead draw planes as translucent using addSphere? We'll try custom first
  // For now we will draw border and also add two triangles as lines
  // Use addCustom for each plane if possible
  // Try to create planes via addCustom; if fails, use lines
  try{
    if(p1a&&p1b&&p1c&&p1d){
      const verts1 = [new $3Dmol.Vector3(p1a.x,p1a.y,p1a.z), new $3Dmol.Vector3(p1b.x,p1b.y,p1b.z), new $3Dmol.Vector3(p1c.x,p1c.y,p1c.z), new $3Dmol.Vector3(p1d.x,p1d.y,p1d.z)];
      const nA = normalize(cross(vecSub(p1b,p1a), vecSub(p1c,p1a)));
      const normals1 = [new $3Dmol.Vector3(nA.x,nA.y,nA.z),new $3Dmol.Vector3(nA.x,nA.y,nA.z),new $3Dmol.Vector3(nA.x,nA.y,nA.z),new $3Dmol.Vector3(nA.x,nA.y,nA.z)];
      const s1 = viewer.addCustom({vertexArr:verts1, normalArr:normals1, faceArr:[0,1,2,0,2,3]});
      window._customPlanes = window._customPlanes||[]; window._customPlanes.push(s1);
    }
    if(p2a&&p2b&&p2c&&p2d){
      const verts2 = [new $3Dmol.Vector3(p2a.x,p2a.y,p2a.z), new $3Dmol.Vector3(p2b.x,p2b.y,p2b.z), new $3Dmol.Vector3(p2c.x,p2c.y,p2c.z), new $3Dmol.Vector3(p2d.x,p2d.y,p2d.z)];
      const nB = normalize(cross(vecSub(p2b,p2a), vecSub(p2c,p2a)));
      const normals2 = [new $3Dmol.Vector3(nB.x,nB.y,nB.z),new $3Dmol.Vector3(nB.x,nB.y,nB.z),new $3Dmol.Vector3(nB.x,nB.y,nB.z),new $3Dmol.Vector3(nB.x,nB.y,nB.z)];
      const s2 = viewer.addCustom({vertexArr:verts2, normalArr:normals2, faceArr:[0,1,2,0,2,3]});
      window._customPlanes = window._customPlanes||[]; window._customPlanes.push(s2);
    }
  } catch(e){
    // add lines as fallback
  }
  // Always add border lines for visibility
  if(p1a&&p1b&&p1c&&p1d){
    planeShapes.push(viewer.addLine({start:{x:p1a.x,y:p1a.y,z:p1a.z}, end:{x:p1b.x,y:p1b.y,z:p1b.z}, color:'#8080ff'}));
    planeShapes.push(viewer.addLine({start:{x:p1b.x,y:p1b.y,z:p1b.z}, end:{x:p1c.x,y:p1c.y,z:p1c.z}, color:'#8080ff'}));
    planeShapes.push(viewer.addLine({start:{x:p1c.x,y:p1c.y,z:p1c.z}, end:{x:p1d.x,y:p1d.y,z:p1d.z}, color:'#8080ff'}));
    planeShapes.push(viewer.addLine({start:{x:p1d.x,y:p1d.y,z:p1d.z}, end:{x:p1a.x,y:p1a.y,z:p1a.z}, color:'#8080ff'}));
  }
  if(p2a&&p2b&&p2c&&p2d){
    planeShapes.push(viewer.addLine({start:{x:p2a.x,y:p2a.y,z:p2a.z}, end:{x:p2b.x,y:p2b.y,z:p2b.z}, color:'#ffb060'}));
    planeShapes.push(viewer.addLine({start:{x:p2b.x,y:p2b.y,z:p2b.z}, end:{x:p2c.x,y:p2c.y,z:p2c.z}, color:'#ffb060'}));
    planeShapes.push(viewer.addLine({start:{x:p2c.x,y:p2c.y,z:p2c.z}, end:{x:p2d.x,y:p2d.y,z:p2d.z}, color:'#ffb060'}));
    planeShapes.push(viewer.addLine({start:{x:p2d.x,y:p2d.y,z:p2d.z}, end:{x:p2a.x,y:p2a.y,z:p2a.z}, color:'#ffb060'}));
  }
  viewer.render();
}
function doPlanes(checked){
  showPlanes = checked;
  if(checked){
    updatePlanes();
  } else {
    clearPlanes();
    viewer.render();
  }
}
function doVDW(checked){
  showVDW = checked;
  const whiteBox = document.getElementById('divwhite');
  if(checked){
    // show white checkbox
    if(whiteBox) whiteBox.style.display="";
    // zoomto 1.0 110 equivalent: viewer zoom
    viewer.zoomTo();
    // set all atoms to spacefill 88%
    // In 3Dmol, sphere scale 1.0 = vdW; 0.88 = 88%
    viewer.setStyle({}, {stick:{hidden:true}, sphere:{scale:0.88, opacity:0.85, colorfunc:undefined}});
    // re-apply per-element colors with opacity
    viewer.setStyle({elem:'C'}, {sphere:{color:'#c8c8c8', opacity:0.6}});
    viewer.setStyle({atom:'CA'}, {sphere:{color:'#707070', opacity:0.6}});
    viewer.setStyle({elem:'N'}, {sphere:{color:'#6580ff', opacity:0.6}});
    viewer.setStyle({elem:'O'}, {sphere:{color:'#ff6060', opacity:0.6}});
    viewer.setStyle({elem:'H'}, {sphere:{color:'white', opacity:0.6, scale:0.88}});
    viewer.setBackgroundColor('#d0d0d0');
    viewer.render();
  } else {
    if(whiteBox){
      whiteBox.style.display="none";
      document.getElementById('idwhite').checked=false;
    }
    viewer.setBackgroundColor('#d0d0d0');
    applyOriginalStyle();
    viewer.render();
  }
}
function doWhite(checked){
  if(checked){
    viewer.setBackgroundColor('white');
    // make second model? In original, select 1.2 color white translucent 6
    // For our single model, set spheres white translucent
    viewer.setStyle({}, {sphere:{color:'white', opacity:0.4}});
    // need to keep element colors? Original sets 1.2 white translucent 6 – that is second model copy for VDW
    // We'll just make background white and keep vdW translucent
    viewer.render();
  } else {
    viewer.setBackgroundColor('#d0d0d0');
    // restore vdW colors
    if(showVDW){
      doVDW(true);
    } else {
      viewer.render();
    }
  }
}
function clearClashes(){
  if(!trailClashes){
    clashShapes.forEach(s=>viewer.removeShape(s));
    clashShapes=[];
  }
}
function updateClashes(){
  if(!showClashes) return;
  if(!trailClashes){
    clashShapes.forEach(s=>viewer.removeShape(s));
    clashShapes=[];
  }
  // compute clashes: check pairs as per original makeClash calls but simplified to all non-bonded <88%
  const atoms = model.selectedAtoms({});
  const radii = {H:1.2, C:1.7, N:1.55, O:1.52, S:1.8};
  const scale = 0.88;
  // list of pairs to check – original checks:
  // "{1.1 and sidechain} {ca1}" where ca1 = {15.o,16.h,16.o,17.h}
  // "{1.1 and 15.o} {1.1 and (17.h or 17.n)}"
  // "{1.1 and 16.o} {1.1 and 15.o}"
  // "{1.1 and 16.h} {1.1 and 17.h}"
  // "{1.1 and 16.ha} {1.1 and (15.o or 17.h)}"
  // We'll just check all pairs distance < sum*0.88 and not directly bonded (distance >1.0)
  for(let i=0;i<atoms.length;i++){
    for(let j=i+1;j<atoms.length;j++){
      const a=atoms[i], b=atoms[j];
      // skip if same residue and bonded? approximate bond length <1.6
      const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
      const d=Math.hypot(dx,dy,dz);
      if(d<1.0) continue; // bonded or same
      if(d>4.0) continue;
      // estimate radii by elem (default 1.5)
      const ra = radii[a.elem] || 1.5;
      const rb = radii[b.elem] || 1.5;
      const sum = (ra+rb)*scale;
      if(d < sum){
        // clash – show orange sphere at midpoint (original orange if clashesOrange else default)
        const mid = {x:(a.x+b.x)/2, y:(a.y+b.y)/2, z:(a.z+b.z)/2};
        const color = clashesOrange ? 'orange' : '#ff8000';
        const s = viewer.addSphere({center:mid, radius:0.25, color:color, opacity:0.9});
        clashShapes.push(s);
        // also add line?
        const l = viewer.addLine({start:{x:a.x,y:a.y,z:a.z}, end:{x:b.x,y:b.y,z:b.z}, color:color, linewidth:2, dashed:true});
        clashShapes.push(l);
      }
    }
  }
  viewer.render();
}
function doClashes(checked){
  showClashes = checked;
  const trailDiv = document.getElementById('divtrailclashes');
  if(checked){
    if(trailDiv && document.body.contains(trailDiv)){
      // offer trail only if original would? Always show?
      trailDiv.style.display="";
    }
    updateClashes();
  } else {
    clearClashes();
    // hide trail
    if(trailDiv) trailDiv.style.display="none";
    document.getElementById('idtrailclashes').checked=false;
    trailClashes=false;
    clashShapes.forEach(s=>viewer.removeShape(s)); clashShapes=[];
    viewer.render();
  }
}
function doTrailClashes(checked){
  trailClashes = checked;
  if(checked){
    // keep existing, just set flag
  } else {
    // remove all but last?
    clashShapes.forEach(s=>viewer.removeShape(s)); clashShapes=[];
    if(showClashes) updateClashes();
    else viewer.render();
  }
}
function resetViewer(){
  // reload tripeptide and reset checkboxes
  document.getElementById('idalanine').checked=false;
  document.getElementById('idpeptidebonds').checked=false;
  document.getElementById('idplanes').checked=false;
  document.getElementById('idvdw').checked=false;
  document.getElementById('idwhite').checked=false;
  document.getElementById('idclashes').checked=false;
  document.getElementById('idtrailclashes').checked=false;
  document.getElementById('divwhite').style.display="none";
  document.getElementById('divtrailclashes').style.display="none";
  showPlanes=false; showVDW=false; showClashes=false; trailClashes=false;
  clearPlanes();
  clashShapes.forEach(s=>viewer.removeShape(s)); clashShapes=[];
  alanineShapes.forEach(s=>viewer.removeShape(s)); alanineShapes=[];
  if(window._pepShapes){window._pepShapes.forEach(s=>viewer.removeShape(s)); window._pepShapes=[];}
  if(window._highlightShapes){window._highlightShapes.forEach(s=>viewer.removeShape(s)); window._highlightShapes=[];}
  if(window._customPlanes){window._customPlanes.forEach(s=>viewer.removeShape(s)); window._customPlanes=[];}
  loadTripeptide();
}

// Ramachandran heatplot
function drawRamaHeat(){
  const size=121;
  const xs=[], ys=[];
  for(let i=0;i<size;i++){ xs.push(-180 + i*3); }
  for(let j=0;j<size;j++){ ys.push(-180 + j*3); }
  const z = [];
  for(let j=0;j<size;j++){
    const row=[];
    for(let i=0;i<size;i++){
      const x = xs[i], y=ys[j];
      let v=0;
      v += gaussian2(x,y,-60,-45,25,25); // alpha
      v += gaussian2(x,y,-120,130,30,28); // beta
      v += gaussian2(x,y,-120,-150,22,18); // beta lower
      v += gaussian2(x,y,60,40,20,20); // left-handed
      v += gaussian2(x,y,-70,140,18,18)*0.5;
      v += gaussian2(x,y,50,-130,15,15)*0.4;
      row.push(v);
    }
    z.push(row);
  }
  const data = [{
    x: xs, y: ys, z: z, type:'heatmap', colorscale:'YlOrRd', reversescale:false, zsmooth:'best', showscale:false, hoverinfo:'skip'
  }];
  const layout = {
    title:{text:'Ramachandran plot — click to animate', font:{size:14}},
    xaxis:{title:'Phi (°)', range:[-180,180], dtick:60, gridcolor:'#ddd', zeroline:true, zerolinecolor:'#999'},
    yaxis:{title:'Psi (°)', range:[-180,180], dtick:60, gridcolor:'#ddd', zeroline:true, zerolinecolor:'#999'},
    margin:{t:40,l:60,r:20,b:50},
    height:520,
    hovermode:'closest',
    plot_bgcolor:'#f8f8f8',
    paper_bgcolor:'white'
  };
  Plotly.newPlot('ramaplot', data, layout, {displayModeBar:false, responsive:true}).then(()=>{
    Plotly.addTraces('ramaplot', [
      {x:[currentPhi], y:[currentPsi], mode:'markers', marker:{color:'red', size:10, line:{color:'black', width:1}}, name:'current', hoverinfo:'skip', showlegend:false},
      {x:[currentPhi, currentPhi], y:[-180,180], mode:'lines', line:{color:'red', width:1, dash:'dot'}, hoverinfo:'none', showlegend:false},
      {x:[-180,180], y:[currentPsi, currentPsi], mode:'lines', line:{color:'red', width:1, dash:'dot'}, hoverinfo:'none', showlegend:false}
    ]);
    const plot = document.getElementById('ramaplot');
    plot.on('plotly_click', function(data){
      // data.points[0] is heatmap point; x is phi, y is psi
      // But heatmap click gives x/y via point's x/y? For heatmap, points have x and y fields
      let phi, psi;
      if(data.points && data.points[0]){
        phi = data.points[0].x;
        psi = data.points[0].y;
        // If clicking on marker traces, also handle
        if(data.points[0].data && data.points[0].data.mode==='markers'){
          phi = data.points[0].x; psi = data.points[0].y;
        }
      } else {
        // fallback from event
        const evt = data.event;
        // not precise
        phi = 0; psi=0;
      }
      if(typeof phi==='number' && typeof psi==='number'){
        // clamp
        phi = Math.max(-180, Math.min(180, phi));
        psi = Math.max(-180, Math.min(180, psi));
        moveToPhiPsi(phi, psi);
      }
    });
    // also allow click on plot area background – use restyle? Plotly's click on heatmap already handled.
  });
}
function gaussian2(x,y,muX,muY,sx,sy){
  return Math.exp(-((x-muX)*(x-muX)/(2*sx*sx) + (y-muY)*(y-muY)/(2*sy*sy)));
}
function updatePlotMarker(phi,psi){
  if(isNaN(phi)||isNaN(psi)) return;
  // check if plot exists
  const plot = document.getElementById('ramaplot');
  if(!plot || !plot.data) return;
  // traces 1,2,3 are marker and lines (indices)
  try{
    Plotly.restyle('ramaplot',{x:[[phi]],y:[[psi]]},[1]);
    Plotly.restyle('ramaplot',{x:[[phi,phi]],y:[[-180,180]]},[2]);
    Plotly.restyle('ramaplot',{x:[[-180,180]],y:[[psi,psi]]},[3]);
  } catch(e){
    // plot not yet ready
  }
}
function moveToPhiPsi(targetPhi,targetPsi){
  if(isNaN(targetPhi)||isNaN(targetPsi)) return;
  if(animating) return;
  const cur = getPhiPsi();
  if(isNaN(cur.phi)||isNaN(cur.psi)) return;
  let dPhi = normalizeAngle(targetPhi - cur.phi);
  let dPsi = normalizeAngle(targetPsi - cur.psi);
  // If delta is large, we animate in small steps; original does gradual
  // We'll do sequential small steps: iterate steps = max(|dPhi|,|dPsi|)/3
  const maxDelta = Math.max(Math.abs(dPhi), Math.abs(dPsi));
  const steps = Math.max(1, Math.ceil(maxDelta/4)); // 4° per big step
  let i=0;
  animating=true;
  function step(){
    if(i<steps){
      const pf = dPhi/steps;
      const ps = dPsi/steps;
      // rotate phi then psi in one combined animation frame
      // Do phi first, then psi
      rotateDirect('phi', pf, ()=>{
        rotateDirect('psi', ps, ()=>{
          const cur2 = getPhiPsi();
          updateDihedralsDisplay();
          // update highlight to latest?
          // highlight based on larger delta?
          if(Math.abs(dPhi) > Math.abs(dPsi)) highlightPhiPsi('phi'); else highlightPhiPsi('psi');
          i++;
          setTimeout(step, 80);
          if(i>=steps){
            animating=false;
            viewer.render();
          }
        });
      });
    } else {
      animating=false;
    }
  }
  step();
}
function normalizeAngle(a){
  while(a>180) a-=360;
  while(a<-180) a+=360;
  return a;
}

document.addEventListener('DOMContentLoaded',()=>{
  initViewer();
  loadTripeptide();
  // button listeners – all 4 must respond
  document.getElementById('phi-plus').addEventListener('click',()=> adjustDihedral('phi', 15));
  document.getElementById('phi-minus').addEventListener('click',()=> adjustDihedral('phi', -15));
  document.getElementById('psi-plus').addEventListener('click',()=> adjustDihedral('psi', 15));
  document.getElementById('psi-minus').addEventListener('click',()=> adjustDihedral('psi', -15));
  // legacy compat: if resi-input exists
  const resiInput = document.getElementById('resi-input');
  if(resiInput) resiInput.addEventListener('change', updateDihedralsDisplay);
  const loadBtn = document.getElementById('load-pdb');
  if(loadBtn) loadBtn.addEventListener('click',()=> loadPDB('1CRN'));
  // checkboxes
  document.getElementById('idalanine').addEventListener('change', e=> doAlanine(e.target.checked));
  document.getElementById('idpeptidebonds').addEventListener('change', e=> doPeptideBonds(e.target.checked));
  document.getElementById('idplanes').addEventListener('change', e=> doPlanes(e.target.checked));
  document.getElementById('idvdw').addEventListener('change', e=> doVDW(e.target.checked));
  document.getElementById('idwhite').addEventListener('change', e=> doWhite(e.target.checked));
  document.getElementById('idclashes').addEventListener('change', e=> doClashes(e.target.checked));
  document.getElementById('idtrailclashes').addEventListener('change', e=> doTrailClashes(e.target.checked));
  document.getElementById('resetBtn').addEventListener('click', resetViewer);
  // also handle window resize for plotly
  window.addEventListener('resize', ()=> {
    try{ Plotly.Plots.resize('ramaplot'); viewer.resize(); }catch(e){}
  });
});

