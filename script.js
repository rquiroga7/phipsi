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
let showPeptideBondsFlag = false;
let atomsData = []; // our source of truth for coordinates
let doubleBondShapes = [];
let vdwShapes = [];
let clashAtomSerials = new Set();
const PDB_FALLBACK = `ATOM      1  CA  LYS    15      -7.384  -4.928   3.141  1.00  0.00           C
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
`;

// helpers to parse / generate PDB
function parsePDB(pdbText){
  const atoms=[];
  pdbText.split('\n').forEach(line=>{
    if(line.startsWith('ATOM') || line.startsWith('HETATM')){
      const serial=parseInt(line.substring(6,11).trim());
      const name=line.substring(12,16).trim();
      const resn=line.substring(17,20).trim();
      const chain=line.substring(21,22).trim();
      const resi=parseInt(line.substring(22,26).trim());
      const x=parseFloat(line.substring(30,38));
      const y=parseFloat(line.substring(38,46));
      const z=parseFloat(line.substring(46,54));
      const elem=line.substring(76,78).trim() || name[0].replace(/[0-9]/g,'');
      if(!isNaN(x) && !isNaN(resi)) atoms.push({serial,name,resn,chain,resi,x,y,z,elem});
    }
  });
  return atoms;
}
function generatePDB(atoms){
  let out='';
  atoms.forEach((a,i)=>{
    const serial= String(a.serial||i+1).padStart(5);
    const name= (a.name.length<4? ' '+a.name : a.name).padEnd(4);
    const resn= (a.resn||'ALA').padEnd(3);
    const chain= (a.chain||' ').padEnd(1);
    const resi= String(a.resi).padStart(4);
    const x= a.x.toFixed(3).padStart(8);
    const y= a.y.toFixed(3).padStart(8);
    const z= a.z.toFixed(3).padStart(8);
    const elem= (a.elem||a.name[0]).padStart(2);
    out+=`ATOM  ${serial} ${name} ${resn} ${chain}${resi}    ${x}${y}${z}  1.00  0.00          ${elem}  \n`;
  });
  out+='END\n';
  return out;
}
function syncModelPositions(){
  if(!model) return;
  try{
    const m = viewer.getModel(0) || model;
    const mAtoms = m.atoms;
    if(mAtoms){
      for(let i=0;i<atomsData.length && i<mAtoms.length;i++){
        mAtoms[i].x=atomsData[i].x;
        mAtoms[i].y=atomsData[i].y;
        mAtoms[i].z=atomsData[i].z;
      }
    }
    if(m.frames && m.frames[0]){
      const f=m.frames[0];
      for(let i=0;i<atomsData.length && i<f.length;i++){
        f[i].x=atomsData[i].x;
        f[i].y=atomsData[i].y;
        f[i].z=atomsData[i].z;
      }
    }
  }catch(e){}
  // 3Dmol caches WebGL geometry built at setStyle time, so re-applying the style
  // (without viewer.clear) regenerates ball-and-stick from the updated coords.
  // This is much cheaper than rebuildModel and keeps all extra shapes.
  applyOriginalStyle();
}
function rebuildModel(preserveView){
  let view=null;
  if(preserveView && viewer) try{ view=viewer.getView(); }catch(e){}
  const hadPlanes=showPlanes, hadClash=showClashes, hadPep=showPeptideBondsFlag, hadAlanine=alanineShapes.length>0, hadVDW=showVDW, hadWhite=document.getElementById('idwhite')?.checked, hadHighlight=window._lastHighlight;
  const hadTrail=trailClashes;
  // keep clashShapes if trail, else clear
  if(!hadTrail) clashShapes=[];
  viewer.clear();
  doubleBondShapes=[]; planeShapes=[]; vdwShapes=[];
  // alanineShapes will be re-added via doAlanine, so clear array but keep flag
  alanineShapes=[];
  // viewer.clear removes all shapes, so need to clear arrays that held shape refs
  if(window._pepShapes) window._pepShapes=[];
  if(window._highlightShapes) window._highlightShapes=[];
  if(window._customPlanes) window._customPlanes=[];
  const pdbText=generatePDB(atomsData);
  model=viewer.addModel(pdbText,'pdb', {keepH:true});
  applyOriginalStyle();
  if(view) try{ viewer.setView(view); }catch(e){ viewer.zoomTo(); }
  else viewer.zoomTo();
  viewer.render();
  // re-add persistent shapes that were visible before
  if(hadAlanine) doAlanine(true);
  if(hadPep) doPeptideBonds(true);
  if(hadPlanes) updatePlanes();
  if(hadClash) updateClashes();
  if(hadVDW){
    // re-add VdW overlay at new positions
    const radii={H:1.20, C:1.70, N:1.55, O:1.52};
    const isWhite=hadWhite;
    vdwShapes=[];
    atomsData.forEach(a=>{
      const ra=radii[a.elem]||1.5;
      let col='#c8c8c8';
      if(isWhite) col='white';
      else {
        if(a.atom==='CA') col='#000000';
        else if(a.elem==='N') col='#3050ff';
        else if(a.elem==='O') col='#ff2020';
        else if(a.elem==='H') col='white';
      }
      const op=isWhite?0.28:0.38;
      const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:ra*0.88, color:col, opacity:op});
      vdwShapes.push(s);
    });
    viewer.render();
  }
  if(hadHighlight) highlightPhiPsi(hadHighlight);
}

function initViewer(){
  const elem=document.getElementById('viewer');
  viewer=$3Dmol.createViewer(elem,{backgroundColor:0xd0d0d0});
  viewer.setBackgroundColor(0xd0d0d0);
  setTimeout(()=>{ try{ viewer.resize(); viewer.render(); }catch(e){} }, 80);
  window.addEventListener('resize',()=>{ try{ viewer.resize(); viewer.render(); }catch(e){} });
}

function loadTripeptide(){
  fetch('phipsi-16atoms.pdb').then(r=>{
    if(!r.ok) throw new Error('no local');
    return r.text();
  }).then(txt=>{
    let pdbText=txt;
    if(txt.includes('MODEL')){
      const m=txt.match(/MODEL\s+1([\s\S]*?)ENDMDL/);
      if(m) pdbText=m[1];
    }
    atomsData=parsePDB(pdbText);
    if(atomsData.length===0) throw new Error('empty parse');
    viewer.clear();
    planeShapes=[]; clashShapes=[]; alanineShapes=[];
    model=viewer.addModel(generatePDB(atomsData),'pdb', {keepH:true});
    applyOriginalStyle();
    viewer.zoomTo();
    // apply initial view similar to original: rotate z -108, y 138, z -85
    try{
      viewer.rotate(-108.34,'z');
      viewer.rotate(138.09,'y');
      viewer.rotate(-85.48,'z');
    }catch(e){}
    viewer.render();
    const v=getPhiPsi();
    if(!isNaN(v.phi)) currentPhi=v.phi;
    if(!isNaN(v.psi)) currentPsi=v.psi;
    updateDihedralsDisplay();
    drawRamaHeat();
    highlightPhiPsi('phi');
  }).catch(e=>{
    console.warn('Failed local PDB, using embedded fallback',e);
    atomsData=parsePDB(PDB_FALLBACK);
    viewer.clear();
    planeShapes=[]; clashShapes=[]; alanineShapes=[];
    model=viewer.addModel(generatePDB(atomsData),'pdb', {keepH:true});
    applyOriginalStyle();
    viewer.zoomTo();
    try{
      viewer.rotate(-108.34,'z');
      viewer.rotate(138.09,'y');
      viewer.rotate(-85.48,'z');
    }catch(e){}
    viewer.render();
    const v=getPhiPsi();
    if(!isNaN(v.phi)) currentPhi=v.phi;
    if(!isNaN(v.psi)) currentPsi=v.psi;
    updateDihedralsDisplay();
    drawRamaHeat();
    highlightPhiPsi('phi');
  });
}

function loadPDB(pdbId='1CRN'){
  fetch(`https://files.rcsb.org/download/${pdbId}.pdb`).then(r=>r.text()).then(txt=>{
    atomsData=parsePDB(txt);
    viewer.clear();
    planeShapes=[]; clashShapes=[]; alanineShapes=[];
    model=viewer.addModel(generatePDB(atomsData),'pdb', {keepH:true});
    viewer.setStyle({}, {cartoon:{color:'spectrum'}});
    viewer.zoomTo();
    viewer.render();
    updateDihedralsDisplay();
    drawRamaHeat();
  }).catch(e=>alert('Failed to load PDB:'+e));
}

function applyOriginalStyle(){
  // base ball-and-stick: stick radius ~0.12, sphere scale 0.30, white bonds
  viewer.setStyle({}, {stick:{radius:0.12, color:'white'}, sphere:{scale:0.30}});
  // now override colors via addStyle to keep scale
  viewer.addStyle({elem:'C'}, {sphere:{color:'#c8c8c8', scale:0.30}});
  viewer.addStyle({atom:'CA'}, {sphere:{color:'#505050', scale:0.32}}); // Calpha darker gray #505050 than other C #c8c8c8
  viewer.addStyle({elem:'N'}, {sphere:{color:'#3050ff'}});
  viewer.addStyle({elem:'O'}, {sphere:{color:'#ff2020'}});
  viewer.setStyle({elem:'H'}, {sphere:{color:'white', scale:0.25, hidden:false}, stick:{color:'white', radius:0.12, hidden:false}});
  viewer.addStyle({atom:'H'}, {sphere:{color:'white', scale:0.25, hidden:false}});
  viewer.addStyle({atom:'HA'}, {sphere:{color:'white', scale:0.25, hidden:false}});
  viewer.addStyle({atom:'1HB'}, {sphere:{color:'white', scale:0.25, hidden:false}});
  viewer.addStyle({atom:'2HB'}, {sphere:{color:'white', scale:0.25, hidden:false}});
  viewer.addStyle({atom:'3HB'}, {sphere:{color:'white', scale:0.25, hidden:false}});
  viewer.addStyle({}, {stick:{color:'white', radius:0.12, hidden:false}});
  viewer.render();
  // add dotted partial double bonds for peptide and carbonyl (always visible, like original)
  updateDoubleBonds();
}
function updateDoubleBonds(){
  doubleBondShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
  doubleBondShapes=[];
  const pairs=[
    [findAtom(15,'C'), findAtom(16,'N')],
    [findAtom(16,'C'), findAtom(17,'N')],
    [findAtom(15,'C'), findAtom(15,'O')],
    [findAtom(16,'C'), findAtom(16,'O')]
  ];
  pairs.forEach(pair=>{
    const a=pair[0], b=pair[1];
    if(!a||!b) return;
    // original Jmol: set multiplebondspacing 0.15, partialdouble -> one solid + one dotted parallel offset by 0.15
    // we add a dashed line offset perpendicular to bond by 0.15
    const mx=(a.x+b.x)/2, my=(a.y+b.y)/2, mz=(a.z+b.z)/2;
    const bx=b.x-a.x, by=b.y-a.y, bz=b.z-a.z;
    // find perpendicular: cross bond with arbitrary up
    let ux=0, uy=0, uz=1;
    // if bond is near vertical, use other up
    if(Math.abs(bz) > 0.9*Math.hypot(bx,by,bz)){ ux=1; uy=0; uz=0; }
    let px=by*uz - bz*uy, py=bz*ux - bx*uz, pz=bx*uy - by*ux;
    const plen=Math.hypot(px,py,pz)||1;
    px=px/plen*0.15; py=py/plen*0.15; pz=pz/plen*0.15;
    const a2={x:a.x+px, y:a.y+py, z:a.z+pz}, b2={x:b.x+px, y:b.y+py, z:b.z+pz};
    const s=viewer.addLine({start:a2, end:b2, dashed:true, dashLength:0.10, gapLength:0.10, color:'white', linewidth:2});
    doubleBondShapes.push(s);
  });
  viewer.render();
}

function highlightPhiPsi(which){
  window._lastHighlight=which;
  if(!viewer) return;
  if(window._highlightShapes){
    window._highlightShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
  }
  window._highlightShapes=[];
  const getAtom=(resi,name)=> atomsData.find(a=>a.resi===resi && a.name===name) || null;
  let a1,b1;
  if(which==='phi'){ a1=getAtom(16,'N'); b1=getAtom(16,'CA'); }
  else { a1=getAtom(16,'CA'); b1=getAtom(16,'C'); }
  if(a1 && b1){
    const shape=viewer.addCylinder({
      start:{x:a1.x,y:a1.y,z:a1.z},
      end:{x:b1.x,y:b1.y,z:b1.z},
      radius:0.18,
      color:'#80ff80',
      fromCap:1,toCap:1
    });
    window._highlightShapes.push(shape);
    viewer.render();
  }
  const phiAngleEl=document.getElementById('phiangle');
  const psiAngleEl=document.getElementById('psiangle');
  const phiLabel=document.getElementById('phi-label');
  const psiLabel=document.getElementById('psi-label');
  if(which==='phi'){
    if(phiAngleEl){phiAngleEl.style.color='#00c800';phiAngleEl.style.fontWeight='bold';}
    if(psiAngleEl){psiAngleEl.style.color='#a0a0a0';psiAngleEl.style.fontWeight='normal';}
    if(phiLabel) phiLabel.style.color='#00c800';
    if(psiLabel) psiLabel.style.color='#a0a0a0';
  } else {
    if(phiAngleEl){phiAngleEl.style.color='#a0a0a0';phiAngleEl.style.fontWeight='normal';}
    if(psiAngleEl){psiAngleEl.style.color='#00c800';psiAngleEl.style.fontWeight='bold';}
    if(phiLabel) phiLabel.style.color='#a0a0a0';
    if(psiLabel) psiLabel.style.color='#00c800';
  }
}

function findAtom(resi,name){
  return atomsData.find(a=>a.resi===resi && a.name===name) || null;
}
function findAtomByChain(chain,resi,name){
  return atomsData.find(a=>a.chain===chain && a.resi===resi && a.name===name) || null;
}
function dihedral(p1,p2,p3,p4){
  const v1=vecSub(p2,p1), v2=vecSub(p3,p2), v3=vecSub(p4,p3);
  const n1=cross(v1,v2), n2=cross(v2,v3);
  const y=dot(cross(n1,n2), normalize(v2));
  const x=dot(n1,n2);
  return Math.atan2(y,x)*180/Math.PI;
}
function vecSub(a,b){return {x:a.x-b.x,y:a.y-b.y,z:a.z-b.z};}
function cross(a,b){return {x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x};}
function dot(a,b){return a.x*b.x+a.y*b.y+a.z*b.z;}
function normalize(a){const m=Math.hypot(a.x,a.y,a.z)||1;return {x:a.x/m,y:a.y/m,z:a.z/m};}
function getPhiPsi(){
  const Cprev=findAtom(15,'C'), N=findAtom(16,'N'), CA=findAtom(16,'CA'), C=findAtom(16,'C'), Nnext=findAtom(17,'N');
  if(!N||!CA||!C) return {phi:NaN,psi:NaN};
  let phi=Cprev?dihedral(Cprev,N,CA,C):NaN;
  let psi=Nnext?dihedral(N,CA,C,Nnext):NaN;
  if(isNaN(phi) && isNaN(psi)){
    const chain=getFirstChain();
    const resi=parseInt(document.getElementById('resi-input')?.value||16,10);
    const N2=findAtomByChain(chain,resi,'N'), CA2=findAtomByChain(chain,resi,'CA'), C2=findAtomByChain(chain,resi,'C');
    if(N2&&CA2&&C2){
      const Cprev2=findAtomByChain(chain,resi-1,'C'), Nnext2=findAtomByChain(chain,resi+1,'N');
      if(Cprev2) phi=dihedral(Cprev2,N2,CA2,C2);
      if(Nnext2) psi=dihedral(N2,CA2,C2,Nnext2);
    }
  }
  return {phi,psi};
}
function getFirstChain(){
  if(!atomsData.length) return 'A';
  return atomsData[0].chain||'A';
}
function updateDihedralsDisplay(){
  const v=getPhiPsi();
  if(!isNaN(v.phi)) currentPhi=v.phi;
  if(!isNaN(v.psi)) currentPsi=v.psi;
  const phiAngleEl=document.getElementById('phiangle'), psiAngleEl=document.getElementById('psiangle');
  const phiValEl=document.getElementById('phi-val'), psiValEl=document.getElementById('psi-val');
  if(phiAngleEl) phiAngleEl.innerHTML=(isNaN(v.phi)?'n/a':v.phi.toFixed(1)+'&deg;');
  if(psiAngleEl) psiAngleEl.innerHTML=(isNaN(v.psi)?'n/a':v.psi.toFixed(1)+'&deg;');
  if(phiValEl) phiValEl.textContent=isNaN(v.phi)?'n/a':v.phi.toFixed(1);
  if(psiValEl) psiValEl.textContent=isNaN(v.psi)?'n/a':v.psi.toFixed(1);
  updatePlotMarker(v.phi,v.psi);
}

function rotateAtomsAboutAxis(axisA, axisB, angleDeg, predicate){
  const angle=angleDeg*Math.PI/180;
  const u={x:axisB.x-axisA.x,y:axisB.y-axisA.y,z:axisB.z-axisA.z};
  const un=normalize(u);
  const cosA=Math.cos(angle), sinA=Math.sin(angle);
  atomsData.forEach(a=>{
    if(!predicate(a)) return;
    const p={x:a.x-axisA.x,y:a.y-axisA.y,z:a.z-axisA.z};
    const dotu=dot(un,p);
    const crossu=cross(un,p);
    const rot={
      x: un.x*dotu*(1-cosA)+p.x*cosA+crossu.x*sinA,
      y: un.y*dotu*(1-cosA)+p.y*cosA+crossu.y*sinA,
      z: un.z*dotu*(1-cosA)+p.z*cosA+crossu.z*sinA
    };
    a.x=rot.x+axisA.x;
    a.y=rot.y+axisA.y;
    a.z=rot.z+axisA.z;
  });
  syncModelPositions();
  updateDoubleBonds();
  if(window._lastHighlight) highlightPhiPsi(window._lastHighlight);
  if(showPlanes) updatePlanes();
  if(showPeptideBondsFlag) doPeptideBonds(true);
  if(showVDW){
    // update VdW spheres to new positions (keep overlay)
    vdwShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    vdwShapes=[];
    const radii={H:1.20, C:1.70, N:1.55, O:1.52};
    const isWhite=document.getElementById('idwhite')?.checked;
    atomsData.forEach(a=>{
      if(showClashes && clashAtomSerials.has(a.serial)) return; // lens replaces this atom's VdW
      const ra=radii[a.elem]||1.5;
      let col='#c8c8c8';
      if(isWhite) col='white';
      else {
        if(a.atom==='CA') col='#000000';
        else if(a.elem==='N') col='#3050ff';
        else if(a.elem==='O') col='#ff2020';
        else if(a.elem==='H') col='white';
      }
      const op=isWhite?0.28:0.38;
      const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:ra*0.88, color:col, opacity:showClashes?(isWhite?0.16:0.22):op});
      vdwShapes.push(s);
    });
    // if clashes also on, recolor after
    if(showClashes) updateClashes();
  } else if(showClashes){
    // VdW off but clashes on: ensure VdW for clash viz is updated
    updateClashes();
  }
  if(alanineShapes.length>0){
    alanineShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    alanineShapes=[];
    const alaAtoms=atomsData.filter(a=>a.resi===16);
    alaAtoms.forEach(a=>{
      const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:0.9, color:'black', opacity:0.18});
      alanineShapes.push(s);
    });
  }
}

function isMovingAtomForPhi(atom){
  // phi N-CA: rotate N-terminal peptide (hPlane blue: 15.CA/15.O/16.H) so phi moves Lys + H, psi moves C-terminal (sPlane orange)
  // This makes phi->blue, psi->orange distinct as per original Jmol (phi hPlane, psi sPlane)
  if(atom.resi===15) return true;
  if(atom.resi===16){
    if(atom.name==='H') return true;
    return false;
  }
  return false;
}
function isMovingAtomForPsi(atom){
  // psi rotates about CA-C, moves C-terminal plane (sPlane: 16.O,17.H,17.CA) plus downstream
  if(atom.resi===15) return false;
  if(atom.resi===16){
    if(atom.name==='O' || atom.name==='C') return true;
    return false;
  }
  if(atom.resi===17) return true;
  return false;
}
function genericMovingPredicate(chain,resi,type){
  return (a)=>{
    if(a.chain!==chain) return false;
    if(type==='phi'){
      if(a.resi<resi) return false;
      if(a.resi===resi && (a.name==='N' || a.name==='H')) return false;
      return true;
    } else {
      if(a.resi>resi) return true;
      if(a.resi===resi && (a.name==='O' || a.name==='OXT')) return true;
      return false;
    }
  };
}

function adjustDihedral(type,deltaDeg){
  if(animating) return;
  const N=findAtom(16,'N'), CA=findAtom(16,'CA'), C=findAtom(16,'C');
  let axisA, axisB, predicate;
  let isTri=!!(N&&CA&&C);
  if(isTri){
    // phi upstream needs inverted sign to get +15 phi
    if(type==='phi'){ axisA={x:N.x,y:N.y,z:N.z}; axisB={x:CA.x,y:CA.y,z:CA.z}; predicate=isMovingAtomForPhi; highlightPhiPsi('phi'); deltaDeg=-deltaDeg; }
    else { axisA={x:CA.x,y:CA.y,z:CA.z}; axisB={x:C.x,y:C.y,z:C.z}; predicate=isMovingAtomForPsi; highlightPhiPsi('psi'); }
    animating=true;
    animateRotation(axisA,axisB,deltaDeg,predicate,()=>{
      animating=false;
      if(showPlanes) updatePlanes();
      if(showClashes) updateClashes();
      updateDihedralsDisplay();
      highlightPhiPsi(type);
      viewer.render();
    });
  } else {
    const chain=getFirstChain();
    const resi=parseInt(document.getElementById('resi-input')?.value||16,10);
    const N2=findAtomByChain(chain,resi,'N'), CA2=findAtomByChain(chain,resi,'CA'), C2=findAtomByChain(chain,resi,'C');
    if(type==='phi'){
      if(!N2||!CA2) return;
      axisA={x:N2.x,y:N2.y,z:N2.z}; axisB={x:CA2.x,y:CA2.y,z:CA2.z}; predicate=genericMovingPredicate(chain,resi,'phi'); highlightPhiPsi('phi');
    } else {
      if(!CA2||!C2) return;
      axisA={x:CA2.x,y:CA2.y,z:CA2.z}; axisB={x:C2.x,y:C2.y,z:C2.z}; predicate=genericMovingPredicate(chain,resi,'psi'); highlightPhiPsi('psi');
    }
    animating=true;
    animateRotation(axisA,axisB,deltaDeg,predicate,()=>{
      animating=false;
      if(showPlanes) updatePlanes();
      if(showClashes) updateClashes();
      updateDihedralsDisplay();
      highlightPhiPsi(type);
      viewer.render();
    });
  }
}

function rotateDirect(type,deltaDeg,cb){
  const N=findAtom(16,'N'), CA=findAtom(16,'CA'), C=findAtom(16,'C');
  let isTri=!!(N&&CA&&C);
  let axisA,axisB,predicate;
  if(isTri){
    if(type==='phi'){ axisA={x:N.x,y:N.y,z:N.z}; axisB={x:CA.x,y:CA.y,z:CA.z}; predicate=isMovingAtomForPhi; deltaDeg=-deltaDeg; }
    else { axisA={x:CA.x,y:CA.y,z:CA.z}; axisB={x:C.x,y:C.y,z:C.z}; predicate=isMovingAtomForPsi; }
    animateRotation(axisA,axisB,deltaDeg,predicate,cb);
  } else {
    const chain=getFirstChain();
    const resi=parseInt(document.getElementById('resi-input')?.value||16,10);
    const N2=findAtomByChain(chain,resi,'N'), CA2=findAtomByChain(chain,resi,'CA'), C2=findAtomByChain(chain,resi,'C');
    if(type==='phi'){ axisA={x:N2.x,y:N2.y,z:N2.z}; axisB={x:CA2.x,y:CA2.y,z:CA2.z}; predicate=genericMovingPredicate(chain,resi,'phi'); }
    else { axisA={x:CA2.x,y:CA2.y,z:CA2.z}; axisB={x:C2.x,y:C2.y,z:C2.z}; predicate=genericMovingPredicate(chain,resi,'psi'); }
    if(!axisA||!axisB){ if(cb) cb(); return; }
    animateRotation(axisA,axisB,deltaDeg,predicate,cb);
  }
}

function animateRotation(axisA,axisB,totalDeg,predicate,cb){
  const steps=10;
  const step=totalDeg/steps;
  let i=0;
  const aA={x:axisA.x,y:axisA.y,z:axisA.z}, aB={x:axisB.x,y:axisB.y,z:axisB.z};
  function stepFn(){
    if(i<steps){
      rotateAtomsAboutAxis(aA,aB,step,predicate);
      // highlight cylinder needs to follow moved atoms – remove old and re-add next?
      // keep highlight updated after each small step for visual feedback
      i++;
      setTimeout(stepFn,35);
    } else {
      if(cb) cb();
    }
  }
  stepFn();
}

// ---------- Controls ----------
function doAlanine(checked){
  alanineShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
  alanineShapes=[];
  if(checked){
    const alaAtoms=atomsData.filter(a=>a.resi===16);
    alaAtoms.forEach(a=>{
      const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:0.9, color:'black', opacity:0.18});
      alanineShapes.push(s);
    });
    viewer.render();
  } else viewer.render();
}
function doPeptideBonds(checked){
  showPeptideBondsFlag=checked;
  if(window._pepShapes){ window._pepShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} }); window._pepShapes=[]; }
  else window._pepShapes=[];
  if(checked){
    const pairs=[[findAtom(15,'C'),findAtom(16,'N')],[findAtom(16,'C'),findAtom(17,'N')],[findAtom(15,'C'),findAtom(15,'O')],[findAtom(16,'C'),findAtom(16,'O')]];
    pairs.forEach(pair=>{
      const a=pair[0], b=pair[1];
      if(a&&b){
        const s=viewer.addCylinder({start:{x:a.x,y:a.y,z:a.z}, end:{x:b.x,y:b.y,z:b.z}, radius:0.14, color:'#ff80ff', fromCap:1,toCap:1});
        window._pepShapes.push(s);
      }
    });
    viewer.render();
  } else viewer.render();
}
function clearPlanes(){
  planeShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
  planeShapes=[];
  if(window._customPlanes){ window._customPlanes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} }); window._customPlanes=[]; }
}
function updatePlanes(){
  if(!showPlanes) return;
  clearPlanes();
  const p1a=findAtom(15,'CA'), p1b=findAtom(16,'H'), p1c=findAtom(16,'CA'), p1d=findAtom(15,'O');
  const p2a=findAtom(16,'CA'), p2b=findAtom(17,'H'), p2c=findAtom(17,'CA'), p2d=findAtom(16,'O');
  try{
    if(p1a&&p1b&&p1c&&p1d){
      const verts1=[new $3Dmol.Vector3(p1a.x,p1a.y,p1a.z), new $3Dmol.Vector3(p1b.x,p1b.y,p1b.z), new $3Dmol.Vector3(p1c.x,p1c.y,p1c.z), new $3Dmol.Vector3(p1d.x,p1d.y,p1d.z)];
      const nA=normalize(cross(vecSub(p1b,p1a), vecSub(p1c,p1a)));
      const normals1=[new $3Dmol.Vector3(nA.x,nA.y,nA.z),new $3Dmol.Vector3(nA.x,nA.y,nA.z),new $3Dmol.Vector3(nA.x,nA.y,nA.z),new $3Dmol.Vector3(nA.x,nA.y,nA.z)];
      const s1=viewer.addCustom({vertexArr:verts1, normalArr:normals1, faceArr:[0,1,2,0,2,3], color:{r:0.6,g:0.6,b:1}});
      window._customPlanes=window._customPlanes||[]; window._customPlanes.push(s1);
    }
    if(p2a&&p2b&&p2c&&p2d){
      const verts2=[new $3Dmol.Vector3(p2a.x,p2a.y,p2a.z), new $3Dmol.Vector3(p2b.x,p2b.y,p2b.z), new $3Dmol.Vector3(p2c.x,p2c.y,p2c.z), new $3Dmol.Vector3(p2d.x,p2d.y,p2d.z)];
      const nB=normalize(cross(vecSub(p2b,p2a), vecSub(p2c,p2a)));
      const normals2=[new $3Dmol.Vector3(nB.x,nB.y,nB.z),new $3Dmol.Vector3(nB.x,nB.y,nB.z),new $3Dmol.Vector3(nB.x,nB.y,nB.z),new $3Dmol.Vector3(nB.x,nB.y,nB.z)];
      const s2=viewer.addCustom({vertexArr:verts2, normalArr:normals2, faceArr:[0,1,2,0,2,3], color:{r:1,g:0.7,b:0.4}});
      window._customPlanes=window._customPlanes||[]; window._customPlanes.push(s2);
    }
  }catch(e){}
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
  showPlanes=checked;
  if(checked) updatePlanes(); else { clearPlanes(); viewer.render(); }
}
function doVDW(checked){
  showVDW=checked;
  const whiteBox=document.getElementById('divwhite');
  if(checked){
    if(whiteBox) whiteBox.style.display='';
    vdwShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    vdwShapes=[];
    const radii={H:1.20, C:1.70, N:1.55, O:1.52};
    atomsData.forEach(a=>{
      const ra=radii[a.elem]||1.5;
      let col='#c8c8c8';
      if(a.atom==='CA') col='#000000';
      else if(a.elem==='N') col='#3050ff';
      else if(a.elem==='O') col='#ff2020';
      else if(a.elem==='H') col='white';
      else if(a.elem==='C') col='#c8c8c8';
      const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:ra*0.88, color:col, opacity:0.38});
      vdwShapes.push(s);
    });
    if(showClashes) updateClashes();
    viewer.render();
  } else {
    if(whiteBox){ whiteBox.style.display='none'; document.getElementById('idwhite').checked=false; }
    vdwShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    vdwShapes=[];
    viewer.setBackgroundColor(0xd0d0d0);
    viewer.render();
    if(showClashes) updateClashes();
  }
}
function doWhite(checked){
  if(checked){
    viewer.setBackgroundColor(0xffffff);
    // make VdW overlay white and more translucent, as original does for model 1.2
    vdwShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    vdwShapes=[];
    const radii={H:1.20, C:1.70, N:1.55, O:1.52};
    atomsData.forEach(a=>{
      const ra=radii[a.elem]||1.5;
      const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:ra*0.88, color:'white', opacity:0.28});
      vdwShapes.push(s);
    });
    viewer.render();
  } else {
    viewer.setBackgroundColor(0xd0d0d0);
    // restore VdW with element colors if still checked, else remove
    vdwShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    vdwShapes=[];
    if(showVDW) doVDW(true);
    else viewer.render();
    if(showClashes) updateClashes();
  }
}
function clearClashes(){
  if(!trailClashes){
    clashShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    clashShapes=[];
  }
}
function makeClashPair(a,b, overlaps){
  // overlaps is optional Map to collect per-atom max overlap for alternate viz; we now use separate spheres for overlap only
  // This version just records overlap for later small-sphere creation; actual viz is small red spheres at midpoint
  if(!a||!b) return;
  const radii={H:1.20, C:1.70, N:1.55, O:1.52, S:1.80};
  const scale=0.88;
  const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
  const d=Math.hypot(dx,dy,dz);
  if(d>4.5) return;
  if(d<1.65){
    const bondedPairs=[
      ['15','CA','15','C'],['15','C','15','O'],['15','C','16','N'],
      ['16','N','16','CA'],['16','CA','16','C'],['16','C','16','O'],['16','CA','16','CB'],['16','CA','16','HA'],
      ['16','CB','16','1HB'],['16','CB','16','2HB'],['16','CB','16','3HB'],['16','N','16','H'],
      ['16','C','17','N'],['17','N','17','CA'],['17','N','17','H']
    ];
    for(const p of bondedPairs){
      if((a.resi==p[0]&&a.name==p[1]&&b.resi==p[2]&&b.name==p[3]) || (b.resi==p[0]&&b.name==p[1]&&a.resi==p[2]&&a.name==p[3])) return;
    }
    if(d<1.1) return;
  }
  const ra=radii[a.elem]||1.5, rb=radii[b.elem]||1.5;
  const sum=(ra+rb)*scale;
  if(d < sum && d > 0.01){
    const raa=ra*scale, rbb=rb*scale;
    // Intersection circle of the two VdW spheres: distance t from center a
    const t=(raa*raa - rbb*rbb + d*d)/(2*d);
    if(t < 0 || t > d) return; // fully contained – no lens
    const lensR=Math.sqrt(Math.max(0, raa*raa - t*t));
    if(!isFinite(lensR) || lensR<=0.01) return;
    // Record clashing atoms so their VdW spheres are skipped (avoid depth occlusion)
    clashAtomSerials.add(a.serial); clashAtomSerials.add(b.serial);
    // axis unit u from a to b, plus perpendicular basis v,w
    const ux=(b.x-a.x)/d, uy=(b.y-a.y)/d, uz=(b.z-a.z)/d;
    let vx,vy,vz;
    if(Math.abs(ux) < 0.9){ vx=1; vy=0; vz=0; }
    else { vx=0; vy=1; vz=0; }
    const dotvu=vx*ux+vy*uy+vz*uz;
    vx-=dotvu*ux; vy-=dotvu*uy; vz-=dotvu*uz;
    const vl=Math.hypot(vx,vy,vz)||1;
    vx/=vl; vy/=vl; vz/=vl;
    const wx=uy*vz-uz*vy, wy=uz*vx-ux*vz, wz=ux*vy-uy*vx;

    const segs=18, rings=5;
    const verts=[], normals=[], faces=[];
    // ---------- Cap of sphere a (points on A's surface inside B) ----------
    const thMaxA=Math.acos(Math.min(1, Math.max(-1, t/raa)));
    const ringA=[];
    verts.push({x:a.x+raa*ux, y:a.y+raa*uy, z:a.z+raa*uz});
    normals.push({x:ux, y:uy, z:uz});
    const poleA=verts.length-1;
    for(let k=1;k<=rings;k++){
      const th=thMaxA*k/rings, st=Math.sin(th), ct=Math.cos(th);
      const ring=[];
      for(let i=0;i<segs;i++){
        const ph=i/segs*2*Math.PI, cp=Math.cos(ph), sp=Math.sin(ph);
        const lx=st*cp, ly=st*sp;
        const rdx=lx*vx+ly*wx+ct*ux, rdy=lx*vy+ly*wy+ct*uy, rdz=lx*vz+ly*wz+ct*uz;
        verts.push({x:a.x+raa*rdx, y:a.y+raa*rdy, z:a.z+raa*rdz});
        normals.push({x:rdx, y:rdy, z:rdz});
        ring.push(verts.length-1);
      }
      ringA.push(ring);
    }
    for(let i=0;i<segs;i++){
      faces.push(poleA, ringA[0][(i+1)%segs], ringA[0][i]);
    }
    for(let k=0;k<rings-1;k++){
      for(let i=0;i<segs;i++){
        const i2=(i+1)%segs;
        faces.push(ringA[k][i], ringA[k][i2], ringA[k+1][i2]);
        faces.push(ringA[k][i], ringA[k+1][i2], ringA[k+1][i]);
      }
    }
    // ---------- Cap of sphere b (points on B's surface inside A) ----------
    const thMaxB=Math.acos(Math.min(1, Math.max(-1, (d-t)/rbb)));
    const ringB=[];
    verts.push({x:b.x-rbb*ux, y:b.y-rbb*uy, z:b.z-rbb*uz});
    normals.push({x:-ux, y:-uy, z:-uz});
    const poleB=verts.length-1;
    for(let k=1;k<=rings;k++){
      const th=thMaxB*k/rings, st=Math.sin(th), ct=Math.cos(th);
      const ring=[];
      for(let i=0;i<segs;i++){
        const ph=i/segs*2*Math.PI, cp=Math.cos(ph), sp=Math.sin(ph);
        const lx=st*cp, ly=st*sp;
        const rdx=lx*vx+ly*wx-ct*ux, rdy=lx*vy+ly*wy-ct*uy, rdz=lx*vz+ly*wz-ct*uz;
        verts.push({x:b.x+rbb*rdx, y:b.y+rbb*rdy, z:b.z+rbb*rdz});
        normals.push({x:rdx, y:rdy, z:rdz});
        ring.push(verts.length-1);
      }
      ringB.push(ring);
    }
    for(let i=0;i<segs;i++){
      faces.push(poleB, ringB[0][(i+1)%segs], ringB[0][i]);
    }
    for(let k=0;k<rings-1;k++){
      for(let i=0;i<segs;i++){
        const i2=(i+1)%segs;
        faces.push(ringB[k][i], ringB[k][i2], ringB[k+1][i2]);
        faces.push(ringB[k][i], ringB[k+1][i2], ringB[k+1][i]);
      }
    }
    // Lens mesh: two spherical caps sharing the exact intersection rim
    const spec={vertexArr:verts, normalArr:normals, faceArr:faces, color:0xffa500, opacity:0.85, side:2};
    try{
      const shape=viewer.addCustom(spec);
      clashShapes.push(shape);
    }catch(e){
      const s=viewer.addSphere({center:{x:a.x+(b.x-a.x)*t/d, y:a.y+(b.y-a.y)*t/d, z:a.z+(b.z-a.z)*t/d}, radius:lensR*0.8, color:'orange', opacity:0.85});
      clashShapes.push(s);
    }
  }
}
function updateClashes(){
  if(!showClashes) return;
  if(!trailClashes){
    clashShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    clashShapes=[];
  }
  clashAtomSerials = new Set();
  // Orange lens at overlap only - do NOT auto-show VdW; keep VdW as is (if VdW off, only lenses show)
  const sidechain=['CB','1HB','2HB','3HB'].map(n=>findAtom(16,n)).filter(Boolean);
  const ca1=[findAtom(15,'O'),findAtom(16,'H'),findAtom(16,'O'),findAtom(17,'H')].filter(Boolean);
  sidechain.forEach(sc=> ca1.forEach(c=> makeClashPair(sc,c)));
  const o15=findAtom(15,'O');
  const h17=findAtom(17,'H'), n17=findAtom(17,'N');
  if(o15&&h17) makeClashPair(o15,h17);
  if(o15&&n17) makeClashPair(o15,n17);
  const o16=findAtom(16,'O');
  if(o15&&o16) makeClashPair(o15,o16);
  const h16=findAtom(16,'H');
  if(h16&&h17) makeClashPair(h16,h17);
  const ha=findAtom(16,'HA');
  if(ha&&o15) makeClashPair(ha,o15);
  if(ha&&h17) makeClashPair(ha,h17);
  viewer.render();
}
function doClashes(checked){
  showClashes=checked;
  const trailDiv=document.getElementById('divtrailclashes');
  if(checked){
    if(trailDiv) trailDiv.style.display='';
    updateClashes();
  } else {
    clashShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    clashShapes=[];
    if(trailDiv){ trailDiv.style.display='none'; document.getElementById('idtrailclashes').checked=false; trailClashes=false; }
    viewer.render();
  }
}
function doTrailClashes(checked){
  trailClashes=checked;
  if(!checked){
    clashShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    clashShapes=[];
    if(showClashes) updateClashes(); else viewer.render();
  }
}
function resetViewer(){
  document.getElementById('idalanine').checked=false;
  document.getElementById('idpeptidebonds').checked=false;
  document.getElementById('idplanes').checked=false;
  document.getElementById('idvdw').checked=false;
  document.getElementById('idwhite').checked=false;
  document.getElementById('idclashes').checked=false;
  document.getElementById('idtrailclashes').checked=false;
  document.getElementById('divwhite').style.display='none';
  document.getElementById('divtrailclashes').style.display='none';
  showPlanes=false; showVDW=false; showClashes=false; trailClashes=false;
  clearPlanes();
  clashShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} }); clashShapes=[];
  alanineShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} }); alanineShapes=[];
  if(window._pepShapes){ window._pepShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} }); window._pepShapes=[]; }
  if(window._highlightShapes){ window._highlightShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} }); window._highlightShapes=[]; }
  if(window._customPlanes){ window._customPlanes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} }); window._customPlanes=[]; }
  loadTripeptide();
}

// Ramachandran heatplot — real backbone phi/psi density
// Grid: 30,394 residues from 60 diverse PDB structures (non-Gly/Pro), 5 deg cells.
// RAMA_DENSITY[row][col] = density at (phi=RAMA_PHI[col], psi=RAMA_PHI[row]).
const RAMA_N=72;
const RAMA_PHI=[]; for(let i=0;i<RAMA_N;i++) RAMA_PHI.push(-177.5+i*5);
const RAMA_DENSITY=[[0.12,0.132,0.143,0.145,0.141,0.15,0.17,0.178,0.176,0.17,0.159,0.145,0.133,0.122,0.115,0.116,0.126,0.135,0.138,0.137,0.135,0.131,0.121,0.106,0.093,0.079,0.06,0.038,0.022,0.023,0.025,0.022,0.014,0.006,0.002,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.009,0.018,0.028,0.034,0.036,0.043,0.048,0.042,0.031,0.021,0.012,0.006,0.002,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.004,0.011,0.026,0.044,0.059,0.067,0.071,0.073],[0.112,0.124,0.135,0.135,0.129,0.135,0.151,0.161,0.162,0.16,0.149,0.131,0.115,0.108,0.108,0.113,0.123,0.131,0.131,0.125,0.121,0.118,0.11,0.097,0.083,0.07,0.053,0.034,0.022,0.023,0.026,0.023,0.014,0.006,0.002,0.001,0.0,0.0,0.0,0.0,0.001,0.003,0.007,0.012,0.017,0.024,0.027,0.031,0.039,0.045,0.044,0.038,0.03,0.02,0.01,0.004,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.009,0.021,0.037,0.053,0.063,0.068,0.071],[0.099,0.107,0.114,0.115,0.111,0.114,0.126,0.137,0.143,0.145,0.136,0.114,0.093,0.089,0.096,0.104,0.111,0.115,0.114,0.106,0.1,0.099,0.095,0.083,0.07,0.058,0.044,0.028,0.018,0.02,0.022,0.019,0.012,0.005,0.002,0.0,0.0,0.0,0.001,0.001,0.003,0.007,0.014,0.022,0.027,0.025,0.021,0.023,0.03,0.038,0.041,0.042,0.037,0.026,0.014,0.006,0.002,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.002,0.006,0.014,0.026,0.04,0.053,0.062,0.068],[0.079,0.082,0.085,0.088,0.089,0.092,0.101,0.112,0.122,0.128,0.122,0.101,0.078,0.074,0.083,0.089,0.092,0.094,0.092,0.087,0.082,0.08,0.076,0.066,0.058,0.051,0.039,0.025,0.015,0.013,0.014,0.012,0.007,0.003,0.001,0.0,0.001,0.001,0.002,0.003,0.006,0.012,0.023,0.036,0.042,0.037,0.027,0.022,0.024,0.027,0.033,0.038,0.036,0.027,0.016,0.008,0.003,0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.007,0.016,0.027,0.042,0.057,0.067],[0.054,0.054,0.056,0.062,0.068,0.072,0.077,0.086,0.099,0.109,0.108,0.091,0.07,0.064,0.069,0.073,0.074,0.075,0.075,0.072,0.067,0.061,0.056,0.051,0.049,0.047,0.039,0.026,0.015,0.008,0.007,0.005,0.003,0.002,0.001,0.001,0.002,0.004,0.006,0.009,0.013,0.02,0.031,0.044,0.05,0.045,0.036,0.034,0.032,0.026,0.024,0.027,0.029,0.026,0.019,0.013,0.008,0.004,0.002,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.009,0.018,0.033,0.05,0.061],[0.03,0.03,0.034,0.044,0.053,0.057,0.06,0.067,0.08,0.092,0.094,0.082,0.065,0.057,0.057,0.058,0.059,0.059,0.059,0.058,0.053,0.047,0.043,0.044,0.045,0.043,0.036,0.026,0.015,0.007,0.003,0.002,0.001,0.001,0.0,0.001,0.004,0.008,0.014,0.019,0.023,0.03,0.037,0.044,0.047,0.045,0.044,0.047,0.046,0.036,0.023,0.021,0.028,0.032,0.031,0.025,0.016,0.009,0.004,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.002,0.005,0.011,0.023,0.036,0.046],[0.013,0.015,0.023,0.035,0.043,0.045,0.047,0.053,0.066,0.079,0.082,0.075,0.063,0.055,0.055,0.057,0.055,0.051,0.047,0.045,0.042,0.039,0.04,0.043,0.042,0.037,0.03,0.022,0.013,0.006,0.003,0.002,0.001,0.001,0.0,0.002,0.006,0.013,0.022,0.029,0.033,0.037,0.04,0.039,0.039,0.042,0.048,0.055,0.053,0.041,0.026,0.022,0.033,0.043,0.045,0.038,0.026,0.014,0.006,0.002,0.0,0.0,0.0,0.0,0.0,0.001,0.001,0.003,0.006,0.013,0.021,0.027],[0.009,0.01,0.017,0.028,0.036,0.039,0.04,0.043,0.053,0.065,0.07,0.066,0.059,0.056,0.059,0.063,0.061,0.053,0.043,0.037,0.034,0.033,0.036,0.04,0.039,0.034,0.028,0.021,0.013,0.008,0.007,0.005,0.003,0.002,0.001,0.002,0.007,0.015,0.025,0.034,0.036,0.036,0.037,0.036,0.036,0.044,0.054,0.058,0.052,0.038,0.023,0.022,0.035,0.048,0.052,0.045,0.031,0.017,0.007,0.002,0.001,0.0,0.0,0.0,0.001,0.002,0.004,0.006,0.007,0.008,0.011,0.014],[0.018,0.014,0.014,0.021,0.031,0.038,0.041,0.041,0.046,0.054,0.058,0.057,0.056,0.056,0.059,0.063,0.062,0.054,0.044,0.038,0.037,0.034,0.034,0.037,0.038,0.037,0.034,0.027,0.018,0.014,0.014,0.012,0.007,0.003,0.001,0.002,0.006,0.013,0.022,0.029,0.031,0.031,0.035,0.039,0.045,0.057,0.068,0.067,0.055,0.036,0.023,0.022,0.031,0.042,0.046,0.04,0.028,0.015,0.006,0.002,0.0,0.0,0.0,0.0,0.001,0.004,0.008,0.013,0.015,0.015,0.016,0.018],[0.028,0.022,0.015,0.018,0.029,0.039,0.042,0.041,0.045,0.05,0.052,0.052,0.054,0.053,0.052,0.054,0.054,0.048,0.043,0.043,0.046,0.043,0.037,0.036,0.036,0.038,0.039,0.032,0.023,0.021,0.022,0.019,0.012,0.006,0.003,0.003,0.004,0.009,0.016,0.022,0.024,0.027,0.036,0.047,0.057,0.07,0.081,0.079,0.063,0.043,0.03,0.028,0.028,0.031,0.032,0.028,0.019,0.01,0.004,0.001,0.0,0.0,0.0,0.001,0.002,0.006,0.014,0.022,0.026,0.025,0.025,0.029],[0.033,0.025,0.016,0.017,0.027,0.035,0.036,0.036,0.042,0.048,0.049,0.047,0.047,0.046,0.044,0.043,0.043,0.041,0.041,0.048,0.052,0.048,0.038,0.033,0.034,0.038,0.04,0.035,0.027,0.026,0.027,0.022,0.014,0.009,0.008,0.008,0.007,0.008,0.016,0.025,0.03,0.03,0.035,0.048,0.061,0.075,0.084,0.082,0.069,0.05,0.038,0.035,0.03,0.022,0.019,0.015,0.01,0.005,0.002,0.001,0.0,0.0,0.0,0.001,0.003,0.008,0.018,0.028,0.033,0.031,0.03,0.034],[0.028,0.022,0.015,0.018,0.029,0.035,0.033,0.03,0.037,0.047,0.049,0.045,0.04,0.04,0.042,0.041,0.041,0.041,0.043,0.048,0.05,0.045,0.035,0.033,0.039,0.046,0.048,0.043,0.033,0.028,0.025,0.02,0.014,0.014,0.018,0.018,0.014,0.012,0.021,0.033,0.039,0.035,0.032,0.042,0.057,0.07,0.077,0.077,0.069,0.052,0.039,0.035,0.029,0.018,0.01,0.007,0.004,0.002,0.001,0.0,0.0,0.0,0.0,0.001,0.003,0.008,0.018,0.028,0.033,0.031,0.028,0.03],[0.018,0.014,0.012,0.02,0.032,0.039,0.037,0.031,0.036,0.048,0.054,0.051,0.045,0.044,0.045,0.046,0.046,0.045,0.044,0.044,0.043,0.038,0.035,0.042,0.052,0.057,0.058,0.05,0.038,0.028,0.02,0.014,0.015,0.022,0.028,0.028,0.022,0.016,0.023,0.035,0.041,0.036,0.028,0.037,0.052,0.064,0.068,0.067,0.06,0.046,0.033,0.027,0.022,0.014,0.007,0.003,0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.002,0.006,0.014,0.022,0.026,0.026,0.023,0.023],[0.008,0.007,0.009,0.019,0.032,0.042,0.042,0.035,0.034,0.045,0.054,0.056,0.054,0.053,0.05,0.048,0.048,0.047,0.044,0.039,0.037,0.037,0.042,0.053,0.062,0.066,0.062,0.051,0.039,0.027,0.018,0.013,0.016,0.026,0.033,0.033,0.026,0.018,0.022,0.031,0.035,0.03,0.026,0.036,0.051,0.059,0.06,0.055,0.046,0.033,0.022,0.017,0.013,0.008,0.004,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.004,0.008,0.014,0.02,0.026,0.028,0.027],[0.003,0.005,0.009,0.019,0.033,0.045,0.048,0.039,0.031,0.037,0.047,0.053,0.056,0.056,0.051,0.046,0.047,0.05,0.049,0.044,0.042,0.048,0.055,0.062,0.07,0.072,0.066,0.054,0.043,0.035,0.028,0.02,0.02,0.026,0.03,0.029,0.023,0.019,0.024,0.031,0.03,0.025,0.025,0.035,0.047,0.053,0.053,0.047,0.037,0.025,0.014,0.008,0.006,0.004,0.002,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.002,0.006,0.012,0.023,0.033,0.038,0.035],[0.004,0.008,0.014,0.024,0.036,0.048,0.05,0.041,0.03,0.029,0.037,0.046,0.051,0.051,0.045,0.043,0.05,0.057,0.058,0.055,0.058,0.066,0.07,0.072,0.078,0.082,0.077,0.066,0.056,0.052,0.044,0.033,0.028,0.029,0.027,0.022,0.016,0.02,0.031,0.038,0.038,0.032,0.032,0.039,0.045,0.049,0.052,0.049,0.04,0.027,0.015,0.006,0.003,0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.001,0.002,0.003,0.007,0.016,0.03,0.042,0.045,0.041],[0.006,0.013,0.021,0.03,0.038,0.047,0.05,0.046,0.038,0.033,0.037,0.045,0.049,0.046,0.043,0.048,0.057,0.064,0.065,0.064,0.07,0.081,0.084,0.084,0.089,0.095,0.092,0.082,0.071,0.066,0.057,0.044,0.036,0.034,0.029,0.02,0.014,0.02,0.034,0.045,0.048,0.043,0.041,0.044,0.046,0.048,0.051,0.05,0.042,0.029,0.016,0.007,0.004,0.003,0.002,0.001,0.001,0.0,0.0,0.0,0.001,0.001,0.002,0.004,0.006,0.007,0.01,0.02,0.035,0.046,0.049,0.043],[0.008,0.017,0.028,0.036,0.04,0.045,0.053,0.055,0.049,0.042,0.044,0.052,0.057,0.054,0.053,0.059,0.066,0.069,0.069,0.07,0.08,0.094,0.101,0.103,0.108,0.114,0.113,0.101,0.085,0.074,0.065,0.051,0.04,0.036,0.03,0.021,0.017,0.02,0.031,0.045,0.051,0.048,0.043,0.044,0.046,0.046,0.049,0.046,0.037,0.025,0.015,0.011,0.011,0.01,0.007,0.004,0.002,0.001,0.0,0.001,0.002,0.004,0.006,0.01,0.013,0.014,0.014,0.021,0.034,0.045,0.046,0.04],[0.012,0.022,0.035,0.043,0.044,0.046,0.057,0.063,0.059,0.049,0.049,0.06,0.069,0.069,0.07,0.075,0.079,0.082,0.083,0.085,0.096,0.112,0.124,0.13,0.136,0.144,0.145,0.131,0.107,0.086,0.073,0.059,0.042,0.033,0.028,0.026,0.024,0.022,0.025,0.036,0.042,0.04,0.037,0.04,0.044,0.047,0.049,0.045,0.034,0.022,0.019,0.023,0.025,0.022,0.016,0.009,0.004,0.001,0.0,0.001,0.003,0.008,0.013,0.018,0.022,0.024,0.02,0.02,0.029,0.039,0.04,0.036],[0.021,0.027,0.04,0.049,0.048,0.048,0.059,0.068,0.066,0.055,0.054,0.068,0.079,0.083,0.086,0.092,0.098,0.101,0.104,0.108,0.119,0.137,0.154,0.165,0.175,0.185,0.187,0.171,0.139,0.107,0.085,0.067,0.047,0.034,0.031,0.032,0.03,0.025,0.021,0.026,0.031,0.031,0.031,0.038,0.047,0.054,0.058,0.053,0.04,0.027,0.03,0.04,0.043,0.037,0.026,0.014,0.006,0.002,0.001,0.002,0.006,0.012,0.021,0.027,0.031,0.03,0.024,0.02,0.027,0.036,0.038,0.036],[0.032,0.034,0.042,0.049,0.047,0.047,0.059,0.07,0.07,0.062,0.061,0.072,0.083,0.089,0.096,0.105,0.111,0.115,0.123,0.132,0.147,0.172,0.198,0.216,0.229,0.238,0.236,0.217,0.181,0.14,0.106,0.08,0.057,0.042,0.038,0.036,0.032,0.026,0.026,0.031,0.035,0.035,0.035,0.041,0.053,0.063,0.068,0.063,0.049,0.035,0.04,0.052,0.055,0.047,0.032,0.017,0.007,0.002,0.001,0.002,0.007,0.015,0.025,0.033,0.035,0.031,0.023,0.019,0.026,0.035,0.04,0.043],[0.042,0.04,0.044,0.048,0.047,0.048,0.06,0.071,0.074,0.071,0.071,0.075,0.081,0.089,0.099,0.108,0.115,0.124,0.138,0.156,0.183,0.222,0.262,0.292,0.311,0.316,0.307,0.279,0.236,0.186,0.14,0.103,0.072,0.05,0.041,0.038,0.032,0.029,0.037,0.045,0.046,0.042,0.039,0.043,0.056,0.067,0.072,0.067,0.052,0.038,0.043,0.054,0.056,0.046,0.03,0.016,0.006,0.002,0.001,0.002,0.006,0.013,0.023,0.031,0.033,0.028,0.019,0.015,0.022,0.032,0.043,0.051],[0.049,0.046,0.047,0.052,0.054,0.056,0.066,0.077,0.081,0.082,0.081,0.081,0.084,0.094,0.105,0.115,0.124,0.14,0.163,0.19,0.231,0.287,0.344,0.39,0.416,0.418,0.397,0.355,0.298,0.236,0.177,0.125,0.083,0.055,0.043,0.04,0.036,0.035,0.044,0.051,0.048,0.042,0.038,0.039,0.05,0.06,0.064,0.059,0.046,0.034,0.037,0.046,0.046,0.037,0.023,0.012,0.004,0.001,0.0,0.001,0.004,0.009,0.017,0.025,0.028,0.024,0.016,0.011,0.016,0.027,0.042,0.053],[0.053,0.052,0.057,0.064,0.067,0.068,0.075,0.085,0.088,0.089,0.091,0.094,0.1,0.11,0.122,0.133,0.147,0.167,0.196,0.235,0.293,0.367,0.443,0.502,0.531,0.524,0.487,0.427,0.353,0.277,0.205,0.141,0.091,0.061,0.048,0.046,0.044,0.041,0.042,0.044,0.04,0.035,0.032,0.032,0.037,0.045,0.047,0.044,0.037,0.033,0.037,0.038,0.034,0.026,0.016,0.007,0.003,0.001,0.0,0.001,0.003,0.007,0.014,0.023,0.027,0.024,0.019,0.014,0.015,0.025,0.037,0.046],[0.05,0.057,0.069,0.078,0.079,0.077,0.081,0.087,0.089,0.092,0.101,0.111,0.121,0.132,0.142,0.157,0.175,0.197,0.231,0.285,0.367,0.47,0.569,0.639,0.659,0.63,0.567,0.484,0.393,0.306,0.224,0.152,0.1,0.068,0.052,0.051,0.051,0.045,0.036,0.03,0.026,0.024,0.024,0.022,0.023,0.028,0.032,0.034,0.036,0.043,0.047,0.042,0.031,0.02,0.012,0.006,0.003,0.001,0.001,0.002,0.004,0.008,0.015,0.024,0.029,0.029,0.026,0.022,0.024,0.031,0.039,0.043],[0.047,0.059,0.075,0.087,0.088,0.086,0.087,0.089,0.089,0.094,0.108,0.124,0.137,0.149,0.162,0.18,0.202,0.226,0.267,0.34,0.453,0.592,0.722,0.8,0.8,0.734,0.633,0.522,0.414,0.317,0.23,0.157,0.105,0.072,0.052,0.047,0.048,0.042,0.03,0.019,0.014,0.013,0.014,0.013,0.012,0.016,0.023,0.03,0.042,0.056,0.061,0.054,0.041,0.029,0.02,0.011,0.005,0.002,0.001,0.004,0.008,0.014,0.019,0.023,0.027,0.029,0.029,0.029,0.034,0.043,0.048,0.05],[0.046,0.056,0.071,0.085,0.091,0.094,0.097,0.096,0.093,0.097,0.11,0.128,0.145,0.161,0.179,0.202,0.227,0.257,0.305,0.395,0.534,0.705,0.861,0.942,0.916,0.809,0.669,0.533,0.411,0.305,0.216,0.147,0.099,0.067,0.046,0.037,0.037,0.031,0.021,0.013,0.008,0.007,0.007,0.006,0.007,0.011,0.018,0.03,0.047,0.064,0.072,0.066,0.055,0.044,0.032,0.019,0.009,0.003,0.003,0.007,0.014,0.023,0.028,0.027,0.024,0.024,0.026,0.032,0.042,0.052,0.056,0.056],[0.045,0.049,0.06,0.075,0.087,0.096,0.102,0.101,0.098,0.101,0.113,0.131,0.15,0.168,0.19,0.22,0.251,0.288,0.343,0.439,0.588,0.769,0.929,1.0,0.952,0.818,0.658,0.512,0.386,0.277,0.187,0.123,0.083,0.057,0.04,0.029,0.024,0.019,0.016,0.015,0.012,0.008,0.004,0.005,0.01,0.017,0.024,0.033,0.05,0.068,0.077,0.075,0.066,0.055,0.041,0.026,0.013,0.006,0.006,0.01,0.019,0.03,0.037,0.035,0.026,0.02,0.024,0.035,0.048,0.056,0.06,0.061],[0.041,0.042,0.048,0.062,0.078,0.092,0.101,0.101,0.1,0.106,0.12,0.138,0.157,0.175,0.199,0.23,0.267,0.31,0.369,0.464,0.604,0.767,0.903,0.954,0.894,0.757,0.602,0.462,0.343,0.24,0.158,0.102,0.068,0.05,0.039,0.03,0.019,0.015,0.02,0.022,0.019,0.012,0.006,0.008,0.016,0.027,0.035,0.042,0.055,0.071,0.08,0.079,0.07,0.058,0.045,0.031,0.019,0.012,0.012,0.016,0.021,0.031,0.039,0.037,0.028,0.019,0.024,0.038,0.05,0.054,0.056,0.059],[0.046,0.046,0.048,0.055,0.069,0.087,0.098,0.1,0.102,0.114,0.132,0.15,0.166,0.183,0.207,0.237,0.273,0.317,0.376,0.463,0.582,0.711,0.811,0.838,0.777,0.655,0.519,0.396,0.291,0.203,0.137,0.093,0.065,0.049,0.041,0.032,0.02,0.016,0.022,0.026,0.022,0.014,0.007,0.009,0.019,0.031,0.04,0.047,0.058,0.071,0.079,0.078,0.068,0.056,0.044,0.034,0.024,0.017,0.02,0.023,0.023,0.027,0.031,0.03,0.023,0.016,0.022,0.035,0.044,0.045,0.046,0.049],[0.062,0.062,0.059,0.057,0.066,0.084,0.096,0.099,0.106,0.124,0.144,0.161,0.175,0.191,0.212,0.238,0.271,0.316,0.373,0.448,0.542,0.635,0.698,0.705,0.646,0.545,0.431,0.326,0.236,0.164,0.115,0.084,0.063,0.047,0.037,0.027,0.017,0.014,0.019,0.022,0.019,0.012,0.006,0.008,0.017,0.028,0.036,0.042,0.052,0.064,0.073,0.073,0.067,0.055,0.043,0.034,0.026,0.02,0.023,0.027,0.026,0.024,0.022,0.02,0.014,0.01,0.016,0.025,0.031,0.031,0.031,0.035],[0.077,0.077,0.071,0.063,0.068,0.082,0.094,0.101,0.113,0.133,0.153,0.168,0.182,0.198,0.216,0.236,0.269,0.314,0.369,0.433,0.502,0.564,0.599,0.588,0.53,0.443,0.349,0.26,0.183,0.125,0.086,0.065,0.051,0.038,0.028,0.019,0.011,0.009,0.012,0.014,0.012,0.007,0.004,0.006,0.013,0.022,0.031,0.038,0.045,0.056,0.065,0.068,0.065,0.054,0.04,0.029,0.021,0.017,0.02,0.025,0.028,0.028,0.023,0.016,0.009,0.006,0.01,0.015,0.018,0.018,0.018,0.02],[0.082,0.083,0.078,0.068,0.069,0.081,0.093,0.103,0.119,0.141,0.162,0.178,0.195,0.213,0.228,0.246,0.276,0.319,0.37,0.422,0.471,0.507,0.52,0.496,0.436,0.357,0.276,0.202,0.138,0.09,0.059,0.044,0.036,0.03,0.023,0.014,0.007,0.004,0.005,0.006,0.005,0.004,0.004,0.009,0.016,0.027,0.038,0.046,0.051,0.058,0.064,0.065,0.061,0.051,0.036,0.022,0.014,0.011,0.014,0.022,0.03,0.034,0.028,0.018,0.009,0.01,0.014,0.017,0.018,0.017,0.015,0.012],[0.073,0.078,0.076,0.069,0.07,0.082,0.095,0.107,0.123,0.146,0.172,0.194,0.214,0.233,0.248,0.266,0.294,0.333,0.375,0.413,0.443,0.458,0.452,0.418,0.359,0.287,0.216,0.153,0.102,0.064,0.042,0.031,0.028,0.028,0.023,0.014,0.007,0.002,0.002,0.002,0.002,0.002,0.006,0.013,0.024,0.037,0.05,0.06,0.066,0.071,0.072,0.067,0.058,0.045,0.03,0.017,0.008,0.006,0.01,0.018,0.029,0.033,0.028,0.018,0.012,0.018,0.025,0.029,0.029,0.028,0.022,0.015],[0.056,0.064,0.066,0.065,0.073,0.088,0.101,0.112,0.127,0.15,0.179,0.205,0.227,0.246,0.263,0.283,0.312,0.346,0.377,0.397,0.407,0.406,0.387,0.347,0.291,0.228,0.168,0.117,0.078,0.05,0.032,0.024,0.022,0.023,0.019,0.012,0.006,0.003,0.003,0.003,0.003,0.004,0.007,0.016,0.029,0.044,0.059,0.072,0.081,0.085,0.083,0.074,0.06,0.044,0.027,0.013,0.005,0.003,0.006,0.014,0.022,0.025,0.022,0.015,0.017,0.029,0.038,0.041,0.04,0.037,0.029,0.019],[0.042,0.049,0.054,0.06,0.074,0.093,0.107,0.116,0.129,0.151,0.179,0.208,0.233,0.253,0.27,0.292,0.32,0.349,0.368,0.371,0.364,0.349,0.321,0.278,0.228,0.177,0.13,0.093,0.066,0.043,0.026,0.016,0.014,0.014,0.012,0.008,0.007,0.009,0.01,0.01,0.01,0.009,0.01,0.016,0.029,0.045,0.063,0.079,0.092,0.098,0.093,0.082,0.065,0.046,0.027,0.013,0.005,0.002,0.004,0.008,0.013,0.015,0.013,0.012,0.022,0.038,0.049,0.051,0.045,0.038,0.029,0.019],[0.043,0.046,0.047,0.053,0.071,0.093,0.108,0.118,0.131,0.151,0.178,0.208,0.236,0.257,0.273,0.292,0.316,0.337,0.345,0.336,0.319,0.295,0.26,0.216,0.173,0.134,0.1,0.073,0.053,0.035,0.02,0.01,0.007,0.007,0.006,0.009,0.014,0.019,0.022,0.022,0.022,0.02,0.019,0.022,0.031,0.048,0.067,0.087,0.104,0.111,0.106,0.092,0.07,0.048,0.029,0.015,0.006,0.002,0.002,0.004,0.006,0.007,0.006,0.011,0.024,0.04,0.053,0.054,0.045,0.034,0.024,0.015],[0.056,0.055,0.051,0.05,0.064,0.085,0.103,0.115,0.131,0.155,0.183,0.211,0.236,0.258,0.273,0.287,0.302,0.312,0.31,0.296,0.274,0.246,0.208,0.168,0.134,0.105,0.079,0.057,0.041,0.027,0.014,0.006,0.003,0.003,0.006,0.013,0.022,0.03,0.035,0.036,0.035,0.033,0.031,0.033,0.042,0.058,0.078,0.1,0.119,0.127,0.121,0.103,0.076,0.05,0.03,0.016,0.007,0.002,0.001,0.001,0.002,0.002,0.003,0.008,0.019,0.034,0.045,0.047,0.038,0.026,0.016,0.009],[0.063,0.062,0.059,0.056,0.063,0.079,0.096,0.112,0.133,0.161,0.191,0.214,0.233,0.252,0.266,0.276,0.28,0.277,0.268,0.252,0.23,0.203,0.171,0.139,0.113,0.091,0.071,0.053,0.038,0.025,0.014,0.006,0.002,0.002,0.007,0.015,0.026,0.035,0.041,0.042,0.041,0.039,0.04,0.045,0.056,0.074,0.096,0.118,0.135,0.139,0.13,0.108,0.079,0.051,0.031,0.017,0.007,0.002,0.001,0.0,0.0,0.001,0.002,0.005,0.012,0.022,0.03,0.032,0.027,0.018,0.011,0.006],[0.058,0.059,0.059,0.063,0.071,0.083,0.099,0.114,0.133,0.159,0.187,0.209,0.225,0.24,0.253,0.256,0.25,0.238,0.225,0.21,0.194,0.173,0.149,0.124,0.103,0.086,0.07,0.053,0.037,0.025,0.015,0.006,0.002,0.002,0.006,0.013,0.022,0.03,0.035,0.036,0.036,0.036,0.042,0.053,0.068,0.091,0.117,0.139,0.149,0.144,0.128,0.104,0.077,0.053,0.034,0.02,0.009,0.004,0.002,0.001,0.001,0.0,0.001,0.002,0.006,0.011,0.016,0.019,0.021,0.021,0.017,0.011],[0.047,0.047,0.052,0.064,0.077,0.09,0.104,0.118,0.132,0.15,0.173,0.193,0.208,0.223,0.232,0.228,0.215,0.2,0.187,0.176,0.166,0.153,0.135,0.114,0.095,0.079,0.063,0.046,0.032,0.021,0.012,0.005,0.002,0.002,0.004,0.008,0.014,0.019,0.022,0.024,0.024,0.029,0.042,0.059,0.08,0.108,0.139,0.159,0.163,0.149,0.125,0.098,0.074,0.054,0.038,0.023,0.012,0.007,0.006,0.003,0.002,0.001,0.0,0.001,0.002,0.005,0.009,0.017,0.026,0.031,0.027,0.018],[0.037,0.037,0.046,0.063,0.078,0.09,0.103,0.116,0.131,0.145,0.159,0.174,0.187,0.199,0.205,0.197,0.182,0.17,0.16,0.15,0.143,0.135,0.121,0.104,0.086,0.068,0.05,0.034,0.022,0.014,0.008,0.003,0.002,0.003,0.006,0.007,0.009,0.012,0.017,0.021,0.022,0.027,0.044,0.069,0.098,0.132,0.162,0.179,0.177,0.158,0.129,0.1,0.075,0.055,0.04,0.027,0.018,0.015,0.012,0.007,0.003,0.001,0.0,0.0,0.001,0.003,0.009,0.019,0.031,0.037,0.033,0.024],[0.037,0.039,0.051,0.068,0.081,0.091,0.101,0.113,0.127,0.139,0.148,0.155,0.165,0.175,0.179,0.169,0.155,0.147,0.142,0.135,0.128,0.121,0.108,0.094,0.078,0.058,0.037,0.021,0.012,0.007,0.004,0.002,0.003,0.007,0.012,0.014,0.013,0.015,0.023,0.031,0.032,0.033,0.048,0.079,0.118,0.156,0.185,0.195,0.187,0.165,0.137,0.107,0.079,0.057,0.042,0.03,0.024,0.023,0.019,0.012,0.005,0.002,0.0,0.001,0.002,0.004,0.009,0.018,0.028,0.035,0.034,0.029],[0.044,0.046,0.058,0.073,0.084,0.093,0.101,0.111,0.122,0.13,0.133,0.134,0.141,0.153,0.158,0.148,0.133,0.129,0.131,0.13,0.125,0.115,0.101,0.088,0.071,0.05,0.029,0.014,0.006,0.003,0.001,0.002,0.005,0.012,0.019,0.022,0.02,0.019,0.028,0.039,0.044,0.045,0.056,0.089,0.133,0.175,0.203,0.208,0.193,0.168,0.139,0.111,0.085,0.064,0.047,0.032,0.026,0.026,0.022,0.014,0.006,0.002,0.001,0.001,0.003,0.008,0.013,0.018,0.024,0.03,0.035,0.035],[0.048,0.049,0.056,0.068,0.079,0.09,0.1,0.108,0.114,0.118,0.118,0.116,0.122,0.135,0.14,0.13,0.118,0.12,0.128,0.131,0.125,0.113,0.098,0.083,0.065,0.045,0.026,0.012,0.005,0.003,0.001,0.002,0.006,0.014,0.022,0.026,0.023,0.02,0.029,0.044,0.056,0.063,0.072,0.101,0.144,0.186,0.211,0.211,0.191,0.161,0.133,0.109,0.089,0.071,0.051,0.032,0.023,0.023,0.019,0.012,0.005,0.002,0.001,0.002,0.005,0.012,0.019,0.023,0.024,0.029,0.037,0.043],[0.048,0.047,0.05,0.056,0.066,0.079,0.095,0.106,0.11,0.111,0.111,0.107,0.111,0.122,0.126,0.118,0.114,0.122,0.132,0.134,0.125,0.112,0.097,0.082,0.065,0.046,0.029,0.017,0.011,0.007,0.004,0.002,0.005,0.012,0.019,0.022,0.02,0.018,0.03,0.05,0.071,0.083,0.091,0.111,0.147,0.184,0.203,0.198,0.174,0.144,0.118,0.1,0.086,0.07,0.05,0.029,0.017,0.014,0.012,0.007,0.003,0.001,0.001,0.002,0.006,0.014,0.023,0.028,0.028,0.031,0.04,0.047],[0.045,0.046,0.044,0.044,0.051,0.068,0.089,0.106,0.112,0.114,0.114,0.108,0.105,0.112,0.115,0.113,0.117,0.13,0.14,0.14,0.129,0.115,0.101,0.086,0.069,0.05,0.035,0.027,0.023,0.016,0.009,0.004,0.004,0.007,0.012,0.014,0.013,0.018,0.034,0.057,0.08,0.095,0.1,0.114,0.141,0.168,0.181,0.171,0.146,0.119,0.099,0.086,0.074,0.06,0.042,0.024,0.011,0.007,0.005,0.003,0.002,0.001,0.0,0.002,0.005,0.012,0.02,0.027,0.032,0.036,0.039,0.042],[0.042,0.044,0.043,0.041,0.05,0.068,0.092,0.111,0.123,0.129,0.126,0.115,0.103,0.101,0.103,0.108,0.12,0.137,0.148,0.146,0.135,0.123,0.109,0.091,0.069,0.051,0.043,0.041,0.037,0.027,0.015,0.007,0.003,0.004,0.006,0.007,0.01,0.02,0.037,0.058,0.078,0.091,0.097,0.107,0.127,0.146,0.153,0.143,0.121,0.101,0.088,0.076,0.063,0.048,0.033,0.019,0.008,0.003,0.002,0.001,0.001,0.0,0.0,0.001,0.003,0.008,0.014,0.023,0.033,0.039,0.038,0.034],[0.039,0.046,0.05,0.053,0.064,0.082,0.105,0.126,0.142,0.148,0.143,0.128,0.109,0.097,0.096,0.106,0.124,0.143,0.153,0.148,0.137,0.127,0.115,0.093,0.069,0.053,0.051,0.053,0.048,0.035,0.02,0.01,0.006,0.007,0.006,0.005,0.01,0.021,0.037,0.054,0.069,0.079,0.085,0.093,0.108,0.123,0.129,0.122,0.108,0.096,0.088,0.076,0.06,0.043,0.027,0.015,0.006,0.002,0.001,0.0,0.0,0.0,0.0,0.001,0.002,0.004,0.01,0.019,0.031,0.04,0.039,0.034],[0.037,0.05,0.061,0.067,0.077,0.095,0.117,0.14,0.157,0.164,0.158,0.143,0.123,0.107,0.104,0.118,0.138,0.154,0.159,0.149,0.134,0.124,0.112,0.09,0.069,0.06,0.058,0.056,0.049,0.037,0.023,0.014,0.013,0.014,0.013,0.01,0.012,0.02,0.034,0.049,0.059,0.066,0.068,0.073,0.084,0.098,0.106,0.104,0.098,0.092,0.087,0.076,0.06,0.041,0.024,0.011,0.004,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.008,0.017,0.031,0.042,0.045,0.04],[0.039,0.055,0.069,0.078,0.088,0.103,0.122,0.143,0.159,0.168,0.167,0.157,0.141,0.126,0.123,0.137,0.155,0.166,0.166,0.154,0.136,0.121,0.106,0.087,0.073,0.069,0.065,0.057,0.047,0.037,0.028,0.022,0.022,0.023,0.021,0.02,0.021,0.023,0.032,0.043,0.052,0.056,0.057,0.057,0.063,0.073,0.082,0.086,0.085,0.083,0.077,0.067,0.054,0.037,0.021,0.009,0.003,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.002,0.004,0.009,0.018,0.031,0.043,0.046,0.042],[0.042,0.058,0.072,0.084,0.097,0.111,0.123,0.137,0.148,0.158,0.165,0.163,0.153,0.141,0.137,0.147,0.164,0.174,0.172,0.161,0.144,0.125,0.106,0.089,0.081,0.08,0.073,0.061,0.05,0.042,0.036,0.031,0.03,0.029,0.027,0.029,0.031,0.029,0.031,0.04,0.049,0.055,0.056,0.052,0.049,0.055,0.066,0.074,0.078,0.074,0.065,0.054,0.043,0.031,0.018,0.008,0.003,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.008,0.013,0.019,0.027,0.035,0.037,0.034],[0.044,0.058,0.072,0.086,0.102,0.115,0.123,0.13,0.135,0.144,0.156,0.161,0.158,0.15,0.144,0.15,0.168,0.18,0.18,0.168,0.151,0.13,0.11,0.097,0.093,0.09,0.079,0.065,0.053,0.047,0.041,0.037,0.032,0.029,0.028,0.033,0.036,0.032,0.028,0.036,0.047,0.054,0.054,0.048,0.044,0.046,0.059,0.071,0.076,0.073,0.063,0.051,0.042,0.032,0.019,0.009,0.004,0.002,0.001,0.001,0.001,0.001,0.001,0.002,0.005,0.012,0.019,0.023,0.024,0.024,0.024,0.021],[0.046,0.062,0.076,0.088,0.103,0.116,0.123,0.128,0.132,0.143,0.157,0.167,0.167,0.162,0.158,0.164,0.181,0.194,0.191,0.175,0.155,0.134,0.117,0.109,0.107,0.099,0.082,0.064,0.052,0.046,0.041,0.036,0.031,0.029,0.027,0.029,0.031,0.027,0.023,0.03,0.042,0.049,0.048,0.044,0.041,0.041,0.051,0.064,0.072,0.074,0.067,0.056,0.047,0.037,0.023,0.012,0.008,0.006,0.004,0.003,0.003,0.002,0.002,0.002,0.006,0.014,0.023,0.027,0.025,0.02,0.018,0.016],[0.053,0.073,0.089,0.1,0.112,0.122,0.129,0.135,0.145,0.16,0.178,0.19,0.193,0.19,0.19,0.197,0.209,0.216,0.207,0.186,0.165,0.146,0.13,0.122,0.118,0.108,0.088,0.067,0.052,0.042,0.035,0.029,0.028,0.029,0.026,0.023,0.021,0.018,0.017,0.025,0.036,0.043,0.043,0.041,0.038,0.035,0.039,0.051,0.063,0.069,0.066,0.055,0.046,0.036,0.026,0.019,0.017,0.013,0.009,0.008,0.009,0.007,0.005,0.003,0.007,0.014,0.023,0.028,0.027,0.026,0.025,0.024],[0.062,0.082,0.1,0.115,0.129,0.138,0.143,0.153,0.169,0.188,0.207,0.222,0.228,0.23,0.233,0.237,0.242,0.241,0.226,0.205,0.186,0.168,0.152,0.139,0.129,0.117,0.099,0.076,0.058,0.044,0.031,0.023,0.025,0.029,0.028,0.025,0.022,0.02,0.023,0.029,0.038,0.044,0.044,0.04,0.035,0.031,0.034,0.043,0.053,0.061,0.058,0.048,0.038,0.033,0.032,0.031,0.029,0.023,0.017,0.017,0.019,0.016,0.01,0.006,0.008,0.015,0.025,0.031,0.033,0.034,0.035,0.034],[0.07,0.086,0.105,0.125,0.144,0.157,0.164,0.177,0.198,0.218,0.237,0.253,0.264,0.27,0.273,0.277,0.276,0.267,0.25,0.231,0.213,0.196,0.179,0.162,0.146,0.13,0.109,0.085,0.067,0.052,0.036,0.027,0.03,0.037,0.04,0.039,0.034,0.031,0.034,0.04,0.045,0.047,0.044,0.039,0.035,0.034,0.038,0.042,0.047,0.051,0.05,0.042,0.033,0.034,0.04,0.041,0.038,0.03,0.023,0.027,0.03,0.026,0.017,0.011,0.015,0.02,0.028,0.035,0.038,0.04,0.044,0.046],[0.082,0.092,0.11,0.132,0.155,0.174,0.189,0.207,0.227,0.247,0.266,0.284,0.297,0.305,0.31,0.312,0.307,0.292,0.275,0.257,0.241,0.227,0.21,0.19,0.169,0.147,0.12,0.093,0.076,0.064,0.051,0.041,0.043,0.052,0.056,0.052,0.042,0.037,0.039,0.045,0.046,0.043,0.037,0.034,0.037,0.043,0.046,0.047,0.048,0.05,0.049,0.043,0.036,0.038,0.045,0.048,0.043,0.033,0.026,0.031,0.035,0.03,0.02,0.018,0.023,0.027,0.029,0.031,0.036,0.044,0.053,0.058],[0.095,0.104,0.122,0.143,0.167,0.191,0.214,0.236,0.256,0.276,0.298,0.32,0.334,0.34,0.344,0.345,0.336,0.318,0.3,0.284,0.272,0.262,0.247,0.224,0.199,0.173,0.141,0.108,0.088,0.079,0.067,0.056,0.055,0.064,0.068,0.059,0.043,0.033,0.034,0.037,0.037,0.031,0.027,0.031,0.042,0.05,0.052,0.053,0.055,0.055,0.053,0.048,0.042,0.043,0.05,0.051,0.043,0.03,0.022,0.027,0.03,0.026,0.02,0.024,0.031,0.032,0.028,0.026,0.031,0.043,0.056,0.062],[0.103,0.116,0.136,0.158,0.182,0.208,0.237,0.263,0.284,0.306,0.334,0.36,0.375,0.379,0.378,0.374,0.36,0.341,0.324,0.313,0.306,0.299,0.284,0.261,0.238,0.212,0.175,0.133,0.105,0.091,0.078,0.064,0.059,0.067,0.07,0.059,0.039,0.024,0.022,0.026,0.026,0.023,0.026,0.037,0.048,0.054,0.053,0.055,0.059,0.058,0.053,0.048,0.043,0.045,0.05,0.049,0.04,0.026,0.016,0.017,0.019,0.017,0.017,0.025,0.032,0.032,0.029,0.029,0.034,0.043,0.053,0.059],[0.104,0.12,0.145,0.172,0.198,0.227,0.261,0.292,0.315,0.339,0.369,0.398,0.411,0.41,0.404,0.394,0.376,0.355,0.343,0.34,0.338,0.333,0.32,0.302,0.282,0.255,0.211,0.16,0.12,0.098,0.08,0.064,0.056,0.059,0.06,0.049,0.03,0.016,0.017,0.024,0.028,0.028,0.035,0.048,0.056,0.057,0.053,0.054,0.059,0.059,0.054,0.046,0.041,0.042,0.044,0.04,0.031,0.019,0.01,0.008,0.009,0.009,0.013,0.021,0.027,0.03,0.034,0.042,0.05,0.054,0.057,0.06],[0.106,0.124,0.153,0.184,0.213,0.246,0.285,0.321,0.347,0.37,0.397,0.421,0.429,0.424,0.413,0.399,0.38,0.363,0.358,0.363,0.366,0.363,0.354,0.342,0.323,0.289,0.236,0.177,0.13,0.099,0.076,0.06,0.054,0.052,0.047,0.035,0.021,0.012,0.018,0.028,0.034,0.034,0.04,0.052,0.059,0.058,0.053,0.052,0.057,0.06,0.057,0.048,0.041,0.04,0.036,0.028,0.02,0.011,0.005,0.003,0.003,0.004,0.008,0.014,0.021,0.029,0.04,0.052,0.063,0.067,0.067,0.067],[0.118,0.137,0.165,0.195,0.225,0.262,0.304,0.343,0.371,0.393,0.413,0.428,0.431,0.424,0.413,0.397,0.38,0.369,0.371,0.379,0.385,0.384,0.38,0.372,0.349,0.304,0.242,0.18,0.131,0.097,0.074,0.061,0.056,0.052,0.042,0.027,0.014,0.009,0.018,0.028,0.034,0.032,0.035,0.043,0.05,0.051,0.049,0.047,0.052,0.059,0.059,0.049,0.041,0.037,0.031,0.02,0.011,0.006,0.002,0.001,0.001,0.002,0.004,0.009,0.018,0.028,0.04,0.052,0.064,0.069,0.069,0.069],[0.134,0.154,0.178,0.203,0.235,0.275,0.317,0.354,0.385,0.406,0.421,0.428,0.426,0.419,0.407,0.391,0.375,0.369,0.373,0.381,0.387,0.39,0.39,0.381,0.352,0.3,0.234,0.172,0.125,0.093,0.074,0.063,0.057,0.049,0.037,0.022,0.01,0.007,0.014,0.022,0.026,0.024,0.023,0.027,0.032,0.036,0.038,0.039,0.046,0.055,0.055,0.047,0.039,0.034,0.025,0.015,0.007,0.003,0.002,0.002,0.001,0.002,0.004,0.009,0.017,0.026,0.035,0.043,0.053,0.059,0.063,0.065],[0.143,0.164,0.188,0.212,0.246,0.288,0.328,0.363,0.393,0.414,0.423,0.423,0.415,0.405,0.392,0.376,0.363,0.359,0.364,0.371,0.375,0.379,0.379,0.367,0.335,0.282,0.22,0.163,0.12,0.089,0.069,0.057,0.048,0.039,0.028,0.016,0.007,0.004,0.008,0.013,0.015,0.013,0.012,0.013,0.018,0.024,0.028,0.032,0.04,0.047,0.047,0.042,0.039,0.032,0.021,0.011,0.005,0.006,0.006,0.006,0.004,0.005,0.009,0.014,0.019,0.023,0.027,0.032,0.041,0.05,0.057,0.063],[0.143,0.169,0.198,0.227,0.262,0.303,0.339,0.37,0.399,0.416,0.419,0.412,0.397,0.379,0.364,0.352,0.344,0.345,0.352,0.358,0.361,0.363,0.36,0.343,0.309,0.259,0.203,0.153,0.115,0.087,0.063,0.045,0.033,0.024,0.016,0.009,0.004,0.002,0.004,0.006,0.007,0.006,0.005,0.008,0.015,0.025,0.032,0.034,0.035,0.038,0.039,0.038,0.036,0.029,0.018,0.009,0.008,0.012,0.014,0.013,0.01,0.011,0.016,0.023,0.026,0.024,0.023,0.027,0.037,0.049,0.061,0.07],[0.147,0.175,0.211,0.245,0.279,0.313,0.342,0.369,0.392,0.405,0.405,0.393,0.373,0.351,0.334,0.326,0.326,0.33,0.338,0.345,0.347,0.346,0.338,0.316,0.279,0.233,0.183,0.138,0.105,0.081,0.059,0.036,0.021,0.012,0.008,0.004,0.002,0.001,0.001,0.002,0.002,0.002,0.003,0.008,0.019,0.031,0.039,0.04,0.035,0.034,0.035,0.033,0.03,0.023,0.014,0.008,0.012,0.019,0.023,0.022,0.018,0.018,0.024,0.031,0.033,0.03,0.029,0.031,0.041,0.058,0.077,0.088],[0.158,0.186,0.224,0.259,0.288,0.311,0.332,0.352,0.368,0.377,0.375,0.363,0.345,0.325,0.309,0.302,0.304,0.309,0.316,0.325,0.328,0.323,0.31,0.282,0.244,0.202,0.158,0.118,0.09,0.071,0.051,0.031,0.015,0.006,0.003,0.001,0.001,0.0,0.0,0.0,0.001,0.002,0.004,0.01,0.019,0.031,0.039,0.04,0.038,0.037,0.035,0.031,0.025,0.018,0.011,0.008,0.014,0.023,0.028,0.028,0.026,0.024,0.027,0.033,0.036,0.037,0.038,0.04,0.049,0.071,0.094,0.107],[0.168,0.193,0.229,0.262,0.285,0.301,0.314,0.327,0.335,0.337,0.334,0.325,0.311,0.297,0.284,0.275,0.275,0.28,0.288,0.297,0.301,0.295,0.275,0.243,0.205,0.168,0.133,0.099,0.073,0.055,0.039,0.023,0.01,0.004,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.004,0.008,0.014,0.02,0.027,0.035,0.04,0.042,0.042,0.038,0.033,0.027,0.02,0.012,0.008,0.012,0.02,0.025,0.028,0.028,0.025,0.025,0.03,0.037,0.044,0.047,0.046,0.055,0.077,0.098,0.11],[0.168,0.189,0.22,0.25,0.271,0.286,0.296,0.302,0.304,0.3,0.292,0.284,0.274,0.264,0.254,0.247,0.246,0.252,0.261,0.269,0.271,0.262,0.239,0.205,0.17,0.14,0.113,0.086,0.061,0.041,0.026,0.014,0.006,0.002,0.0,0.0,0.0,0.0,0.0,0.001,0.002,0.006,0.014,0.022,0.028,0.03,0.036,0.044,0.049,0.047,0.041,0.034,0.029,0.023,0.014,0.007,0.008,0.013,0.018,0.022,0.023,0.02,0.018,0.025,0.036,0.046,0.05,0.05,0.056,0.07,0.085,0.095],[0.16,0.177,0.202,0.227,0.247,0.266,0.277,0.28,0.278,0.271,0.26,0.249,0.24,0.232,0.224,0.219,0.219,0.226,0.238,0.244,0.241,0.227,0.202,0.171,0.141,0.118,0.1,0.08,0.056,0.033,0.017,0.008,0.003,0.001,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.009,0.019,0.03,0.036,0.036,0.039,0.049,0.055,0.055,0.045,0.032,0.025,0.019,0.012,0.006,0.004,0.006,0.009,0.013,0.014,0.012,0.011,0.018,0.03,0.039,0.046,0.051,0.056,0.06,0.066,0.071],[0.152,0.167,0.187,0.206,0.227,0.25,0.263,0.262,0.259,0.252,0.241,0.229,0.218,0.209,0.204,0.201,0.201,0.209,0.221,0.224,0.217,0.2,0.174,0.146,0.122,0.104,0.092,0.076,0.054,0.031,0.014,0.005,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.003,0.01,0.022,0.035,0.041,0.039,0.039,0.049,0.059,0.061,0.05,0.032,0.02,0.013,0.008,0.004,0.002,0.002,0.004,0.006,0.007,0.006,0.006,0.012,0.02,0.029,0.039,0.051,0.057,0.055,0.054,0.057]];
const RAMA_COLORSCALE=[[0,'#050200'],[0.06,'#3a1a00'],[0.15,'#7a3600'],[0.28,'#b35a00'],[0.42,'#e8821a'],[0.55,'#ffa62b'],[0.7,'#ffc94d'],[0.85,'#ffe38a'],[1,'#fff6cf']];
function drawRamaHeat(){
  const data=[{x:RAMA_PHI,y:RAMA_PHI,z:RAMA_DENSITY,type:'heatmap',colorscale:RAMA_COLORSCALE,reversescale:false,zsmooth:'best',showscale:false,hoverinfo:'skip',zmin:0,zmax:1}];
  const layout={
    title:{text:'Ramachandran plot — click to animate',font:{size:14, color:'white'}},
    xaxis:{title:{text:'Phi (°)',font:{color:'white'}},range:[-180,180],dtick:60,gridcolor:'#333',zeroline:true,zerolinecolor:'#666',tickcolor:'#888',tickfont:{color:'white'},linecolor:'#888'},
    yaxis:{title:{text:'Psi (°)',font:{color:'white'}},range:[-180,180],dtick:60,gridcolor:'#333',zeroline:true,zerolinecolor:'#666',tickcolor:'#888',tickfont:{color:'white'},linecolor:'#888'},
    annotations:[
      {x:-60,y:-45,text:'α helix',showarrow:true,xref:'x',yref:'y',axref:'x',ayref:'y',ax:-15,ay:-100,arrowhead:2,arrowsize:1,arrowwidth:1.2,arrowcolor:'white',font:{color:'white',size:12},bgcolor:'rgba(0,0,0,0.55)',borderpad:2},
      {x:-120,y:130,text:'β sheet',showarrow:true,xref:'x',yref:'y',axref:'x',ayref:'y',ax:-160,ay:170,arrowhead:2,arrowsize:1,arrowwidth:1.2,arrowcolor:'white',font:{color:'white',size:12},bgcolor:'rgba(0,0,0,0.55)',borderpad:2},
      {x:60,y:40,text:'left-handed',showarrow:true,xref:'x',yref:'y',axref:'x',ayref:'y',ax:110,ay:95,arrowhead:2,arrowsize:1,arrowwidth:1.2,arrowcolor:'white',font:{color:'white',size:10},bgcolor:'rgba(0,0,0,0.55)',borderpad:2}
    ],
    margin:{t:40,l:60,r:20,b:50}, height:520, hovermode:'closest', plot_bgcolor:'black', paper_bgcolor:'black'
  };
  Plotly.newPlot('ramaplot',data,layout,{displayModeBar:false,responsive:true}).then(()=>{
    Plotly.addTraces('ramaplot',[
      {x:[currentPhi],y:[currentPsi],mode:'markers',marker:{color:'red',size:10,line:{color:'white',width:1}},name:'current',hoverinfo:'skip',showlegend:false},
      {x:[currentPhi,currentPhi],y:[-180,180],mode:'lines',line:{color:'red',width:1,dash:'dot'},hoverinfo:'none',showlegend:false},
      {x:[-180,180],y:[currentPsi,currentPsi],mode:'lines',line:{color:'red',width:1,dash:'dot'},hoverinfo:'none',showlegend:false}
    ]);
    // 5 white contour levels as a dedicated contour trace (heatmap contours
    // don't render reliably with zsmooth:'best'), added last so the red
    // marker/crosshair trace indices (1,2,3) stay intact.
    Plotly.addTraces('ramaplot',[
      {x:RAMA_PHI,y:RAMA_PHI,z:RAMA_DENSITY,type:'contour',showscale:false,hoverinfo:'skip',zmin:0,zmax:1,
       contours:{coloring:'lines',showlines:true,showlabels:true,start:0.1,end:0.9,size:0.2,labelfont:{color:'white',size:10}},
       line:{color:'white',width:1.2},opacity:0.9}
    ]);
    const plot=document.getElementById('ramaplot');
    // helper to convert pixel event to phi/psi using Plotly axis conversion (correct y flip)
    function pixelToPhiPsi(evt){
      try{
        const xaxis=plot._fullLayout.xaxis, yaxis=plot._fullLayout.yaxis;
        if(!xaxis || !yaxis || !xaxis.p2d || !yaxis.p2d) return null;
        const rect=plot.getBoundingClientRect();
        const ex=evt.clientX - rect.left, ey=evt.clientY - rect.top;
        if(ex < xaxis._offset || ex > xaxis._offset + xaxis._length) return null;
        if(ey < yaxis._offset || ey > yaxis._offset + yaxis._length) return null;
        const phi=xaxis.p2d(ex - xaxis._offset);
        const psi=yaxis.p2d(ey - yaxis._offset);
        if(phi<-180||phi>180||psi<-180||psi>180) return null;
        return {phi,psi};
      }catch(e){ return null; }
    }
    // single handler for plotly_click: use pixel conversion for continuous, accurate phi+psi
    let lastClickTime=0;
    plot.on('plotly_click',function(data){
      // debounce duplicate with manual click
      lastClickTime=Date.now();
      const evt=data.event;
      let coords=pixelToPhiPsi(evt);
      // fallback to heatmap point if pixel conversion fails (e.g., p2d not ready)
      if(!coords && data.points && data.points.length){
        const pt=data.points.find(p=>p.curveNumber===0) || data.points[0];
        if(pt && typeof pt.x==='number' && typeof pt.y==='number'){
          coords={phi:pt.x, psi:pt.y};
        }
      }
      if(coords){
        coords.phi=Math.max(-180,Math.min(180,coords.phi));
        coords.psi=Math.max(-180,Math.min(180,coords.psi));
        moveToPhiPsi(coords.phi, coords.psi);
      }
    });
    // also handle clicks that don't trigger plotly_click (e.g., on empty margin) – but only if not recently handled
    plot.addEventListener('click',function(evt){
      if(animating) return;
      if(Date.now()-lastClickTime < 400) return; // already handled by plotly_click
      const coords=pixelToPhiPsi(evt);
      if(coords){
        moveToPhiPsi(coords.phi, coords.psi);
      }
    });
  });
}
function updatePlotMarker(phi,psi){
  if(isNaN(phi)||isNaN(psi)) return;
  const plot=document.getElementById('ramaplot');
  if(!plot || !plot.data) return;
  try{
    Plotly.restyle('ramaplot',{x:[[phi]],y:[[psi]]},[1]);
    Plotly.restyle('ramaplot',{x:[[phi,phi]],y:[[-180,180]]},[2]);
    Plotly.restyle('ramaplot',{x:[[-180,180]],y:[[psi,psi]]},[3]);
  }catch(e){}
}
// Single smooth gradual move: one continuous interpolation, not many macro steps
function moveToPhiPsi(targetPhi,targetPsi){
  if(isNaN(targetPhi)||isNaN(targetPsi)) return;
  if(animating) return;
  const cur=getPhiPsi();
  if(isNaN(cur.phi)||isNaN(cur.psi)) return;
  let dPhi=normalizeAngle(targetPhi-cur.phi);
  let dPsi=normalizeAngle(targetPsi-cur.psi);
  const maxDelta=Math.max(Math.abs(dPhi),Math.abs(dPsi));
  if(maxDelta < 0.2){ updateDihedralsDisplay(); return; }
  // choose frames for smooth 60fps-like: ~1.5° per frame gives 15°=>10 frames, 180°=>120 frames capped
  // request is single smooth move, so we do e.g., 30-50 frames total
  const frames=Math.max(18, Math.min(55, Math.ceil(maxDelta/2.5)));
  const dPhiStep=dPhi/frames;
  const dPsiStep=dPsi/frames;
  animating=true;
  let i=0;
  function step(){
    if(i>=frames){
      animating=false;
      // final correction to exact target (avoid drift)
      const cur2=getPhiPsi();
      const remPhi=normalizeAngle(targetPhi-cur2.phi);
      const remPsi=normalizeAngle(targetPsi-cur2.psi);
      if(Math.abs(remPhi)>0.3 || Math.abs(remPsi)>0.3){
        if(Math.abs(remPhi)>0.3){
          const N=findAtom(16,'N'), CA=findAtom(16,'CA');
          if(N&&CA) rotateAtomsAboutAxis({x:N.x,y:N.y,z:N.z},{x:CA.x,y:CA.y,z:CA.z}, -remPhi, isMovingAtomForPhi);
        }
        if(Math.abs(remPsi)>0.3){
          const CA=findAtom(16,'CA'), C=findAtom(16,'C');
          if(CA&&C) rotateAtomsAboutAxis({x:CA.x,y:CA.y,z:CA.z},{x:C.x,y:C.y,z:C.z}, remPsi, isMovingAtomForPsi);
        }
      }
      updateDihedralsDisplay();
      viewer.render();
      if(showPlanes) updatePlanes();
      if(showClashes) updateClashes();
      return;
    }
    // rotate a little phi and a little psi in same frame – do both without intermediate highlight
    const N=findAtom(16,'N'), CA=findAtom(16,'CA'), C=findAtom(16,'C');
    if(Math.abs(dPhiStep) > 0.001 && N && CA){
      const axisA={x:N.x,y:N.y,z:N.z}, axisB={x:CA.x,y:CA.y,z:CA.z};
      const ang=-dPhiStep*Math.PI/180; // invert for upstream phi (N-terminal)
      const u={x:axisB.x-axisA.x,y:axisB.y-axisA.y,z:axisB.z-axisA.z};
      const un=normalize(u); const cosA=Math.cos(ang), sinA=Math.sin(ang);
      atomsData.forEach(a=>{
        if(!isMovingAtomForPhi(a)) return;
        const p={x:a.x-axisA.x,y:a.y-axisA.y,z:a.z-axisA.z};
        const dotu=dot(un,p), crossu=cross(un,p);
        const rot={x:un.x*dotu*(1-cosA)+p.x*cosA+crossu.x*sinA, y:un.y*dotu*(1-cosA)+p.y*cosA+crossu.y*sinA, z:un.z*dotu*(1-cosA)+p.z*cosA+crossu.z*sinA};
        a.x=rot.x+axisA.x; a.y=rot.y+axisA.y; a.z=rot.z+axisA.z;
      });
    }
    if(Math.abs(dPsiStep) > 0.001 && CA && C){
      const CA2=findAtom(16,'CA'), C2=findAtom(16,'C');
      if(CA2&&C2){
        const axisA={x:CA2.x,y:CA2.y,z:CA2.z}, axisB={x:C2.x,y:C2.y,z:C2.z};
        const ang=dPsiStep*Math.PI/180;
        const u={x:axisB.x-axisA.x,y:axisB.y-axisA.y,z:axisB.z-axisA.z};
        const un=normalize(u); const cosA=Math.cos(ang), sinA=Math.sin(ang);
        atomsData.forEach(a=>{
          if(!isMovingAtomForPsi(a)) return;
          const p={x:a.x-axisA.x,y:a.y-axisA.y,z:a.z-axisA.z};
          const dotu=dot(un,p), crossu=cross(un,p);
          const rot={x:un.x*dotu*(1-cosA)+p.x*cosA+crossu.x*sinA, y:un.y*dotu*(1-cosA)+p.y*cosA+crossu.y*sinA, z:un.z*dotu*(1-cosA)+p.z*cosA+crossu.z*sinA};
          a.x=rot.x+axisA.x; a.y=rot.y+axisA.y; a.z=rot.z+axisA.z;
        });
      }
    }
    syncModelPositions();
    updateDoubleBonds();
    if(showPlanes) updatePlanes();
    if(showPeptideBondsFlag) doPeptideBonds(true);
    if(showVDW){
      vdwShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
      vdwShapes=[];
      const radii={H:1.20, C:1.70, N:1.55, O:1.52};
      const isWhite=document.getElementById('idwhite')?.checked;
      atomsData.forEach(a=>{
        if(showClashes && clashAtomSerials.has(a.serial)) return; // lens replaces this atom's VdW
        const ra=radii[a.elem]||1.5;
        let col=isWhite?'white':(a.atom==='CA'?'#000000':a.elem==='N'?'#3050ff':a.elem==='O'?'#ff2020':a.elem==='H'?'white':'#c8c8c8');
        const op=isWhite?0.28:0.38;
        const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:ra*0.88, color:col, opacity:showClashes?(isWhite?0.16:0.22):op});
        vdwShapes.push(s);
      });
      if(showClashes) updateClashes();
    } else if(showClashes){
      updateClashes();
    }
    if(alanineShapes.length>0){
      alanineShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
      alanineShapes=[];
      const alaAtoms=atomsData.filter(a=>a.resi===16);
      alaAtoms.forEach(a=>{
        const s=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:0.9, color:'black', opacity:0.18});
        alanineShapes.push(s);
      });
    }
    const curStep=getPhiPsi();
    updatePlotMarker(curStep.phi, curStep.psi);
    if(i % 3 ===0){
      if(Math.abs(dPhi) > Math.abs(dPsi)) highlightPhiPsi('phi'); else highlightPhiPsi('psi');
    }
    i++;
    requestAnimationFrame(step);
  }
  requestAnimationFrame(step);
}
function normalizeAngle(a){ while(a>180) a-=360; while(a<-180) a+=360; return a; }

document.addEventListener('DOMContentLoaded',()=>{
  initViewer();
  loadTripeptide();
  document.getElementById('phi-plus').addEventListener('click',()=> adjustDihedral('phi', 15));
  document.getElementById('phi-minus').addEventListener('click',()=> adjustDihedral('phi', -15));
  document.getElementById('psi-plus').addEventListener('click',()=> adjustDihedral('psi', 15));
  document.getElementById('psi-minus').addEventListener('click',()=> adjustDihedral('psi', -15));
  const resiInput=document.getElementById('resi-input');
  if(resiInput) resiInput.addEventListener('change',updateDihedralsDisplay);
  const loadBtn=document.getElementById('load-pdb');
  if(loadBtn) loadBtn.addEventListener('click',()=> loadPDB('1CRN'));
  document.getElementById('idalanine').addEventListener('change',e=> doAlanine(e.target.checked));
  document.getElementById('idpeptidebonds').addEventListener('change',e=> doPeptideBonds(e.target.checked));
  document.getElementById('idplanes').addEventListener('change',e=> doPlanes(e.target.checked));
  document.getElementById('idvdw').addEventListener('change',e=> doVDW(e.target.checked));
  document.getElementById('idwhite').addEventListener('change',e=> doWhite(e.target.checked));
  document.getElementById('idclashes').addEventListener('change',e=> doClashes(e.target.checked));
  document.getElementById('idtrailclashes').addEventListener('change',e=> doTrailClashes(e.target.checked));
  document.getElementById('resetBtn').addEventListener('click',resetViewer);
  window.addEventListener('resize',()=>{ try{ Plotly.Plots.resize('ramaplot'); viewer.resize(); }catch(e){} });
});
