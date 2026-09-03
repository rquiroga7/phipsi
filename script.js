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
  viewer.render();
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
  model=viewer.addModel(pdbText,'pdb');
  applyOriginalStyle(); // adds doubleBondShapes
  if(view) try{ viewer.setView(view); }catch(e){ viewer.zoomTo(); }
  else viewer.zoomTo();
  viewer.render();
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
    model=viewer.addModel(generatePDB(atomsData),'pdb');
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
    model=viewer.addModel(generatePDB(atomsData),'pdb');
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
    model=viewer.addModel(generatePDB(atomsData),'pdb');
    viewer.setStyle({}, {cartoon:{color:'spectrum'}});
    viewer.zoomTo();
    viewer.render();
    updateDihedralsDisplay();
    drawRamaHeat();
  }).catch(e=>alert('Failed to load PDB:'+e));
}

function applyOriginalStyle(){
  // base ball-and-stick: stick radius ~0.12, sphere scale 0.30, white bonds
  viewer.setStyle({}, {stick:{radius:0.12, color:'white'}, sphere:{scale:0.30, colorscheme:'whiteCarbon'}});
  // now override colors via addStyle to keep scale
  viewer.addStyle({elem:'C'}, {sphere:{color:'#c8c8c8'}});
  viewer.addStyle({atom:'CA'}, {sphere:{color:'#000000'}}); // Calpha black as requested (original key #383838, now pure black)
  viewer.addStyle({elem:'N'}, {sphere:{color:'#3050ff'}});
  viewer.addStyle({elem:'O'}, {sphere:{color:'#ff2020'}});
  viewer.addStyle({elem:'H'}, {sphere:{color:'white', scale:0.20}});
  viewer.addStyle({}, {stick:{color:'white', radius:0.12}});
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
    // Overlay translucent VdW spheres without hiding ball-and-stick (original keeps both models)
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
function makeClashPair(a,b){
  if(!a||!b) return;
  const radii={H:1.20, C:1.70, N:1.55, O:1.52, S:1.80};
  const scale=0.88;
  const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
  const d=Math.hypot(dx,dy,dz);
  if(d>4.5) return;
  // skip directly bonded (1-2) – distance ~1.0-1.5 but they are bonded; original excludes bonded via contact definition
  // Use a simple bonded check: if they share a bond (distance <1.65 and are 1-2), skip? For tripeptide, bonded pairs are: 15 CA-C, C-O, C-N, N-CA, CA-C, C-O, C-N, N-CA, CA-CB, etc. We skip d<1.6 for bonded
  if(d<1.65){
    // check if they are known bonded pairs in tripeptide – skip
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
  if(d < sum){
    const overlap=sum-d;
    const mid={x:(a.x+b.x)/2, y:(a.y+b.y)/2, z:(a.z+b.z)/2};
    // red lens at overlap – non-transparent as requested
    const r=Math.min(0.85, 0.30 + overlap*0.7);
    const s=viewer.addSphere({center:mid, radius:r, color:'red', opacity:0.98});
    clashShapes.push(s);
    // also add a small red cap on each atom's vdW surface for emphasis
    const s2=viewer.addSphere({center:{x:a.x,y:a.y,z:a.z}, radius:ra*scale*0.42, color:'red', opacity:0.92});
    const s3=viewer.addSphere({center:{x:b.x,y:b.y,z:b.z}, radius:rb*scale*0.42, color:'red', opacity:0.92});
    clashShapes.push(s2,s3);
  }
}
function updateClashes(){
  if(!showClashes) return;
  if(!trailClashes){
    clashShapes.forEach(s=>{ try{ viewer.removeShape(s);}catch(e){} });
    clashShapes=[];
  }
  // Use original's 5 clash checks (88%)
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

// Ramachandran heatplot
function drawRamaHeat(){
  const size=121, xs=[], ys=[];
  for(let i=0;i<size;i++) xs.push(-180+i*3);
  for(let j=0;j<size;j++) ys.push(-180+j*3);
  const z=[];
  for(let j=0;j<size;j++){
    const row=[];
    for(let i=0;i<size;i++){
      const x=xs[i], y=ys[j];
      let v=0;
      v+=gaussian2(x,y,-60,-45,25,25);
      v+=gaussian2(x,y,-120,130,30,28);
      v+=gaussian2(x,y,-120,-150,22,18);
      v+=gaussian2(x,y,60,40,20,20);
      v+=gaussian2(x,y,-70,140,18,18)*0.5;
      v+=gaussian2(x,y,50,-130,15,15)*0.4;
      row.push(v);
    }
    z.push(row);
  }
  const data=[{x:xs,y:ys,z:z,type:'heatmap',colorscale:'YlOrRd',reversescale:false,zsmooth:'best',showscale:false,hoverinfo:'skip'}];
  const layout={
    title:{text:'Ramachandran plot — click to animate',font:{size:14}},
    xaxis:{title:'Phi (°)',range:[-180,180],dtick:60,gridcolor:'#ddd',zeroline:true,zerolinecolor:'#999'},
    yaxis:{title:'Psi (°)',range:[-180,180],dtick:60,gridcolor:'#ddd',zeroline:true,zerolinecolor:'#999'},
    margin:{t:40,l:60,r:20,b:50}, height:520, hovermode:'closest', plot_bgcolor:'#f8f8f8', paper_bgcolor:'white'
  };
  Plotly.newPlot('ramaplot',data,layout,{displayModeBar:false,responsive:true}).then(()=>{
    Plotly.addTraces('ramaplot',[
      {x:[currentPhi],y:[currentPsi],mode:'markers',marker:{color:'red',size:10,line:{color:'black',width:1}},name:'current',hoverinfo:'skip',showlegend:false},
      {x:[currentPhi,currentPhi],y:[-180,180],mode:'lines',line:{color:'red',width:1,dash:'dot'},hoverinfo:'none',showlegend:false},
      {x:[-180,180],y:[currentPsi,currentPsi],mode:'lines',line:{color:'red',width:1,dash:'dot'},hoverinfo:'none',showlegend:false}
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
function gaussian2(x,y,muX,muY,sx,sy){ return Math.exp(-((x-muX)*(x-muX)/(2*sx*sx)+(y-muY)*(y-muY)/(2*sy*sy))); }
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
    // update marker smoothly (model and all extras already re-added via rebuild)
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
