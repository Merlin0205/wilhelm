// 3D prohlížeč velké U-47 1:50: rozložení na díly, řez, skrývání dílů, výběr kliknutím.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';

const box = document.getElementById('v3d');
const canvas = box.querySelector('canvas');
const status = box.querySelector('.v3d-status');
const $ = (s) => document.querySelector(s);

// ---------- popisy dílů (jméno uzlu v GLB -> název, skupina, poznámka, stav)
const GROUPS = [
  ['spodek', 'Trup – spodek'], ['vrch', 'Trup – horní záď a příď'], ['poklop', 'Poklop s věží (snímatelný)'],
  ['vyzbroj', 'Výzbroj a výstroj paluby'], ['vnitrek', 'Vnitřek: přepážky, sedla, vedení'], ['pohon', 'Pohon a kormidla'], ['modul', 'Modul Holland (máte)'],
];
const RULES = [
  [/^spodek_1/, 'Spodek trupu 1 (záď)', 'spodek', 'Lepí se epoxidem, na přední straně je vnitřní límec pro další díl.', 'lepí se'],
  [/^spodek_6/, 'Spodek trupu 6 (příď)', 'spodek', 'Poslední díl spodku, nasune se na límec dílu 5.', 'lepí se'],
  [/^spodek_(\d)/, 'Spodek trupu $1', 'spodek', 'Díl spodku trupu, slepený přes vnitřní límec.', 'lepí se'],
  [/^vrch_zad_(\d)/, 'Horní záď $1', 'vrch', 'Pevná horní část zádě, lem zapadne do spodku.', 'lepí se'],
  [/^vrch_prid_(\d)/, 'Horní příď $1', 'vrch', 'Pevná horní část přídě s řadami zatopovacích otvorů.', 'lepí se'],
  [/^poklop_(\d)/, 'Poklop $1', 'poklop', 'Poklop nad modulem. Tři díly se slepí do jednoho, drží na 4 magnetech a zvedá se kolmo nahoru.', 'snímá se'],
  [/^vez/, 'Věž', 'poklop', 'Věž podle stavu z října 1939: vysoká přední část, zadní můstek se zábradlím, býk na obou bocích.', 'lepí se'],
  [/^podlaha_mustku/, 'Podlaha předního můstku', 'poklop', 'Periskopy, UZO, kompas a rámová anténa. Vkládá se do věže.', 'lepí se'],
  [/^podlaha_zadni/, 'Podlaha zadního můstku', 'poklop', 'Otevřený zadní můstek.', 'lepí se'],
  [/^periskop/, 'Periskop', 'poklop', 'Mosazná trubička Ø4 zasunutá do podlahy můstku.', 'lepí se'],
  [/^delo_lafeta/, 'Dělo 8,8 cm – lafeta', 'vyzbroj', 'Na palubě před věží, kolík z filamentu 1,75 do otvoru.', 'lepí se'],
  [/^delo_hlaven/, 'Dělo 8,8 cm – hlaveň', 'vyzbroj', 'Zasune se do lafety.', 'lepí se'],
  [/^lafeta_flak/, 'Lafeta kulometu 2 cm', 'vyzbroj', 'Na zadní palubě. V říjnu 1939 jen lafeta, zbraň se nevozila.', 'lepí se'],
  [/^rezac_siti/, 'Řezač sítí', 'vyzbroj', 'Na přídi, drží přední konec lana proti sítím.', 'lepí se'],
  [/^sloupek_lana/, 'Sloupek lana', 'vyzbroj', 'Zadní konec lana proti sítím.', 'lepí se'],
  [/^lano/, 'Lano proti sítím', 'vyzbroj', 'Rybářská šňůra 0,5 mm přes otvory ve věži.', 'koupit'],
  [/^uvaznik/, 'Úvazník', 'vyzbroj', 'Pár sloupků na palubě, 4 kusy.', 'lepí se'],
  [/^navijak/, 'Kotevní naviják', 'vyzbroj', 'Na přídi.', 'lepí se'],
  [/^prid_kormidlo/, 'Příďové kormidlo', 'vyzbroj', 'Pevné, jazýček do štěrbiny v přídi.', 'lepí se'],
  [/^prepazka_zad/, 'Přepážka záď', 'vnitrek', 'Vede táhla u zádě. Tiskne se naplocho, lepí se do drážky mezi žebírka.', 'lepí se'],
  [/^prepazka_hridel/, 'Přepážka hřídelí', 'vnitrek', 'Nese ložiska 624ZZ obou hřídelí a trubičky táhel.', 'lepí se'],
  [/^prepazka_motor(\d)/, 'Přepážka motoru $1', 'vnitrek', 'Drží motory, ty se z ní vysunou dopředu.', 'lepí se'],
  [/^prepazka_doraz/, 'Doraz zadního víka', 'vnitrek', 'Modul se o ni opře zadním víkem.', 'lepí se'],
  [/^prepazka_predni/, 'Doraz předního víka', 'vnitrek', 'Drží modul vpředu.', 'lepí se'],
  [/^sedlo(\d)/, 'Sedlo modulu $1', 'vnitrek', 'Modul v něm leží, O‑kroužek 78 × 3,5 se zahákne za kolíky.', 'lepí se'],
  [/^sani/, 'Sání čerpadla', 'vnitrek', 'Zevnitř nad otvorem ve dně, hadička 6 mm z vývodu modulu.', 'lepí se'],
  [/^trubicka/, 'Vedení táhla', 'vnitrek', 'Mosazná trubička 4 × 3,1 mm, v ní drát 1,5 mm.', 'koupit'],
  [/^stabilizator/, 'Stabilizátor', 'pohon', 'Vodorovná ploutev se žebrem, nese čepy svislých kormidel.', 'lepí se'],
  [/^kormidlo_horni/, 'Svislé kormidlo – horní list', 'pohon', 'Na mosazném čepu Ø3, otáčí se.', 'hýbe se'],
  [/^kormidlo_dolni/, 'Svislé kormidlo – dolní list', 'pohon', 'Na mosazném čepu Ø3, otáčí se.', 'hýbe se'],
  [/^paka_kormidla/, 'Páka svislého kormidla', 'pohon', 'Obě páky spojuje drát 1,5 mm, ovládá je spodní táhlo.', 'hýbe se'],
  [/^hloubkove_kormidlo/, 'Hloubkové kormidlo', 'pohon', 'List na hřídeli Ø4 napříč zádí.', 'hýbe se'],
  [/^paka_hloubkovych/, 'Páka hloubkových kormidel', 'pohon', 'Uvnitř zádě na hřídeli, šroubek M2, ovládá ji horní táhlo.', 'hýbe se'],
  [/^konzola_hridele/, 'Konzola hřídele', 'pohon', 'Ložisko 624ZZ venku před vrtulí.', 'lepí se'],
  [/^vrtule/, 'Lodní šroub Ø32', 'pohon', 'Třílistý, levý a pravý. Na hřídeli stavěcím šroubkem M3.', 'hýbe se'],
  [/^hridel/, 'Hřídel Ø4', 'pohon', 'Mosazná tyčka, dvě ložiska 624ZZ.', 'hýbe se'],
  [/^motor/, 'Motor LICHIFIT 7,4 V', 'pohon', 'Ø31 × 50 mm, vysune se dopředu pod poklop.', 'snímá se'],
  [/^modul_trubka/, 'Modul Holland – plexi trubka', 'modul', 'Ø90 × 440 mm, uvnitř elektronika, serva a čerpadlo s vakem.', 'snímá se'],
  [/^modul_vickoAR/, 'Zadní víko modulu', 'modul', 'Vývody: 2 táhla v měších, kabel motorů, sání čerpadla, anténa.', 'snímá se'],
  [/^modul_vickoAV/, 'Přední víko modulu', 'modul', 'Hladké víko s O‑kroužkem.', 'snímá se'],
];
function describe(name) {
  for (const [re, label, g, note, state] of RULES) { const m = name.match(re); if (m) return { label: label.replace('$1', m[1] || ''), group: g, note, state }; }
  return { label: name, group: 'vnitrek', note: '', state: '' };
}
// rozložení: posun v mm (scéna: x = délka, y = nahoru, z = -bok)
function explodeVec(name, c) {
  const spread = (c.x - 665) * 0.16;
  const d = describe(name).group;
  if (d === 'spodek') return new THREE.Vector3(spread, -70, 0);
  if (d === 'vrch') return new THREE.Vector3(spread, 70, 0);
  if (d === 'poklop') return new THREE.Vector3(spread * 0.6, /^(vez|podlaha|periskop)/.test(name) ? 230 : 150, 0);
  if (d === 'vyzbroj') return new THREE.Vector3(spread, /^(delo|lafeta)/.test(name) ? 200 : 110, 0);
  if (d === 'pohon') return new THREE.Vector3(-90, -15, c.z * 0.8);
  if (d === 'modul') return new THREE.Vector3(spread * 0.4, 40, 0);
  return new THREE.Vector3(spread, -20, 0);
}

// ---------- scéna
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2)); renderer.localClippingEnabled = true;
renderer.toneMapping = THREE.ACESFilmicToneMapping; renderer.outputColorSpace = THREE.SRGBColorSpace;
const scene = new THREE.Scene(); scene.background = new THREE.Color(0x0a1a28);
scene.environment = new THREE.PMREMGenerator(renderer).fromScene(new RoomEnvironment(), 0.04).texture;
const sun = new THREE.DirectionalLight(0xffffff, 1.6); sun.position.set(400, 900, 700); scene.add(sun);
scene.add(new THREE.HemisphereLight(0xbfd8ff, 0x1a2430, 0.6));
const camera = new THREE.PerspectiveCamera(32, 1, 5, 20000);
const controls = new OrbitControls(camera, canvas); controls.enableDamping = true; controls.autoRotate = true; controls.autoRotateSpeed = 0.8;
controls.minDistance = 60; controls.maxDistance = 6000;
const CENTER = new THREE.Vector3(665, 75, 0);
const VIEWS = {
  celek: [new THREE.Vector3(1650, 620, 1350), CENTER], bok: [new THREE.Vector3(665, 80, 2600), CENTER], shora: [new THREE.Vector3(665, 2700, 1), CENTER],
  zad: [new THREE.Vector3(-260, 170, 330), new THREE.Vector3(90, 60, 0)], vez: [new THREE.Vector3(880, 330, 430), new THREE.Vector3(690, 160, 0)],
  vnitrek: [new THREE.Vector3(700, 520, 900), new THREE.Vector3(430, 50, 0)],
};
let camAnim = null;
function viewPos(k) {           // na úzké obrazovce kameru odsunout, aby se vešla celá délka
  const [p, t] = VIEWS[k]; const f = Math.max(1, 1.45 / Math.max(0.3, camera.aspect));
  return [t.clone().add(p.clone().sub(t).multiplyScalar(f)), t];
}
function goView(k) {
  const [p, t] = viewPos(k); camAnim = { p0: camera.position.clone(), t0: controls.target.clone(), p1: p.clone(), t1: t.clone(), s: 0 };
  setAuto(false);
}

// ---------- model
const parts = []; // {mesh, name, info, base: Vector3, ex: Vector3, mats}
const clip = new THREE.Plane(new THREE.Vector3(0, 0, -1), 1e5); const clipOn = [clip];
const pick = new THREE.Raycaster();
let explode = 0, explodeTarget = 0, selected = null;
new GLTFLoader().load('velka/models/u47_50.glb', (gltf) => {
  gltf.scene.traverse((o) => {
    if (!o.isMesh) return;
    const name = (o.name || o.parent?.name || '').replace(/\.\d+$/, '');
    o.material = o.material.clone(); o.material.side = THREE.DoubleSide; o.material.clippingPlanes = clipOn;
    if (/^modul_trubka/.test(name)) { o.material.transparent = true; o.material.opacity = 0.28; o.material.depthWrite = false; }
    const c = new THREE.Box3().setFromObject(o).getCenter(new THREE.Vector3());
    parts.push({ mesh: o, name, info: describe(name), base: o.position.clone(), ex: explodeVec(name, c), center: c });
  });
  scene.add(gltf.scene); buildList(); status.hidden = true; resize(); goViewNow('celek');
}, (e) => { if (e.total) status.textContent = `Načítám 3D model… ${Math.round(100 * e.loaded / e.total)} %`; }, () => { status.textContent = 'Model se nepodařilo načíst.'; });
function goViewNow(k) { const [p, t] = viewPos(k); camera.position.copy(p); controls.target.copy(t); controls.update(); }

// rodičovské transformace: posun v souřadnicích scény převedeme do lokálních souřadnic rodiče
function applyExplode() {
  for (const p of parts) {
    const w = p.ex.clone().multiplyScalar(explode); const parent = p.mesh.parent;
    const inv = new THREE.Matrix4().copy(parent.matrixWorld).invert(); const o = new THREE.Vector3(0, 0, 0).applyMatrix4(parent.matrixWorld);
    const target = o.clone().add(w); const l0 = new THREE.Vector3(0, 0, 0).applyMatrix4(inv), l1 = target.applyMatrix4(inv);
    p.mesh.position.copy(p.base).add(l1.sub(l0));
  }
}

// ---------- ovládání
function setAuto(on) { controls.autoRotate = on; $('#v3dRot').classList.toggle('on', on); $('#v3dRot').textContent = on ? 'Otáčení: zap' : 'Otáčení: vyp'; }
controls.addEventListener('start', () => { setAuto(false); camAnim = null; });
$('#v3dRot').onclick = () => setAuto(!controls.autoRotate);
$('#v3dExplode').onclick = () => { explodeTarget = explodeTarget > 0.5 ? 0 : 1; $('#v3dExplode').textContent = explodeTarget ? 'Složit' : 'Rozložit na díly'; };
$('#v3dExSlider').oninput = (e) => { explodeTarget = explode = e.target.value / 100; $('#v3dExplode').textContent = explode > 0.5 ? 'Složit' : 'Rozložit na díly'; applyExplode(); };
document.querySelectorAll('[data-view]').forEach((b) => b.onclick = () => goView(b.dataset.view));
$('#v3dZoomIn').onclick = () => zoomBy(0.75); $('#v3dZoomOut').onclick = () => zoomBy(1.33);
function zoomBy(f) { const d = camera.position.clone().sub(controls.target).multiplyScalar(f); camera.position.copy(controls.target).add(d); setAuto(false); }
const CUTS = { zadny: null, podelny: [new THREE.Vector3(0, 0, -1), -60, 60, 0], pricny: [new THREE.Vector3(-1, 0, 0), 0, 1330, 560], vodorovny: [new THREE.Vector3(0, -1, 0), 0, 200, 60] };
let cutKind = 'zadny';
function setCut(k) {
  cutKind = k; document.querySelectorAll('[data-cut]').forEach((b) => b.classList.toggle('on', b.dataset.cut === k));
  const s = $('#v3dCutSlider'); const c = CUTS[k];
  if (!c) { clip.constant = 1e5; clip.normal.set(0, 0, -1); s.disabled = true; return; }
  s.disabled = false; s.min = c[1]; s.max = c[2]; s.value = c[3]; clip.normal.copy(c[0]); clip.constant = +c[3];
}
document.querySelectorAll('[data-cut]').forEach((b) => b.onclick = () => setCut(b.dataset.cut));
$('#v3dCutSlider').oninput = (e) => { if (CUTS[cutKind]) clip.constant = +e.target.value; };
$('#v3dAll').onclick = () => { parts.forEach((p) => p.mesh.visible = true); syncList(); };
setCut('zadny'); setAuto(true);

// výběr kliknutím (bez tažení)
let down = null;
canvas.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
canvas.addEventListener('pointerup', (e) => {
  if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 5) return;
  const r = canvas.getBoundingClientRect(); pick.setFromCamera(new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1), camera);
  const hits = pick.intersectObjects(parts.filter((p) => p.mesh.visible).map((p) => p.mesh), false).filter((h) => clip.distanceToPoint(h.point) >= -0.01 || cutKind === 'zadny');
  select(hits.length ? parts.find((p) => p.mesh === hits[0].object) : null);
});
function select(p) {
  if (selected) selected.mesh.material.emissive?.setHex(0x000000);
  selected = p; const card = $('#v3dInfo');
  if (!p) { card.hidden = true; return; }
  p.mesh.material.emissive?.setHex(0x553300);
  card.hidden = false; card.querySelector('h3').textContent = p.info.label;
  card.querySelector('.st').textContent = p.info.state; card.querySelector('.st').dataset.s = p.info.state;
  card.querySelector('p').textContent = p.info.note;
  card.querySelector('.grp').textContent = GROUPS.find((g) => g[0] === p.info.group)?.[1] || '';
}
$('#v3dHide').onclick = () => { if (selected) { selected.mesh.visible = false; select(null); syncList(); } };
$('#v3dOnly').onclick = () => { if (selected) { parts.forEach((p) => p.mesh.visible = p === selected || p.name === selected.name); syncList(); } };

// seznam dílů
function buildList() {
  const host = $('#v3dList'); host.innerHTML = '';
  for (const [g, label] of GROUPS) {
    const items = parts.filter((p) => p.info.group === g); if (!items.length) continue;
    const det = document.createElement('details'); det.open = g === 'poklop' || g === 'pohon';
    const sum = document.createElement('summary'); const gcb = document.createElement('input'); gcb.type = 'checkbox'; gcb.checked = true; gcb.dataset.group = g;
    gcb.onclick = (e) => { e.stopPropagation(); items.forEach((p) => p.mesh.visible = gcb.checked); syncList(); };
    sum.append(gcb, ' ', label, Object.assign(document.createElement('span'), { className: 'cnt', textContent: ` ${items.length}` })); det.append(sum);
    const seen = new Map(); for (const p of items) { if (!seen.has(p.info.label)) seen.set(p.info.label, []); seen.get(p.info.label).push(p); }
    for (const [lab, ps] of seen) {
      const row = document.createElement('label'); row.className = 'row';
      const cb = document.createElement('input'); cb.type = 'checkbox'; cb.checked = true; cb.dataset.label = lab;
      cb.onchange = () => { ps.forEach((p) => p.mesh.visible = cb.checked); syncList(); };
      const only = document.createElement('button'); only.type = 'button'; only.textContent = 'jen tento'; only.className = 'only';
      only.onclick = (e) => { e.preventDefault(); parts.forEach((p) => p.mesh.visible = ps.includes(p)); select(ps[0]); syncList(); goFit(ps); };
      row.append(cb, ' ', lab + (ps.length > 1 ? ` ×${ps.length}` : ''), only); det.append(row);
    }
    host.append(det);
  }
}
function syncList() {
  document.querySelectorAll('#v3dList input[data-label]').forEach((cb) => { const ps = parts.filter((p) => p.info.label === cb.dataset.label); cb.checked = ps.some((p) => p.mesh.visible); });
  document.querySelectorAll('#v3dList input[data-group]').forEach((cb) => { const ps = parts.filter((p) => p.info.group === cb.dataset.group); const n = ps.filter((p) => p.mesh.visible).length; cb.checked = n > 0; cb.indeterminate = n > 0 && n < ps.length; });
}
function goFit(ps) {
  const b = new THREE.Box3(); ps.forEach((p) => b.expandByObject(p.mesh)); const c = b.getCenter(new THREE.Vector3()); const r = Math.max(30, b.getSize(new THREE.Vector3()).length());
  const dir = camera.position.clone().sub(controls.target).normalize();
  camAnim = { p0: camera.position.clone(), t0: controls.target.clone(), p1: c.clone().add(dir.multiplyScalar(r * 1.6)), t1: c, s: 0 }; setAuto(false);
}

// ---------- smyčka
function resize() { const w = box.clientWidth, h = canvas.clientHeight; renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix(); }
new ResizeObserver(resize).observe(box);
renderer.setAnimationLoop(() => {
  if (Math.abs(explodeTarget - explode) > 1e-3) { explode += (explodeTarget - explode) * 0.08; $('#v3dExSlider').value = Math.round(explode * 100); applyExplode(); }
  if (camAnim) {
    camAnim.s = Math.min(1, camAnim.s + 0.035); const k = 1 - Math.pow(1 - camAnim.s, 3);
    camera.position.lerpVectors(camAnim.p0, camAnim.p1, k); controls.target.lerpVectors(camAnim.t0, camAnim.t1, k); if (camAnim.s >= 1) camAnim = null;
  }
  controls.update(); renderer.render(scene, camera);
});
