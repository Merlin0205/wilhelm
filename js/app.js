// Wilhelm – 3D prezentace (Three.js r160, GSAP ScrollTrigger). Souřadnice modelu: x délka od zádi, y pravobok, z nahoru (mm).
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { STLLoader } from 'three/addons/loaders/STLLoader.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';

const ST = window.WILHELM_STATS;                 // z v2_stats.json (vloženo do HTML)
const L = ST.length_mm, WL = ST.waterline_mm, DECK = ST.deck_mm;
const CENTER = new THREE.Vector3(L / 2, 0, 16);   // střed otáčení v souřadnicích modelu
const reduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
const isTouch = window.matchMedia('(pointer: coarse)').matches;

// model (x,y,z) -> three (x, z, -y), posunuto do středu
const toThree = (x, y, z) => new THREE.Vector3(x - CENTER.x, z - CENTER.z, -(y - CENTER.y));
function fixGeom(g) { g.rotateX(-Math.PI / 2); g.translate(-CENTER.x, -CENTER.z, CENTER.y); g.computeVertexNormals(); return g; }

// ---------- načtení geometrií (sdílené mezi scénami)
const loader = new STLLoader();
const FILES = ['hull_stern', 'hull_bow', 'tower', 'gun', 'flak', 'lid', 'planes', 'rudders', 'cup', 'siphon', 'void'];
const geoms = {};
await Promise.all(FILES.map(n => new Promise((res, rej) => loader.load(`models/${n}.stl`, g => { geoms[n] = fixGeom(g); res(); }, undefined, rej))))
  .catch(e => { document.querySelectorAll('.overlay').forEach(o => o.textContent = 'Model se nepodařilo načíst. Stránku je nutné otevřít přes http server, ne jako soubor.'); throw e; });
document.querySelectorAll('.overlay').forEach(o => o.remove());

const COLORS = { hull: 0x8a98a4, tower: 0x66737f, brass: 0xc9a24a, cup: 0x9fd3ff, siphon: 0xffa94d, void: 0x7fb3e0, lid: 0xd3b35a };
const PART_INFO = {
  hull_stern: { name: 'Trup, záďový díl', color: COLORS.hull, explode: [0, 0, 0], file: 'hull_stern_print.stl' },
  hull_bow:   { name: 'Trup, příďový díl', color: COLORS.hull, explode: [46, 0, 0], file: 'hull_bow_print.stl' },
  tower:      { name: 'Věž', color: COLORS.tower, explode: [0, 0, 38], file: 'tower.stl' },
  gun:        { name: 'Dělo 8,8 cm', color: COLORS.brass, explode: [0, 0, 28], file: 'gun_88.stl' },
  flak:       { name: 'Flak 2 cm', color: COLORS.brass, explode: [0, 0, 22], file: 'flak_20.stl' },
  lid:        { name: 'Zátka s mřížkou', color: COLORS.lid, explode: [0, 0, -30], file: 'lid_normal.stl' },
  planes:     { name: 'Hloubková kormidla (otočná)', color: COLORS.brass, explode: [0, 0, -24], file: 'plane_bow.stl, plane_stern.stl' },
  rudders:    { name: 'Kormidla (otočná)', color: COLORS.brass, explode: [0, 0, -30], file: 'rudder.stl ×2' },
};
const INNER = { cup: COLORS.cup, siphon: COLORS.siphon, void: COLORS.void };

function makeEnv(renderer) {
  const pmrem = new THREE.PMREMGenerator(renderer);
  const env = pmrem.fromScene(new RoomEnvironment(renderer), 0.04).texture; pmrem.dispose(); return env;
}
function buildSub(clipPlane) {
  const root = new THREE.Group(); const parts = {}; const inner = {};
  for (const [k, info] of Object.entries(PART_INFO)) {
    const m = new THREE.MeshStandardMaterial({ color: info.color, metalness: k === 'gun' || k === 'flak' || k === 'planes' || k === 'lid' ? 0.75 : 0.35, roughness: k === 'lid' ? 0.35 : 0.5, side: THREE.DoubleSide, clippingPlanes: [clipPlane], clipShadows: true });
    const mesh = new THREE.Mesh(geoms[k], m); mesh.userData = { key: k, explode: new THREE.Vector3(info.explode[0], info.explode[2], -info.explode[1]), home: new THREE.Vector3() };
    mesh.castShadow = mesh.receiveShadow = true; root.add(mesh); parts[k] = mesh;
  }
  for (const [k, c] of Object.entries(INNER)) {
    const m = new THREE.MeshStandardMaterial({ color: c, emissive: c, emissiveIntensity: 0.25, transparent: true, opacity: 0, depthWrite: false, side: THREE.DoubleSide, clippingPlanes: [clipPlane] });
    const mesh = new THREE.Mesh(geoms[k], m); mesh.visible = false; mesh.renderOrder = 2; root.add(mesh); inner[k] = mesh;
  }
  return { root, parts, inner };
}
function setExplode(sub, f) { for (const m of Object.values(sub.parts)) m.position.copy(m.userData.explode).multiplyScalar(f); }
function setInner(sub, on, fade) {
  for (const m of Object.values(sub.inner)) { m.visible = on || fade > 0; m.material.opacity = Math.max(0, Math.min(0.55, fade)); }
}
function addLights(scene) {
  scene.add(new THREE.HemisphereLight(0xcfe0f0, 0x101820, 0.5));
  const key = new THREE.DirectionalLight(0xffffff, 1.6); key.position.set(180, 260, 220); key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048); key.shadow.camera.left = -260; key.shadow.camera.right = 260; key.shadow.camera.top = 160; key.shadow.camera.bottom = -160; key.shadow.camera.near = 50; key.shadow.camera.far = 900; key.shadow.bias = -0.0005;
  scene.add(key);
  const rim = new THREE.DirectionalLight(0x9fc3e6, 0.6); rim.position.set(-260, 80, -200); scene.add(rim);
  return key;
}
function makeRenderer(canvas) {
  const r = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true, powerPreference: 'high-performance' });
  r.setPixelRatio(Math.min(window.devicePixelRatio, 2)); r.shadowMap.enabled = true; r.shadowMap.type = THREE.PCFSoftShadowMap;
  r.toneMapping = THREE.ACESFilmicToneMapping; r.toneMappingExposure = 1.05; r.localClippingEnabled = true; return r;
}
// ---------- voda a bubliny (sdílené pomocné tvary)
const WATER_VS = `uniform float t; varying vec3 vPos; varying vec3 vN;
  void main(){ vec3 p = position; float w = sin(p.x*0.05 + t*1.3)*1.6 + sin(p.y*0.07 - t*0.9)*1.1 + sin((p.x+p.y)*0.03 + t*0.6)*1.4; p.z += w;
  vPos = (modelMatrix*vec4(p,1.0)).xyz; vN = normalize(vec3(-cos(p.x*0.05+t*1.3)*0.08, -cos(p.y*0.07-t*0.9)*0.077, 1.0));
  gl_Position = projectionMatrix*modelViewMatrix*vec4(p,1.0);}`;
const WATER_FS = `uniform vec3 cShallow; uniform vec3 cDeep; uniform float alpha; varying vec3 vPos; varying vec3 vN;
  void main(){ vec3 V = normalize(cameraPosition - vPos); float f = pow(1.0 - max(dot(V, normalize(vN)), 0.0), 2.2);
  vec3 c = mix(cDeep, cShallow, f*0.9 + 0.1); gl_FragColor = vec4(c, alpha + f*0.3);}`;
function makeWater() {
  const g = new THREE.PlaneGeometry(1400, 900, 90, 60);
  const m = new THREE.ShaderMaterial({ vertexShader: WATER_VS, fragmentShader: WATER_FS, transparent: true, depthWrite: false, side: THREE.DoubleSide,
    uniforms: { t: { value: 0 }, cShallow: { value: new THREE.Color(0x6fb6ff) }, cDeep: { value: new THREE.Color(0x0f3f6e) }, alpha: { value: 0.42 } } });
  const mesh = new THREE.Mesh(g, m); mesh.rotation.x = -Math.PI / 2; mesh.renderOrder = 3; return mesh;
}
function bubbleTexture() {
  const c = document.createElement('canvas'); c.width = c.height = 64; const x = c.getContext('2d');
  const gr = x.createRadialGradient(32, 32, 4, 32, 32, 30); gr.addColorStop(0, 'rgba(255,255,255,.95)'); gr.addColorStop(0.35, 'rgba(190,230,255,.55)'); gr.addColorStop(1, 'rgba(160,215,255,0)');
  x.fillStyle = gr; x.fillRect(0, 0, 64, 64); const t = new THREE.CanvasTexture(c); return t;
}
class Bubbles {
  constructor(n = 400) {
    this.n = n; this.pos = new Float32Array(n * 3); this.vel = new Float32Array(n); this.life = new Float32Array(n); this.i = 0;
    const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.BufferAttribute(this.pos, 3));
    this.points = new THREE.Points(g, new THREE.PointsMaterial({ size: 3.2, map: bubbleTexture(), transparent: true, depthWrite: false, opacity: 0.9, sizeAttenuation: true }));
    this.points.frustumCulled = false; this.points.renderOrder = 4; for (let k = 0; k < n; k++) this.pos[k * 3 + 1] = -9999;
  }
  spawn(p, spread, speed) { const k = this.i = (this.i + 1) % this.n; this.pos[k * 3] = p.x + (Math.random() - .5) * spread; this.pos[k * 3 + 1] = p.y; this.pos[k * 3 + 2] = p.z + (Math.random() - .5) * spread; this.vel[k] = speed * (0.7 + Math.random() * 0.6); this.life[k] = 1; }
  update(dt, surfaceY) {
    for (let k = 0; k < this.n; k++) { if (this.life[k] <= 0) continue; this.pos[k * 3 + 1] += this.vel[k] * dt; this.pos[k * 3] += Math.sin(this.pos[k * 3 + 1] * 0.3 + k) * 8 * dt; if (this.pos[k * 3 + 1] > surfaceY) { this.life[k] = 0; this.pos[k * 3 + 1] = -9999; } }
    this.points.geometry.attributes.position.needsUpdate = true;
  }
}

// ---------- fyzika cyklu (společná pro scénu ponoření)
class Dive {
  constructor() { this.reset(); this.speed = 1; this.playing = true; this.status = ''; }
  reset() { this.gas = 0; this.phase = 'fill'; this.y = -1; this.v = 0; this.dumpT = 0; }
  step(dt, depth) {
    const NEUTRAL = Math.abs(ST.net_g) / ST.gas_capacity_cm3 + 0.1;   // podíl plynu pro neutrální vztlak (+prášek)
    if (this.phase === 'fill') { this.gas = Math.min(1, this.gas + dt / 9); if (this.gas >= 1) { this.phase = 'dump'; this.dumpT = 0; } }
    else { this.dumpT += dt; this.gas = Math.max(0, 1 - this.dumpT / 1.4); if (this.gas <= 0) this.phase = 'fill'; }
    const lift = (this.gas - NEUTRAL) * 70; this.v += lift * dt; this.v *= Math.pow(0.3, dt); this.y += this.v * dt;
    if (this.y > 0) { this.y = 0; this.v = Math.min(this.v, 0); } if (this.y < -depth) { this.y = -depth; this.v = Math.max(this.v, 0); }
    const atTop = this.y > -0.5, atBottom = this.y < -depth + 0.5;
    this.status = this.phase === 'dump' ? 'sifon vypouští plyn, ponorka klesá' : atTop ? 'na hladině, komora se doplňuje' : this.v > 3 ? 'stoupá' : this.v < -3 ? 'klesá' : atBottom ? (this.gas < NEUTRAL ? 'leží na dně, hromadí plyn' : 'odlepuje se ode dna') : 'vznáší se';
  }
}

// =====================================================================  STORY (sticky plátno)
(function story() {
  const canvas = document.getElementById('storyCanvas'); if (!canvas) return;
  const stage = canvas.parentElement;
  const renderer = makeRenderer(canvas);
  const scene = new THREE.Scene(); scene.environment = makeEnv(renderer); scene.environmentIntensity = 0.9;
  const camera = new THREE.PerspectiveCamera(30, 16 / 9, 1, 4000);
  const clip = new THREE.Plane(new THREE.Vector3(0, 0, -1), 1e6);
  const sub = buildSub(clip); scene.add(sub.root);
  addLights(scene);
  // podlaha pro stín (jen v suchých scénách)
  const floor = new THREE.Mesh(new THREE.PlaneGeometry(1400, 900), new THREE.ShadowMaterial({ opacity: 0.35 })); floor.rotation.x = -Math.PI / 2; floor.position.y = -CENTER.z - 2; floor.receiveShadow = true; scene.add(floor);
  // voda, dno, bubliny (scéna ponoření)
  const water = makeWater(); water.visible = false; scene.add(water);
  const DEPTH = 150;
  const tank = new THREE.Mesh(new THREE.PlaneGeometry(1400, 900), new THREE.MeshStandardMaterial({ color: 0x1b2a36, roughness: 1 })); tank.rotation.x = -Math.PI / 2; tank.visible = false; tank.receiveShadow = true; scene.add(tank);
  const bubbles = new Bubbles(500); bubbles.points.visible = false; scene.add(bubbles.points);
  // plyn v komoře
  const cupH = ST.cup.top - 2.4, gas = new THREE.Mesh(new THREE.BoxGeometry(ST.cup.l - 0.5, 1, ST.cup.w - 0.5), new THREE.MeshStandardMaterial({ color: 0xbfe6ff, emissive: 0x9fd3ff, emissiveIntensity: 0.6, transparent: true, opacity: 0.75, clippingPlanes: [clip] }));
  gas.visible = false; sub.root.add(gas);
  const siphonMat = sub.inner.siphon.material;

  const controls = new OrbitControls(camera, canvas); controls.enableZoom = false; controls.enablePan = false; controls.enableDamping = true; controls.dampingFactor = 0.07; controls.autoRotateSpeed = 0.9;
  controls.enabled = !isTouch;
  const labelsEl = document.getElementById('storyLabels'); const capEl = document.getElementById('diveCaption');
  const LABELS = {
    tower: ['věž (Turm 0)', ST.tower_cx, 0, DECK + 16], gun: ['dělo 8,8 cm SK C/35', ST.gun_x, 0, DECK + 7], flak: ['flak 2 cm', ST.flak_x, 0, DECK + 5],
    bow: ['příďový díl 118 mm', ST.split_x + 60, 0, DECK + 3], stern: ['záďový díl 215 mm', 90, 0, DECK + 3], split: ['spoj: 2 kolíky Ø4', ST.split_x, -16, 14],
    lid: ['zátka s mřížkou', ST.cup.cx, 0, -4], planes: ['otočná hloubková kormidla, západky 0 / ±10 / ±20°', 60.5 * 5, 20, 16], rud: ['dvě kormidla', 13, 0, 8],
    cup: ['komora na prášek 12 × 60 mm, 25 mm před těžištěm', ST.cup.cx, 0, 5], siphon: ['plynový sifon Ø 9 mm', ST.siphon.xb, 0, -7], voidl: ['uzavřená dutina plováku ' + ST.void_cm3 + ' cm³', ST.cup.cx - 60, 0, 26], slots: ['zaplavovací otvory nástavby', 260, 9, DECK - 2.6], rudders: ['kormidla na čepech, stabilizátor', 10, 0, 6], name: ['reliéf WILHELM', 77, -9, 20],
    screws: ['otvory M4 pro dolaďovací šroubky', ST.screw_x[1], 0, 3], exit: ['výstup plynu, pravobok', ST.siphon.xc - 14, 14, ST.siphon.exit_z + 12],
  };
  const labelEls = {}; for (const k of Object.keys(LABELS)) { const d = document.createElement('div'); d.className = 'lbl'; d.textContent = LABELS[k][0]; labelsEl.appendChild(d); labelEls[k] = d; }

  const state = { explode: 0, inner: 0, clip: 0, clipOn: 0, water: 0, dive: 0, auto: 1 };
  const camState = { px: 300, py: 150, pz: 360, tx: 0, ty: 0, tz: 0 };
  const SCENES = {
    hero:    { cam: [380, 190, 470, -50, 0, 0], explode: 0, inner: 0, clipOn: 0, water: 0, dive: 0, auto: 1, labels: [] },
    side:    { cam: [10, 40, 470, 0, 0, 0], explode: 0, inner: 0, clipOn: 0, water: 0, dive: 0, auto: 0, labels: ['tower', 'gun', 'flak', 'planes', 'rudders', 'slots', 'name'] },
    explode: { cam: [330, 190, 330, 10, 0, 0], explode: 1, inner: 0, clipOn: 0, water: 0, dive: 0, auto: 0, labels: ['bow', 'stern', 'split', 'tower', 'gun', 'flak', 'lid', 'planes', 'rudders'] },
    inside:  { cam: [80, 120, 330, 0, -6, 0], explode: 0, inner: 1, clipOn: 1, water: 0, dive: 0, auto: 0, labels: ['cup', 'siphon', 'voidl', 'exit', 'screws', 'lid'] },
    dive:    { cam: [220, 60, 420, 0, -60, 0], explode: 0, inner: 0.6, clipOn: 1, water: 1, dive: 1, auto: 0, labels: [] },
  };
  let current = 'hero';
  function setScene(name) {
    if (!SCENES[name] || current === name) return; current = name; const s = SCENES[name];
    const d = reduced ? 0 : 1.4;
    gsap.to(camState, { px: s.cam[0], py: s.cam[1], pz: s.cam[2], tx: s.cam[3], ty: s.cam[4], tz: s.cam[5], duration: d, ease: 'power2.inOut', overwrite: true });
    gsap.to(state, { explode: s.explode, inner: s.inner, clipOn: s.clipOn, water: s.water, dive: s.dive, auto: s.auto, duration: d, ease: 'power2.inOut', overwrite: true });
    for (const k of Object.keys(labelEls)) labelEls[k].classList.toggle('on', s.labels.includes(k));
    capEl.classList.toggle('on', name === 'dive');
    document.getElementById('scrollHint').style.opacity = name === 'hero' ? 1 : 0;
    if (name === 'dive') dive.reset();
  }
  const dive = new Dive();
  document.getElementById('divePlay')?.addEventListener('click', e => { dive.playing = !dive.playing; e.currentTarget.textContent = dive.playing ? 'Pauza' : 'Přehrát'; });
  document.getElementById('diveSpeed')?.addEventListener('input', e => { dive.speed = parseFloat(e.target.value); });

  function resize() { const w = stage.clientWidth, h = stage.clientHeight; if (!w || !h) return; renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix(); }
  new ResizeObserver(resize).observe(stage); resize();

  const tmp = new THREE.Vector3(); let last = performance.now(); let userRotating = false;
  controls.addEventListener('start', () => { userRotating = true; }); controls.addEventListener('end', () => { userRotating = false; });
  function frame(now) {
    requestAnimationFrame(frame); const dt = Math.min(0.05, (now - last) / 1000); last = now;
    // kamera
    const cf = camera.aspect < 0.8 ? 1.75 : camera.aspect < 1.2 ? 1.3 : 1;   // na výšku dál od modelu
    if (!userRotating) { camera.position.lerp(tmp.set(camState.px * cf, camState.py * cf, camState.pz * cf), 0.12); }
    controls.target.lerp(tmp.set(camState.tx, camState.ty, camState.tz), 0.12);
    controls.autoRotate = state.auto > 0.5 && !reduced; controls.update();
    // stav
    setExplode(sub, state.explode); setInner(sub, state.inner > 0.01, state.inner * 0.55);
    clip.constant = state.clipOn > 0.5 ? 0.0 : 1e6;
    gas.material.clippingPlanes = [clip];
    water.visible = state.water > 0.02; water.material.uniforms.alpha.value = 0.42 * state.water; water.material.uniforms.t.value += dt;
    tank.visible = state.dive > 0.02; floor.visible = state.dive < 0.5;
    // ponoření
    const subY = state.dive > 0.02 ? dive.y : 0;
    if (state.dive > 0.02 && dive.playing && !reduced) dive.step(dt * dive.speed, DEPTH);
    water.position.y = THREE.MathUtils.lerp(300, WL - CENTER.z, state.dive);
    tank.position.y = -CENTER.z - DEPTH - 2;
    sub.root.position.y = subY * state.dive;
    sub.root.rotation.z = state.dive > 0.02 ? THREE.MathUtils.degToRad(-14 + 32 * dive.gas) * state.dive : 0;   // kluzák: příď nahoru s plynem
    gas.visible = state.dive > 0.02;
    const gh = Math.max(0.01, cupH * dive.gas); gas.scale.y = gh; gas.position.set(ST.cup.cx - CENTER.x, (ST.cup.top - gh / 2) - CENTER.z, 0);
    siphonMat.emissiveIntensity = dive.phase === 'dump' ? 1.4 : 0.25; sub.inner.siphon.material.opacity = state.dive > 0.02 ? 0.85 : state.inner * 0.55;
    bubbles.points.visible = state.dive > 0.02;
    if (state.dive > 0.5 && dive.playing && !reduced) {
      if (dive.phase === 'fill' && Math.random() < 0.6) bubbles.spawn(toThree(ST.cup.cx + (Math.random() - .5) * 50, 0, 3).add(sub.root.position), 6, 25);
      if (dive.phase === 'dump') for (let k = 0; k < 3; k++) bubbles.spawn(toThree(ST.siphon.xc, 16, ST.siphon.exit_z).add(sub.root.position), 4, 70);
    }
    bubbles.update(dt, WL - CENTER.z);
    if (current === 'dive') capEl.textContent = `${dive.status} · plyn ${(dive.gas * ST.gas_capacity_cm3).toFixed(1).replace('.', ',')} cm³`;
    // popisky
    const w = stage.clientWidth, h = stage.clientHeight;
    const LABEL_PART = { tower: 'tower', gun: 'gun', flak: 'flak', lid: 'lid', planes: 'planes', rudders: 'rudders', bow: 'hull_bow' };
    for (const k of Object.keys(labelEls)) { if (!labelEls[k].classList.contains('on')) continue; const L3 = LABELS[k];
      const base = toThree(L3[1], L3[2], L3[3]); const part = sub.parts[LABEL_PART[k]]; if (part) base.add(part.position); base.add(sub.root.position); base.project(camera);
      labelEls[k].style.left = `${(base.x * 0.5 + 0.5) * w}px`; labelEls[k].style.top = `${(-base.y * 0.5 + 0.5) * h - 22}px`; }
    renderer.render(scene, camera);
  }
  requestAnimationFrame(frame);

  // ScrollTrigger kroky
  gsap.registerPlugin(ScrollTrigger);
  document.querySelectorAll('.step[data-scene]').forEach(step => {
    ScrollTrigger.create({ trigger: step, start: 'top 55%', end: 'bottom 55%', onEnter: () => setScene(step.dataset.scene), onEnterBack: () => setScene(step.dataset.scene) });
  });
})();

// =====================================================================  samostatný prohlížeč
(function viewer() {
  const canvas = document.getElementById('viewerCanvas'); if (!canvas) return;
  const stage = canvas.parentElement; const renderer = makeRenderer(canvas);
  const scene = new THREE.Scene(); scene.environment = makeEnv(renderer);
  const camera = new THREE.PerspectiveCamera(32, 16 / 9, 1, 4000); camera.position.set(260, 150, 330);
  const clip = new THREE.Plane(new THREE.Vector3(0, 0, -1), 1e6); const sub = buildSub(clip); scene.add(sub.root); addLights(scene);
  const floor = new THREE.Mesh(new THREE.PlaneGeometry(1400, 900), new THREE.ShadowMaterial({ opacity: 0.3 })); floor.rotation.x = -Math.PI / 2; floor.position.y = -CENTER.z - 2; floor.receiveShadow = true; scene.add(floor);
  const water = new THREE.Mesh(new THREE.PlaneGeometry(600, 160), new THREE.MeshBasicMaterial({ color: 0x3b8bd1, transparent: true, opacity: 0.28, side: THREE.DoubleSide, depthWrite: false })); water.rotation.x = -Math.PI / 2; water.position.y = WL - CENTER.z; water.visible = false; scene.add(water);
  const grid = new THREE.GridHelper(600, 60, 0x24405a, 0x16293a); grid.position.y = -CENTER.z - 1.9; scene.add(grid);
  const controls = new OrbitControls(camera, canvas); controls.enableDamping = true; controls.minDistance = 40; controls.maxDistance = 1200; controls.autoRotate = !reduced; controls.autoRotateSpeed = 0.8;
  controls.addEventListener('start', () => { controls.autoRotate = false; pressView(null); });
  function resize() { const w = stage.clientWidth, h = stage.clientHeight; if (!w || !h) return; renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix(); }
  new ResizeObserver(resize).observe(stage); resize();
  const views = { iso: [260, 150, 330], side: [0, 10, 520], top: [0, 520, 0.01], bottom: [0, -520, 0.01], bow: [420, 60, 120], stern: [-420, 60, 120] };
  const vbtns = document.querySelectorAll('#viewer .btn[data-view]');
  function pressView(v) { vbtns.forEach(b => b.setAttribute('aria-pressed', String(b.dataset.view === v))); }
  vbtns.forEach(b => b.addEventListener('click', () => { const p = views[b.dataset.view]; gsap.to(camera.position, { x: p[0], y: p[1], z: p[2], duration: reduced ? 0 : 0.9, ease: 'power2.inOut' }); controls.target.set(0, 0, 0); controls.autoRotate = false; pressView(b.dataset.view); }));
  const vs = { explode: 0, inner: 0 };
  function tog(id, fn) { const b = document.getElementById(id); b?.addEventListener('click', () => { const on = b.getAttribute('aria-pressed') !== 'true'; b.setAttribute('aria-pressed', String(on)); fn(on); }); }
  tog('vWater', on => { water.visible = on; });
  tog('vExplode', on => gsap.to(vs, { explode: on ? 1 : 0, duration: reduced ? 0 : 0.8, ease: 'power2.inOut' }));
  tog('vInner', on => gsap.to(vs, { inner: on ? 1 : 0, duration: reduced ? 0 : 0.6 }));
  const cutRange = document.getElementById('vCutPos');
  tog('vCut', on => { cutRange.disabled = !on; clip.constant = on ? parseFloat(cutRange.value) : 1e6; });
  cutRange?.addEventListener('input', () => { clip.constant = parseFloat(cutRange.value); });
  tog('vAuto', on => { controls.autoRotate = on; });
  // seznam dílů: zobrazit/skrýt
  const list = document.getElementById('partList');
  if (list) for (const [k, info] of Object.entries(PART_INFO)) {
    const b = document.createElement('button'); b.setAttribute('aria-pressed', 'true');
    b.innerHTML = `<div><i class="sw" style="background:#${info.color.toString(16).padStart(6, '0')}"></i>${info.name}</div><span>${info.file} · ${ST.mass_g[k] ?? ''} g</span>`;
    b.addEventListener('click', () => { const on = b.getAttribute('aria-pressed') !== 'true'; b.setAttribute('aria-pressed', String(on)); sub.parts[k].visible = on; });
    list.appendChild(b);
  }
  let last = performance.now();
  function frame(now) { requestAnimationFrame(frame); const dt = Math.min(0.05, (now - last) / 1000); last = now; setExplode(sub, vs.explode); setInner(sub, vs.inner > 0.01, vs.inner * 0.55); controls.update(); renderer.render(scene, camera); }
  requestAnimationFrame(frame);
})();

// =====================================================================  drobnosti: reveal, aktivní odkaz v menu
(function ui() {
  if (!window.gsap) return; gsap.registerPlugin(ScrollTrigger);
  document.querySelectorAll('.reveal').forEach(el => { ScrollTrigger.create({ trigger: el, start: 'top 88%', once: true, onEnter: () => el.classList.add('in') }); });
  const links = [...document.querySelectorAll('.nav a.l')];
  document.querySelectorAll('section[id]').forEach(sec => {
    ScrollTrigger.create({ trigger: sec, start: 'top 45%', end: 'bottom 45%', onToggle: self => { if (self.isActive) links.forEach(a => a.classList.toggle('active', a.getAttribute('href') === '#' + sec.id)); } });
  });
})();
