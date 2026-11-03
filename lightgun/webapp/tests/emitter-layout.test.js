'use strict';
const test=require('node:test'),assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const context={};
vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../boards/OpenFIREshared.js'),'utf8')+';this.shared=OpenFIREshared;',context);
vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../js/core/lensfit.js'),'utf8'),context);
vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../js/core/emitter-layout.js'),'utf8'),context);
const L=context.OF.EmitterLayout,records=context.shared.cameraOptics;
// Shapes the firmware tells apart reliably: upright up to 0.85, wide from 1.3 (corners) or 1.5.
const nonSquare=r=>r<=.85+1e-12||r>=1.3-1e-12;
test('manufacturer fields and firmware optical coefficients are generated together',()=>{
 assert.equal(records.length,3);assert.equal(records[0].vertical,23);assert.equal(records[0].horizontal,33);
 assert.equal(records[2].horizontal,111.3);assert.equal(records[2].k1,-.128946);assert.equal(records[2].sensorX,98);
});
test('physical dimensions use diagonal AND aspect ratio',()=>{
 const d=L.dimensions(55,16/9);assert.ok(Math.abs(d.width-1.2176)<.001);assert.ok(Math.abs(d.height-.6849)<.001);
 const b=L.dimensions(55,4/3);assert.ok(b.height>d.height);assert.ok(Math.abs(Math.hypot(b.width,b.height)-55*.0254)<1e-10);
});
test('narrow field cannot promise four LEDs at a short distance',()=>{
 const x=L.calculate({inches:55,distance:1,aspect:16/9},records[0]);assert.equal(x.status,'close');assert.ok(x.minDistance>3);
});
test('wide field can recommend a horizontal rectangle without a camera-ID branch',()=>{
 const camera={...records[2],model:12345};const x=L.calculate({inches:55,distance:2},camera);
 assert.equal(x.status,'good');assert.ok(x.ratio>1.2);assert.ok(x.base<=x.width);assert.ok(x.minDistance<2);
});
test('Diamond geometry stays at the centres of its four edges',()=>{
 const a=L.calculate({diamond:true,inches:55,distance:.6},records[2]),b=L.calculate({diamond:true,inches:55,distance:3},records[2]);
 assert.equal(a.ratio,b.ratio);assert.equal(a.base,a.width);assert.equal(b.base,b.width);
 const p=L.leds(a.settings,a.ratio);assert.equal(p[0][0],0);assert.equal(p[1][1],0);assert.equal(p[2][0],0);assert.equal(p[3][1],0);
});
test('one play profile: old playing space inputs are ignored',()=>{
 const a=L.calculate({inches:55,distance:3},records[1]);
 const b=L.calculate({inches:55,distance:3,scenario:'wide',separation:1.5},records[1]);
 assert.equal(JSON.stringify([a.ratio,a.minDistance,a.status]),JSON.stringify([b.ratio,b.minDistance,b.status]));
 assert.ok(!('scenario' in a.settings)&&!('separation' in a.settings));
});
test('add-on factor changes tangent of FOV, not degrees by linear division',()=>{
 const x=L.optics(records[0],.45);assert.ok(x.horizontal>66&&x.horizontal<68);assert.ok(x.estimated);
 assert.equal(L.optics(records[2]).estimated,false);
});
test('refinement is local, never overwrites the nominal record',()=>{
 const original=JSON.stringify(records[2]);const x=L.optics(records[2],1,{focal:records[2].focal*1.1});
 assert.ok(x.horizontal<records[2].horizontal);assert.equal(JSON.stringify(records[2]),original);assert.ok(x.refined);
});
test('refined radial bounds widen as in firmware; invalid curves and focal values are ignored',()=>{
 const a=L.optics(records[2],1,{k1:0,k2:.25,focal:records[2].focal*.1});
 assert.equal(a.focal,records[2].focal);assert.ok(a.max>=2);assert.equal(a.k2,.25);
 const b=L.optics(records[2],1,{k1:-2,k2:-2});assert.equal(b.k1,records[2].k1);assert.equal(b.k2,records[2].k2);
});
test('invalid inputs are bounded and proposed layouts are clearly non-square',()=>{
 const a=L.calculate({inches:NaN,distance:-5,aspect:Infinity,scenario:'bad'},records[2]);
 assert.ok(a.settings.distance>=.3);assert.ok(nonSquare(a.ratio));assert.ok(Number.isFinite(a.base));
});
test('position quotes are symmetric for every supported camera and monitor size',()=>{
 for(const r of records)for(const inches of [20,27,55,85,120]){
  const a=L.calculate({inches,distance:3},r,1,null,false);
  assert.ok(a.base>0&&a.base<=a.width+1e-10);assert.ok(a.ratio<=a.settings.aspect);
  const p=L.leds(a.settings,a.ratio);assert.equal(p[0][0],-p[1][0]);assert.equal(p[2][1],-p[0][1]);
 }
});
test('direct radial boundary test agrees with projected coverage',()=>{
 for(const record of records)for(const ratio of [.5,.711,1.3,16/9])for(const distance of [.5,.8,1.2,2,3,5,8]){
  const s=L.settings({inches:55,distance}),o=L.optics(record),p=L.poses(s,distance);
  const c=L.coverage(s,o,ratio,p);
  if(Math.abs(c.margin-.05)>1e-4)assert.equal(L.fits(s,o,ratio,p),c.margin>=.05);
 }
});
test('projection reaches the declared vertical FOV boundary for every native lens',()=>{
 for(const record of records){
  const o=L.optics(record),distance=1;
  const pose={eye:[0,0,-distance],f:[0,0,1],right:[1,0,0],down:[0,1,0]};
  const point=L.project(o,pose,[0,distance*o.ty,0]);
  assert.ok(Math.abs(point[0])<1e-10);assert.ok(Math.abs(point[1]-1)<2e-5);
 }
});
test('all preference formats and ultrawide monitor formats use their actual physical proportions',()=>{
 for(const aspect of [16/9,16/10,3/2,5/4,4/3,21/9,21/10,32/9]){
  const r=L.calculate({inches:55,distance:3,aspect},records[2],1,null,false);
  assert.ok(Math.abs(r.width/r.height-aspect)<1e-12);
  assert.ok(r.base<=r.width+1e-10);assert.ok(nonSquare(r.ratio));
 }
});
test('candidates avoid near-square shapes: upright up to 0.85, wide from 1.5 or at the corners',()=>{
 for(const aspect of [16/9,16/10,3/2,5/4,4/3,21/9,32/9]){
  const c=L.candidates(L.settings({aspect}));
  assert.ok(c.includes(.85));assert.ok(c.every(r=>r<=.85||r>=1.5||(r===c[c.length-1]&&r>=1.3)));
  assert.ok(c.every(r=>r>=L.DEFAULT_RATIO),'never narrower than the proportion recommended so far');
  assert.ok(c.every(r=>r<=aspect+1e-12));
 }
 assert.ok(!L.candidates(L.settings({aspect:5/4})).includes(5/4));
 // 1.3 and 1.4 only when even one step forward the player stays 4 screen heights away.
 const s=L.settings({inches:55,aspect:16/9}),h=L.dimensions(55,16/9).height,edge=L.FAR_HEIGHTS*h/L.FAR_STEP;
 for(const d of [edge*.99,edge*1.01]){
  const c=L.candidates(s,d),far=d>=edge;
  assert.equal(c.includes(1.3),far);assert.equal(c.includes(1.4),far);
  assert.ok(c.every(r=>r<=.85||r>=1.3));
 }
 // A 5:4 screen (1.25) never gets its corners, not even from far away.
 const q=L.settings({inches:55,aspect:5/4});
 assert.ok(!L.candidates(q,20).includes(5/4));
 for(const r of records)assert.ok(nonSquare(L.calculate({inches:55,aspect:5/4,distance:5},r,1,null,false).ratio));
 assert.equal(JSON.stringify(L.candidates(L.settings({aspect:16/9,diamond:true}))),JSON.stringify([16/9]));
});
test('the suggested minimum distance uses the same step-forward margin as the layout choice',()=>{
 for(const r of records)for(const inches of [32,55,85])for(const distance of [1.5,2.5,3.5,5,7]){
  const x=L.calculate({inches,distance},r);
  if(x.minDistance===null)continue;
  if(x.status==='good')assert.ok(x.minDistance<=distance+.05,`${r.label} ${inches} ${distance}: ${x.minDistance}`);
  if(x.moveTo!==null)assert.ok(x.minDistance<=x.moveTo+.05,`${r.label} ${inches} ${distance}: ${x.minDistance} > ${x.moveTo}`);
  const there=L.calculate({inches,distance:x.minDistance},r,1,null,false);
  assert.notEqual(there.status,'close',`${r.label} ${inches} ${distance}: close at ${x.minDistance}`);
 }
});
test('lean-forward reserve is 10%, while near-square wide shapes keep the 15% reserve',()=>{
 assert.equal(L.STEP_FORWARD,.9);assert.equal(L.FAR_STEP,.85);
});
test('moving away, the advice never goes back to a layout it already left',()=>{
 for(const r of records)for(const inches of [27,55]){
  const seen=[];
  for(let d=1;d<=7.001;d+=.05){const x=L.calculate({inches,distance:d},r,1,null,false);
   if(seen[seen.length-1]!==x.ratio){assert.ok(!seen.includes(x.ratio),`${r.label} ${inches} ${d.toFixed(2)}: back to ${x.ratio}`);seen.push(x.ratio);}}
 }
});
test('too close: the advice does not change while moving back, and the suggested distance works',()=>{
 for(const r of records)for(const inches of [32,55,85]){
  const near=[.5,.8,1.1].map(distance=>L.calculate({inches,distance},r,1,null,false)).filter(x=>x.status==='close'&&x.moveTo!==null);
  for(const x of near){
   assert.equal(x.ratio,near[0].ratio);assert.equal(x.moveTo,near[0].moveTo);
   const there=L.calculate({inches,distance:x.moveTo},r,1,null,false);
   assert.notEqual(there.status,'close');assert.equal(there.ratio,x.ratio);
  }
 }
});
test('play profile: gun straight needs four LEDs, rolled 20 or 45 degrees may lose one',()=>{
 const s=L.settings({inches:55,distance:3}),p=L.poses(s,3);
 assert.ok(p.every(q=>q.spare===0||q.spare===1));
 assert.equal(p.filter(q=>q.spare===0).length*4,p.filter(q=>q.spare===1).length);
 assert.ok(p.filter(q=>q.spare===0).every(q=>Math.abs(q.right[1])<1e-12),'straight poses have no roll');
});
