'use strict';
const test=require('node:test'),assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
let pending;
const context={setTimeout:fn=>{pending=fn;return 1;},clearTimeout:()=>{pending=null;}};
const source=file=>fs.readFileSync(path.join(__dirname,'..',file),'utf8');
vm.runInNewContext(source('boards/OpenFIREshared.js')+';this.shared=OpenFIREshared;',context);
for(const file of['js/core/lensfit.js','js/core/emitter-layout.js'])vm.runInNewContext(source(file),context);
context.OF.Boards={shared:context.shared};
vm.runInNewContext(source('js/app/emitter-alignment.js'),context);
const records=context.shared.cameraOptics,L=context.OF.EmitterLayout;
function alignment(record=records[0],custom=null){
 const a=new context.OF.EmitterAlignment({closed:false,render(){}},{cameraModel:record.model,lensValues:custom});
 a.fields.distance={number:{min:.5,max:L.MAX_DISTANCE,step:L.DISTANCE_STEP,value:2},slider:{max:L.MAX_DISTANCE,value:2}};
 a.fields.lens={number:{value:1},slider:{value:1},preset:{}};
 a.ideal={disabled:true};
 a.update=()=>{};
 return a;
}
function flush(){assert.ok(pending,'recalculation scheduled');const fn=pending;pending=null;fn();}

test('setup edits preserve the distance and the ideal action chooses the first green 10 cm step',()=>{
 for(const record of records)for(const factor of record.addOnLens?[1,.45,.67]:[1])for(const inches of[27,55,120])for(const aspect of[16/9,4/3,32/9]){
  const a=alignment(record);Object.assign(a.values,{inches,aspect,lens:factor,distance:9});a.schedule();flush();
  assert.equal(a.values.distance,9);assert.equal(a.result.settings.distance,9);
  a.schedule(true);flush();
  assert.ok(a.ideal.disabled,'the ideal action is disabled once its distance is selected');
  const d=a.values.distance;assert.ok(d>=.5&&d<=10);
  assert.equal(Number(a.fields.distance.number.value),d);assert.equal(Number(a.fields.distance.slider.value),d);
  assert.ok(Number(a.fields.distance.slider.max)>=d);
  if(L.calculate({...a.values,distance:10},record,factor,null,false).status==='good'){
   assert.equal(a.result.status,'good');assert.ok(L.fits(a.result.settings,a.result.optics,a.result.ratio,L.poses(a.result.settings,d*.9)));
   if(d>.5)assert.notEqual(L.calculate({...a.values,distance:(Math.round(d*10)-1)/10},record,factor,null,false).status,'good',`${record.label} ${factor} ${inches} ${aspect}: the step before ${d} is not yet green`);
  }else{assert.equal(d,10);assert.notEqual(a.result.status,'good');}
 }
});

test('manual distance edits remain unchanged across setup changes',()=>{
 const a=alignment();a.schedule(true);flush();
 assert.ok(a.ideal.disabled);
 a.values.distance=4;a.schedule();flush();assert.equal(a.result.settings.distance,4);
 assert.equal(a.ideal.disabled,false);
 a.values.inches=27;a.schedule();flush();assert.equal(a.values.distance,4);assert.equal(a.result.settings.distance,4);
 a.schedule(true);flush();assert.equal(a.result.status,'good');assert.ok(a.values.distance<4);assert.ok(a.ideal.disabled);
});

test('a manual distance edit cancels an earlier pending ideal-distance action',()=>{
 const a=alignment();a.values.inches=27;a.schedule(true);a.values.distance=1.5;a.schedule();flush();
 assert.equal(a.values.distance,1.5);assert.equal(a.result.settings.distance,1.5);
});

test('camera changes use the current optics and keep the wide-camera converter disabled',()=>{
 const a=alignment();a.schedule(true);flush();const before=a.values.distance;
 a.record=records[2];a.values.lens=.45;a.schedule();flush();
 assert.equal(a.factor,1);assert.equal(a.values.distance,before);
 assert.equal(a.ideal.disabled,false,'new camera re-enables the action for its different ideal distance');
 assert.ok(Object.values(a.fields.lens).every(f=>f.disabled));
 a.schedule(true);flush();assert.equal(a.result.status,'good');assert.ok(a.values.distance<before);
});

test('the ideal distance accounts for saved lens refinement',()=>{
 const a=alignment(records[1],{focal:records[1].focal*1.5});a.schedule(true);flush();
 const b=alignment(records[1]);b.schedule(true);flush();
 assert.ok(a.result.optics.refined);assert.ok(a.values.distance>b.values.distance);
});

test('a setup needing more than ten metres stops at ten with an honest status',()=>{
 const a=alignment();Object.assign(a.values,{inches:120,aspect:5/4});a.schedule(true);flush();
 assert.equal(a.values.distance,10);assert.notEqual(a.result.status,'good');
 assert.ok(a.ideal.disabled,'the capped distance is already selected, even without green');
 assert.equal(Number(a.fields.distance.slider.max),10);
 assert.equal(Number(a.fields.distance.slider.value),Number(a.fields.distance.number.value));
 a.values.distance=9.9;a.schedule();flush();assert.equal(a.ideal.disabled,false);
 a.values.distance=10;a.schedule();flush();assert.ok(a.ideal.disabled);
});

test('unreachable green advice stays honest at the supported distance limit',()=>{
 const a=alignment(records[0],{focal:records[0].focal*2});Object.assign(a.values,{inches:120,aspect:5/4});a.schedule(true);flush();
 assert.equal(a.values.distance,10);assert.notEqual(a.result.status,'good');assert.equal(a.result.minDistance,null);
});

test('switching Square and Diamond keeps manual distance until the ideal action',()=>{
 const a=alignment();a.schedule(true);flush();const square=a.values.distance;
 a.values.distance=1;a.values.diamond=true;a.schedule();flush();
 assert.ok(a.result.settings.diamond);assert.equal(a.values.distance,1);
 a.schedule(true);flush();assert.notEqual(a.values.distance,1);assert.equal(a.result.status,'good');
 a.values.distance=1;a.values.diamond=false;a.schedule();flush();assert.equal(a.values.distance,1);
 a.schedule(true);flush();assert.equal(a.values.distance,square);
});
