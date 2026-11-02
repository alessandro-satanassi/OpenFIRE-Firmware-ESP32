/* Fullscreen emitter placement planner. EN/IT: Browser only / solo browser.
   This UI owns no serial session and never changes camera or calibration settings. */
(function(root){
    'use strict';
    const OF=root.OF=root.OF||{},doc=root.document;
    const t=(key,...args)=>OF.i18n.t(key,...args);
    const node=(tag,cls,text)=>{const e=doc.createElement(tag);if(cls)e.className=cls;if(text)e.textContent=text;return e;};
    const cm=x=>(x*100).toLocaleString(OF.i18n.currentLang,{maximumFractionDigits:1});
    class EmitterAlignment {
        constructor(win,options={}) {
            this.win=win;this.records=OF.Boards.shared.cameraOptics||[];this.fields={};this.values={inches:55,distance:2,aspect:16/9,scenario:'normal',separation:.8,diamond:!!options.diamond};
            this.model=options.cameraModel;this.custom=options.lensValues;this.factor=1;
            this.record=this.records.find(r=>r.model===this.model)||this.records[0];
        }
        mount() {
            const overlay=this.win.overlay;overlay.classList.add('alignment-planner');
            this.header=node('header','alignment-header');
            const titles=node('div');titles.append(node('small','alignment-kicker','OPENFIRE'),node('h1','',t('IR emitter placement')));
            const actions=node('div','alignment-actions');
            this.hide=node('button','',t('Hide controls'));this.hide.type='button';this.hide.onclick=()=>{this.body.hidden=!this.body.hidden;this.hide.textContent=t(this.body.hidden?'Show controls':'Hide controls');};
            const close=node('button','',t('Close'));close.type='button';close.onclick=()=>this.win.shutdown();actions.append(this.hide,close);this.header.append(titles,actions);
            this.body=node('div','alignment-body');
            const controls=node('section','alignment-card alignment-controls');
            controls.append(node('h2','',t('Your setup')));
            const modes=node('div','alignment-layouts');
            for(const [key,label] of [['square','Square'],['diamond','Diamond']]) {
                const b=node('button','',label);b.type='button';b.dataset.layout=key;b.setAttribute('aria-pressed',(key==='diamond')===this.values.diamond?'true':'false');
                b.onclick=()=>{this.values.diamond=key==='diamond';for(const n of modes.children)n.setAttribute('aria-pressed',n===b?'true':'false');this.schedule();};modes.append(b);
            }
            controls.append(modes);
            this.range(controls,'inches','Screen diagonal','″',20,120,1,[20,24,27,32,42,43,48,50,55,65,75,77,83,85,98,100,120]);
            this.range(controls,'distance','Playing distance','m',.5,10,.05,[.5,.75,1,1.5,2,2.5,3,4,5,6,8,10]);
            // EN: Include every format in the preferences; ultrawide additions are
            // planner-only and do not add or change firmware profile enum values.
            // IT: Tutti i formati delle preferenze; gli ultrawide aggiuntivi valgono
            // solo per l'assistente, senza modificare gli enum dei profili firmware.
            const formats=[...new Set([...OF.Maps.aspectRatios,'21:9','21:10','32:9'])];
            this.select(controls,'aspect','Screen format',formats.map(label=>{const [w,h]=label.split(':').map(Number);return[w/h,label];}));
            this.select(controls,'camera','Camera',this.records.map(r=>[r.model,r.label]),this.record?.model);
            this.select(controls,'lens','Additional lens',[[1,t('None')],[.67,t('0.67x converter (estimate)')],[.45,t('0.45x converter (estimate)')]],1);
            this.select(controls,'scenario','Playing space',[['normal',t('Normal movement')],['wide',t('Wide movement (up to 75° roll)')],['two',t('Two players')]]);
            this.separation=this.range(controls,'separation','Distance between lightguns','m',.5,1.5,.05,[.6,.8,1,1.2,1.5]);this.separation.hidden=true;
            this.results=node('section','alignment-card alignment-results');
            this.results.setAttribute('aria-live','polite');this.results.setAttribute('aria-atomic','true');
            this.body.append(controls,this.results);
            this.footer=node('footer','alignment-footer',t('Measure between LED centres. After moving the emitters, calibrate again.'));
            overlay.append(this.header,this.body,this.footer);
            if(!this.record){this.results.append(node('p','',t('Camera optics are unavailable. Rebuild the WebApp.')));return;}
            this.schedule();
        }
        range(parent,key,label,unit,min,max,step,ticks) {
            const box=node('div','alignment-field'),id='alignment-'+key;
            const title=node('label','',t(label));title.htmlFor=id;
            const line=node('div','alignment-value'),number=node('input');number.type='number';number.min=min;number.max=key==='distance'?20:max;number.step=step;number.id=id;number.value=this.values[key];number.setAttribute('aria-label',t(label));
            const slider=node('input');slider.type='range';slider.min=min;slider.max=max;slider.step=step;slider.value=this.values[key];slider.setAttribute('aria-label',t(label)+' — '+t('Slider'));
            const list=node('datalist');list.id=id+'-ticks';for(const tick of ticks){const o=node('option');o.value=tick;list.append(o);}slider.setAttribute('list',list.id);
            line.append(title,number,node('span','',unit));box.append(line,slider,list);parent.append(box);this.fields[key]={number,slider};
            const update=e=>{const value=Number(e.target.value);if(!e.target.value||!Number.isFinite(value))return;this.values[key]=Math.max(min,Math.min(Number(number.max),value));if(e.target!==number)number.value=this.values[key];slider.value=this.values[key];this.schedule();};
            number.oninput=slider.oninput=update;number.onchange=()=>{number.value=this.values[key];};
            if(key==='inches'){
                const preset=node('select','alignment-preset');preset.setAttribute('aria-label',t('Common screen sizes'));
                const prompt=node('option','',t('Choose a common screen size…'));prompt.value='';preset.append(prompt);
                for(const tick of ticks){const option=node('option','',tick+'″');option.value=tick;preset.append(option);}
                preset.onchange=()=>{if(!preset.value)return;number.value=preset.value;update({target:number});preset.value='';};box.append(preset);
            }
            return box;
        }
        select(parent,key,label,entries,initial) {
            const box=node('div','alignment-select'),title=node('label','',t(label)),select=node('select');select.id='alignment-'+key;title.htmlFor=select.id;
            for(const [value,label] of entries){const option=node('option','',label);option.value=value;select.append(option);}select.value=initial===undefined?this.values[key]:initial;
            box.append(title,select);parent.append(box);this.fields[key]=select;
            select.onchange=()=>{
                if(key==='camera'){this.record=this.records.find(r=>r.model===Number(select.value));this.factor=1;this.fields.lens.value=1;}
                else if(key==='lens')this.factor=Number(select.value);
                else this.values[key]=key==='aspect'?Number(select.value):select.value;
                this.schedule();
            };
        }
        schedule() {
            if(!this.record)return;
            clearTimeout(this.timer);
            // R3 already has a very wide native field. This is metadata driven,
            // not a camera ID check; future wide records get the same behaviour.
            const wide=this.record.addOnLens===false;
            this.fields.lens.disabled=wide;
            if(wide){this.factor=1;this.fields.lens.value=1;}
            this.separation.hidden=this.values.scenario!=='two';
            this.timer=setTimeout(()=>{
                if(this.win.closed)return;
                const custom=this.record.model===this.model?this.custom:null;
                this.result=OF.EmitterLayout.calculate(this.values,this.record,this.factor,custom);
                this.update();this.win.render();
            },70);
        }
        update() {
            const r=this.result,s=r.settings;
            this.results.replaceChildren(node('h2','',t('Suggested placement')));
            const badge=node('p','alignment-status '+(r.optics.estimated&&r.status==='good'?'limited':r.status),t(r.status==='good'?(r.optics.estimated?'Estimated movement margin':'Good movement margin'):r.status==='limited'?'Limited movement margin':'Move farther back'));
            this.results.append(badge);
            const title=s.diamond?t('Diamond: centre of each screen edge'):t('Square: two emitters above and two below');
            this.results.append(node('p','alignment-placement',title));
            const grid=node('dl','alignment-measures');
            const entry=(key,value)=>grid.append(node('dt','',t(key)),node('dd','',value));
            entry('Screen size',cm(r.width)+' × '+cm(r.height)+' cm');
            if(!s.diamond){entry('Horizontal LED spacing',cm(r.base)+' cm');entry('From the screen centre','± '+cm(r.base/2)+' cm');entry('From each side edge',cm((r.width-r.base)/2)+' cm');entry('Base / height',r.ratio.toLocaleString(OF.i18n.currentLang,{maximumFractionDigits:2}));}
            entry('Vertical LED spacing',cm(r.height)+' cm');
            entry('Suggested minimum distance',r.minDistance===null?t('More than 20 m'):r.minDistance.toLocaleString(OF.i18n.currentLang,{maximumFractionDigits:2})+' m');
            if(r.moveTo!==null)entry('Move back to at least',r.moveTo.toLocaleString(OF.i18n.currentLang,{maximumFractionDigits:2})+' m');
            this.results.append(grid);
            const notes=node('div','alignment-notes');
            notes.append(node('p','',t(s.scenario==='wide'?'Includes side movement and roll up to 75°; the most demanding poses need more distance.':s.scenario==='two'?'Checked from both player positions, with small side movements and roll up to 15°.':'Includes small side movements and roll up to 15°.')));
            notes.append(node('p','',t('When you tilt the gun, one LED may briefly leave the view: the gun fills it in, with slightly less precise aim.')));
            if(r.optics.estimated)notes.append(node('p','alignment-estimate',t('Unmeasured add-on lens: approximate geometry only. Check the real field of view with the IR camera test.')));
            if(this.factor>=.5&&this.factor<1)notes.append(node('p','alignment-estimate',t('After fitting the lens, use Refine lens in the App, so the gun knows the new lens.')));
            else if(this.factor<.5)notes.append(node('p','alignment-estimate',t('This lens may be too strong for Refine lens in the App: aiming may be less precise.')));
            if(r.optics.refined)notes.append(node('p','',t('Uses the lens refinement currently selected in the App.')));
            notes.append(node('p','',t('Check all four screen corners from your playing positions in the IR camera test.')));
            notes.append(node('p','alignment-small',t('Coverage is a geometric estimate, not a guarantee of detection range. LED power, mounting and reflections also matter.')));
            this.results.append(notes);
        }
        draw(ctx) {
            const w=this.win.width,h=this.win.height;
            ctx.fillStyle='#0b1421';ctx.fillRect(0,0,w,h);
            const r=this.result;if(!r)return;
            const aspect=r.settings.aspect,screenAspect=w/h;
            const rw=screenAspect>aspect?h*aspect:w,rh=rw/aspect,x=(w-rw)/2,y=(h-rh)/2;
            this.footer.textContent=t(Math.abs(screenAspect/aspect-1)<.025?'Measure between LED centres. After moving the emitters, calibrate again.':'Scaled preview: use the centimetre measurements to position the LEDs.');
            ctx.strokeStyle='#324259';ctx.lineWidth=1;ctx.strokeRect(x+.5,y+.5,rw-1,rh-1);
            const fraction=r.base/r.width;
            const points=r.settings.diamond?[[.5,0],[1,.5],[.5,1],[0,.5]]:[[.5-fraction/2,0],[.5+fraction/2,0],[.5-fraction/2,1],[.5+fraction/2,1]];
            const color=r.settings.diamond?'#ffbd66':'#58d6e8';
            ctx.strokeStyle=color+'50';ctx.setLineDash([6,8]);ctx.beginPath();
            for(const [i,p] of (r.settings.diamond?[0,1,2,3,0]:[0,1,3,2,0]).map((k,i)=>[i,points[k]])){if(i)ctx.lineTo(x+p[0]*rw,y+p[1]*rh);else ctx.moveTo(x+p[0]*rw,y+p[1]*rh);}ctx.stroke();ctx.setLineDash([]);
            for(let i=0;i<points.length;++i){const [u,v]=points[i],px=x+u*rw,py=y+v*rh;
                ctx.fillStyle=color+'25';ctx.beginPath();ctx.arc(px,py,24,0,Math.PI*2);ctx.fill();ctx.strokeStyle=color;ctx.lineWidth=2;ctx.beginPath();ctx.arc(px,py,13,0,Math.PI*2);ctx.stroke();ctx.fillStyle=color;ctx.beginPath();ctx.arc(px,py,5,0,Math.PI*2);ctx.fill();
                const lx=Math.max(35,Math.min(w-35,px)),ly=Math.max(38,Math.min(h-46,py));ctx.font='600 12px system-ui';ctx.textAlign='center';ctx.fillText('LED '+(i+1),lx,ly);
            }
            if(!r.settings.diamond){
                const left=x+points[0][0]*rw,right=x+points[1][0]*rw,dy=Math.max(y+55,64);
                ctx.strokeStyle=color+'80';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(left,dy-5);ctx.lineTo(left,dy+5);ctx.moveTo(left,dy);ctx.lineTo(right,dy);ctx.moveTo(right,dy-5);ctx.lineTo(right,dy+5);ctx.stroke();
                const text=cm(r.base)+' cm';ctx.font='600 13px system-ui';ctx.textAlign='center';const tw=ctx.measureText(text).width;
                ctx.fillStyle='#0b1421';ctx.fillRect((left+right-tw)/2-7,dy-10,tw+14,20);ctx.fillStyle=color;ctx.fillText(text,(left+right)/2,dy+5);
            }
        }
        destroy(){clearTimeout(this.timer);}
    }
    OF.EmitterAlignment=EmitterAlignment;
})(typeof globalThis!=='undefined'?globalThis:this);
