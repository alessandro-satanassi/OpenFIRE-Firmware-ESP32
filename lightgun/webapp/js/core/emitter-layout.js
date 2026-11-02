/* EN: Offline emitter planner. Geometry and uncertainty only; no device writes.
   IT: Assistente emettitori offline. Solo geometria e stime; nessuna scrittura sul dispositivo.
   The camera is a record, never a model-specific algorithm. Add-on optics are
   explicitly approximate: a conversion factor is not a measured IR lens profile. */
(function (root) {
    'use strict';
    const OF = root.OF = root.OF || {};
    const R = Math.PI / 180;
    const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
    const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const norm = (v) => { const n = Math.hypot(...v); return v.map((x) => x / n); };
    const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
    const DEFAULT_RATIO = 0.711;
    function optics(record, factor = 1, custom = null) {
        if (!record || !(record.horizontal > 0 && record.vertical > 0 && record.focal > 0))
            throw new Error('Missing camera optics');
        const o = { ...record, estimated: factor !== 1, refined: false };
        if (custom && factor === 1) {
            const k1=Number.isFinite(custom.k1)?custom.k1:record.k1,k2=Number.isFinite(custom.k2)?custom.k2:record.k2;
            if(Number.isFinite(custom.k1)||Number.isFinite(custom.k2)){
                const curve=OF.LensFit?.curveValid(k1,k2);
                if(curve?.ok){o.k1=k1;o.k2=k2;o.min=Math.min(o.min,curve.min);o.max=Math.max(o.max,curve.max);o.refined=true;}
            }
            if(Number.isFinite(custom.focal)&&custom.focal>=record.focal*.5&&custom.focal<=record.focal*2){o.focal=custom.focal;o.refined=true;}
            const nativeGain = 1 + record.k1 + record.k2, gain = 1 + o.k1 + o.k2;
            const scale = gain / nativeGain * record.focal / o.focal;
            if (!(scale > 0)) throw new Error('Invalid refined optics');
            o.horizontal = 2 * Math.atan(Math.tan(record.horizontal * R / 2) * scale) / R;
            o.vertical = 2 * Math.atan(Math.tan(record.vertical * R / 2) * scale) / R;
        }
        if (factor !== 1) {
            if (!(factor >= .3 && factor <= 1)) throw new Error('Invalid converter');
            o.horizontal = 2 * Math.atan(Math.tan(record.horizontal * R / 2) / factor) / R;
            o.vertical = 2 * Math.atan(Math.tan(record.vertical * R / 2) / factor) / R;
            o.k1 = o.k2 = 0; // Unknown added distortion: never pretend to have measured it.
        }
        o.tx = Math.tan(o.horizontal * (1 - o.tolerance) * R / 2);
        o.ty = Math.tan(o.vertical * (1 - o.tolerance) * R / 2);
        o.sideGain = 1 + o.k1 + o.k2;
        return o;
    }
    function dimensions(inches, aspect) {
        return { height: inches * .0254 / Math.sqrt(1 + aspect * aspect),
            width: inches * .0254 * aspect / Math.sqrt(1 + aspect * aspect) };
    }
    function settings(input) {
        const finite = (v, a, b, fallback) => Number.isFinite(+v) ? clamp(+v, a, b) : fallback;
        return { inches: finite(input.inches, 20, 120, 55), aspect: finite(input.aspect, 1, 4, 16/9),
            distance: finite(input.distance, .3, 20, 2), diamond: !!input.diamond };
    }
    function poses(s, distance) {
        const { width: w, height: h } = dimensions(s.inches, s.aspect);
        // EN: One play profile: side moves up to 20% of the screen width, the gun held straight
        //     (all four LEDs) or rolled up to 20° and 45°: then one LED may leave the view, the
        //     firmware rebuilds it from the other three (spare). Larger rolls are rare: not checked.
        // IT: Un solo profilo di gioco: spostamenti laterali fino al 20% della larghezza schermo,
        //     arma dritta (tutti e 4 i LED) o ruotata fino a 20° e 45°: allora un LED puo' uscire
        //     dall'inquadratura, il firmware lo ricostruisce dagli altri tre (spare). Oltre: raro.
        const sides = [-w*.2, 0, w*.2];
        const rolls = [-45, -20, 0, 20, 45];
        const list = [];
        for (const side of sides) for (const up of [-h*.1, h*.1]) for (const ax of [-.5, 0, .5]) for (const ay of [-.5, 0, .5]) {
            const f = norm([ax*w-side, ay*h-up, distance]), right = norm([f[2], 0, -f[0]]), down = cross(f, right);
            for (const roll of rolls) {
                const c = Math.cos(roll*R), sn = Math.sin(roll*R);
                list.push({ eye: [side, up, -distance], f, spare: roll ? 1 : 0,
                    right: right.map((v, i) => v*c+down[i]*sn), down: down.map((v, i) => v*c-right[i]*sn) });
            }
        }
        return list;
    }
    function leds(s, ratio) {
        const { width: w, height: h } = dimensions(s.inches, s.aspect);
        return s.diamond ? [[0,-h/2,0],[w/2,0,0],[0,h/2,0],[-w/2,0,0]]
            : [[-ratio*h/2,-h/2,0],[ratio*h/2,-h/2,0],[-ratio*h/2,h/2,0],[ratio*h/2,h/2,0]];
    }
    function project(o, pose, led) {
        const v = led.map((x, i) => x-pose.eye[i]), z = dot(v, pose.f);
        if (!(z > 0)) return [1e6, 1e6];
        let x = dot(v, pose.right) / (z * o.tx) * o.sideGain;
        let y = dot(v, pose.down) / (z * o.ty) * o.sideGain;
        const r = Math.hypot(x, y);
        if (r && (o.k1 || o.k2)) {
            let lo = 0, hi = Math.max(2, r / o.min + 1);
            for (let n=0; n<18; ++n) {
                const m=(lo+hi)/2, q=m*m, gain=clamp(1+o.k1*q+o.k2*q*q, o.min, o.max);
                if (m*gain < r) lo=m; else hi=m;
            }
            const scale=(lo+hi)/(2*r); x*=scale; y*=scale;
        }
        return [x,y];
    }
    function coverage(s, o, ratio, scenePoses) {
        const points=leds(s,ratio); let full=0, visible=0, margin=1;
        for (const pose of scenePoses) {
            const edge=points.map(led=>{const q=project(o,pose,led);return 1-Math.max(Math.abs(q[0]),Math.abs(q[1]));}).sort((a,b)=>a-b);
            const need=edge[pose.spare||0];
            margin=Math.min(margin,need); if(need>=0) ++full; visible+=edge.filter(m=>m>=0).length;
        }
        return { full: full/scenePoses.length, visible: visible/(4*scenePoses.length), margin };
    }
    // EN: Test the radial boundary directly for the distance search (no inversion).
    // IT: Verifica diretta del bordo radiale nella ricerca distanza (nessuna inversione).
    function fits(s,o,ratio,scenePoses,border=.05) {
        const points=leds(s,ratio),limit=1-border;
        for(const pose of scenePoses){
            let out=0;
            for(const led of points){
                const v=led.map((x,i)=>x-pose.eye[i]),z=dot(v,pose.f);
                let bad=!(z>0);
                if(!bad){
                    const x=dot(v,pose.right)/(z*o.tx)*o.sideGain,y=dot(v,pose.down)/(z*o.ty)*o.sideGain,m=Math.max(Math.abs(x),Math.abs(y));
                    if(m){const r=(x*x+y*y)*limit*limit/(m*m),gain=clamp(1+o.k1*r+o.k2*r*r,o.min,o.max);bad=m>limit*gain;}
                }
                if(bad&&++out>(pose.spare||0))return false;
            }
        }
        return true;
    }
    // Conditioning proxy: evaluate the inverse homography at the camera centre.
    // Sensor jitter is perturbed BEFORE lens correction, as in the real pipeline.
    function query(o, raw, aspect, ratio) {
        const q=raw.map(([x,y])=>{const r=x*x+y*y,g=clamp(1+o.k1*r+o.k2*r*r,o.min,o.max);return[x*g,y*g];});
        const [[x0,y0],[x1,y1],[x3,y3],[x2,y2]]=q;
        const dx1=x1-x2,dy1=y1-y2,dx2=x3-x2,dy2=y3-y2,sx=x0-x1+x2-x3,sy=y0-y1+y2-y3;
        const det=dx1*dy2-dx2*dy1; if(Math.abs(det)<1e-15) return null;
        const g=(sx*dy2-dx2*sy)/det,h=(dx1*sy-sx*dy1)/det;
        const a=x1-x0+g*x1,b=x3-x0+h*x3,c=x0,d=y1-y0+g*y1,e=y3-y0+h*y3,f=y0;
        const den=a*e-b*d; if(Math.abs(den)<1e-15) return null;
        return [((b*f-c*e)/den-.5)*ratio/aspect, (c*d-a*f)/den-.5];
    }
    function conditioning(s,o,ratio,scenePoses) {
        if(s.diamond) return 1;
        const points=leds(s,ratio); let sum=0,n=0;
        // A representative subset suffices to rank geometry; coverage uses every pose.
        for(let i=0;i<scenePoses.length;i+=Math.max(1,Math.floor(scenePoses.length/36))) {
            const raw=points.map(p=>project(o,scenePoses[i],p));
            if(raw.some(q=>Math.max(Math.abs(q[0]),Math.abs(q[1]))>1)) continue;
            let cost=0;
            for(let k=0;k<4;++k) for(let axis=0;axis<2;++axis) {
                const value=raw[k][axis],eps=1e-4;raw[k][axis]=value+eps;const a=query(o,raw,s.aspect,ratio);
                raw[k][axis]=value-eps;const b=query(o,raw,s.aspect,ratio);raw[k][axis]=value;
                if(!a||!b) return Infinity;
                const pixel=2/(axis?o.sensorY:o.sensorX);
                // Equal physical aiming error on X/Y (both measured in screen-height units).
                cost+=(((a[0]-b[0])*s.aspect)**2+(a[1]-b[1])**2)*pixel*pixel/(4*eps*eps);
            }
            sum+=cost; ++n;
        }
        return n ? Math.sqrt(sum/n) : Infinity;
    }
    // EN: The layout must keep working when the player steps forward: it is checked at this
    //     share of the playing distance. IT: La disposizione deve reggere anche un passo avanti:
    //     si verifica a questa frazione della distanza di gioco.
    const STEP_FORWARD = .85;
    // EN: Shapes close to a square let the firmware swap the LEDs when seen up close (tracking
    //     simulations): upright up to 0.85, wide from 1.5 or at the screen corners (4:3 and wider).
    //     1.3 and 1.4 only from 4 screen heights, even one step forward: no swaps seen there.
    //     Never narrower than 0.711, the proportion recommended so far.
    // IT: Rettangoli quasi quadrati fanno scambiare i LED al firmware visti da vicino (simulazioni
    //     del tracking): verticali fino a 0,85, larghi da 1,5 o agli angoli (dal 4:3 in su).
    //     1,3 e 1,4 solo da 4 altezze schermo, anche col passo avanti: li' nessuno scambio.
    //     Mai piu' stretto di 0,711, la proporzione consigliata finora.
    const FAR_HEIGHTS = 4;
    function candidates(s,distance=0) {
        if(s.diamond) return [s.aspect];
        const far=distance*STEP_FORWARD>=FAR_HEIGHTS*dimensions(s.inches,s.aspect).height;
        return [DEFAULT_RATIO,.75,.8,.85,1.3,1.4,1.5,1.6,1.68,s.aspect]
            .filter((r,i,a)=>r<=s.aspect && (r<=.85 || r>=1.5 || (r===s.aspect && r>=1.3) || (far && (r===1.3 || r===1.4))) && a.indexOf(r)===i);
    }
    function rowsAt(s,o,distance) {
        const check=poses(s,distance*STEP_FORWARD);
        return candidates(s,distance).map(r=>({ratio:r,...coverage(s,o,r,check)}));
    }
    // EN/IT: 'good', 'limited' or null (no layout keeps the LEDs in view) / nessuna disposizione va.
    // Same thresholds as choose(), with the faster direct boundary test.
    function viable(s,o,distance) {
        const check=poses(s,distance*STEP_FORWARD),list=candidates(s,distance);
        return list.some(r=>fits(s,o,r,check))?'good':list.some(r=>fits(s,o,r,check,0))?'limited':null;
    }
    // EN/IT: Best layout at this playing distance, or null when no layout keeps the LEDs in view.
    function choose(s,o,distance) {
        const rows=rowsAt(s,o,distance);
        const safe=rows.filter(r=>r.margin>=.05),complete=rows.filter(r=>r.full===1);
        const pool=safe.length?safe:complete;
        if(!pool.length) return null;
        const scenePoses=poses(s,distance);
        for(const r of pool) r.noise=conditioning(s,o,r.ratio,scenePoses);
        pool.sort((a,b)=>(Number.isFinite(a.noise)&&Number.isFinite(b.noise)?a.noise-b.noise:0)||b.margin-a.margin);
        return {...pool[0],status:safe.length?'good':'limited'};
    }
    function calculate(input, record, factor=1, custom=null, withDistance=true) {
        const s=settings(input),o=optics(record,factor,custom),dim=dimensions(s.inches,s.aspect);
        let chosen=choose(s,o,s.distance),moveTo=null;
        if(!chosen) {
            // EN: Too close: suggest the layout of the nearest distance that works, and that distance,
            //     so the advice stays the same while the player moves back.
            // IT: Troppo vicino: si suggerisce la disposizione della distanza piu' vicina che funziona,
            //     e quella distanza, cosi' il consiglio resta lo stesso mentre il giocatore arretra.
            const want=viable(s,o,20);
            if(want) {
                let lo=.3,hi=20; // fixed interval: the same answer from every closer distance
                for(let n=0;n<11;++n){const mid=(lo+hi)/2;if(viable(s,o,mid)===want)hi=mid;else lo=mid;}
                moveTo=Math.ceil(hi*20)/20;chosen=choose(s,o,moveTo);
            }
            if(!chosen) {
                let pool=rowsAt(s,o,s.distance);
                const best=Math.max(...pool.map(r=>r.full));pool=pool.filter(r=>r.full===best);
                const v=Math.max(...pool.map(r=>r.visible));chosen=pool.find(r=>r.visible===v);
            }
            chosen={...chosen,status:'close'};
        }
        const result={...chosen,...dim,settings:s,optics:o,base:s.diamond?dim.width:chosen.ratio*dim.height,
            minDistance:null,moveTo};
        if(withDistance) {
            // EN/IT: Same step-forward margin as the layout choice / Stesso margine della scelta.
            const accepts=d=>fits(s,o,chosen.ratio,poses(s,d*STEP_FORWARD))&&candidates(s,d).includes(chosen.ratio);
            let lo=.3,hi=20;
            if(accepts(hi)) {for(let n=0;n<12;++n){const mid=(lo+hi)/2;if(accepts(mid))hi=mid;else lo=mid;}result.minDistance=Math.ceil(hi*20)/20;}
        }
        return result;
    }
    OF.EmitterLayout={calculate,candidates,dimensions,optics,settings,poses,leds,project,coverage,fits,DEFAULT_RATIO,STEP_FORWARD,FAR_HEIGHTS};
})(typeof globalThis!=='undefined'?globalThis:this);
