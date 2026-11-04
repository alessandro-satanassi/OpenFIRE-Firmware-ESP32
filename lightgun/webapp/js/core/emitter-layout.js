/* EN: Offline emitter planner. Geometry and uncertainty only; no device writes.
   IT: Assistente emettitori offline. Solo geometria e stime; nessuna scrittura sul dispositivo.
   The camera is a record, never a model-specific algorithm. Add-on optics use
   the entered conversion factor combined with the native field of view. */
(function (root) {
    'use strict';
    const OF = root.OF = root.OF || {};
    const R = Math.PI / 180;
    const PLAY_ROLLS = [-15,-10,-5,0,5,10,15].map(roll=>({
        c:Math.cos(roll*R),sn:Math.sin(roll*R),spare:roll?1:0
    }));
    const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
    const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const norm = (v) => { const n = Math.hypot(...v); return v.map((x) => x / n); };
    const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
    const DEFAULT_RATIO = 0.711;
    const MAX_DISTANCE = 10, DISTANCE_STEP = .1;
    // EN: The layout must keep working when the player leans forward: it is checked at this
    //     share of the playing distance (gun to screen). IT: La disposizione deve reggere anche
    //     sporgendosi in avanti: si verifica a questa frazione della distanza pistola-schermo.
    const STEP_FORWARD = .9;
    // EN: Near-square wide shapes keep the stricter 15% reserve: stability comes first.
    // IT: I rettangoli larghi quasi quadrati tengono la riserva piu' severa del 15%: prima la stabilita'.
    const FAR_STEP = .85;
    // EN: Shapes close to a square let the firmware swap the LEDs when seen up close (tracking
    //     simulations): upright up to 0.85, wide from 1.5 or at the screen corners (4:3 and wider).
    //     1.3 and 1.4 keep the extra distance reserve of 4 screen heights.
    //     0.711 is the minimum: narrower shapes trade precision for a few centimetres.
    // IT: Rettangoli quasi quadrati fanno scambiare i LED al firmware visti da vicino (simulazioni
    //     del tracking): verticali fino a 0,85, larghi da 1,5 o agli angoli (dal 4:3 in su).
    //     1,3 e 1,4 mantengono la riserva aggiuntiva di 4 altezze schermo.
    //     0,711 e' il minimo: rettangoli piu' stretti cedono precisione per pochi centimetri.
    const FAR_HEIGHTS = 4;
    const distanceCeil = value => Math.ceil(value * (1 / DISTANCE_STEP)) / (1 / DISTANCE_STEP);
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
            o.focal *= factor; o.k1 = o.k2 = 0; // Ideal converter model; the factor describes the combined FOV.
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
            distance: finite(input.distance, .3, MAX_DISTANCE, 2), diamond: !!input.diamond };
    }
    function poses(s, distance) {
        const { width: w, height: h } = dimensions(s.inches, s.aspect);
        // EN: One play profile: side moves up to 20% of the screen width, the gun held straight
        //     (all four LEDs) or rolled up to 15°: then one LED may leave the view, the
        //     firmware rebuilds it from the other three (spare). Larger rolls are rare: not checked.
        // IT: Un solo profilo di gioco: spostamenti laterali fino al 20% della larghezza schermo,
        //     arma dritta (tutti e 4 i LED) o ruotata fino a 15°: allora un LED puo' uscire
        //     dall'inquadratura, il firmware lo ricostruisce dagli altri tre (spare). Oltre: raro.
        const sides = [-w*.2, 0, w*.2];
        const list = [];
        for (const side of sides) for (const up of [-h*.1, h*.1]) for (const ax of [-.5, 0, .5]) for (const ay of [-.5, 0, .5]) {
            const f = norm([ax*w-side, ay*h-up, distance]), right = norm([f[2], 0, -f[0]]), down = cross(f, right);
            for (const {c,sn,spare} of PLAY_ROLLS) {
                list.push({ eye: [side, up, -distance], f, spare,
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
    // EN: Pointing sensitivity in screen-height units, using the firmware radial map.
    // IT: Sensibilita della mira in unita di altezza schermo, con la stessa mappa radiale.
    function pointer(s,o,ratio,raw) {
        const q=raw.map(([x,y])=>{const r=x*x+y*y,g=clamp(1+o.k1*r+o.k2*r*r,o.min,o.max);return[x*g,y*g];});
        const [[x0,y0],[x1,y1],[x3,y3],[x2,y2]]=q;
        const dx1=x1-x2,dy1=y1-y2,dx2=x3-x2,dy2=y3-y2,sx=x0-x1+x2-x3,sy=y0-y1+y2-y3;
        const det=dx1*dy2-dx2*dy1;if(Math.abs(det)<1e-12)return null;
        const g=(sx*dy2-dx2*sy)/det,h=(dx1*sy-sx*dy1)/det;
        const a=x1-x0+g*x1,b=x3-x0+h*x3,c=x0,d=y1-y0+g*y1,e=y3-y0+h*y3,f=y0;
        const den=a*e-b*d;if(Math.abs(den)<1e-12)return null;
        return [((b*f-c*e)/den-.5)*ratio,(c*d-a*f)/den-.5];
    }
    // EN: Approximate hidden points with the same virtual-camera geometry as Anchor.
    // IT: Stima dei punti nascosti con la stessa geometria di camera virtuale di Anchor.
    function anchorPlane(points,focal) {
        if(!(focal>0))return {to:p=>p,from:p=>p};
        const mean=[0,0,0];
        for(const [x,y] of points){const n=Math.hypot(x,y,focal);mean[0]+=x/n;mean[1]+=y/n;mean[2]+=focal/n;}
        const [mx,my,mz]=norm(mean),vx=my,vy=-mx,c=1/(1+mz);
        const r=[1-vy*vy*c,vx*vy*c,vy,vx*vy*c,1-vx*vx*c,-vx,-vy,vx,1-(vx*vx+vy*vy)*c];
        return {
            to:([x,y])=>{const z=r[6]*x+r[7]*y+r[8]*focal;return z>0?[(r[0]*x+r[1]*y+r[2]*focal)*focal/z,(r[3]*x+r[4]*y+r[5]*focal)*focal/z]:null;},
            from:([x,y])=>{const z=r[2]*x+r[5]*y+r[8]*focal;return z>0?[(r[0]*x+r[3]*y+r[6]*focal)*focal/z,(r[1]*x+r[4]*y+r[7]*focal)*focal/z]:null;}
        };
    }
    // EN: Keep the complete reference fixed while perturbing current measurements.
    // IT: Il riferimento completo resta fisso mentre si perturbano le misure correnti.
    function reconstruction(o,reference,visible,focal) {
        if(visible.length===4)return current=>current;
        if(!visible.length)return ()=>reference;
        const sensor=p=>[p[0]*o.sensorX/2,p[1]*o.sensorY/2];
        const old=visible.map(i=>sensor(reference[i]));
        if(visible.length===1)return current=>{
            const now=sensor(current[visible[0]]),dx=now[0]-old[0][0],dy=now[1]-old[0][1];
            return reference.map(p=>[p[0]+2*dx/o.sensorX,p[1]+2*dy/o.sensorY]);
        };
        const P=anchorPlane(old,focal),p=old.map(P.to);
        if(p.some(x=>!x))return null;
        const a=[p[1][0]-p[0][0],p[1][1]-p[0][1]],b=visible.length===2?[-a[1],a[0]]:[p[2][0]-p[0][0],p[2][1]-p[0][1]];
        const det=a[0]*b[1]-a[1]*b[0];if(Math.abs(det)<1e-12)return null;
        const hidden=[];
        for(let i=0;i<4;++i)if(!visible.includes(i)){
            const h=P.to(sensor(reference[i]));if(!h)return null;
            const x=h[0]-p[0][0],y=h[1]-p[0][1];
            hidden.push({i,s:(x*b[1]-y*b[0])/det,t:(a[0]*y-a[1]*x)/det});
        }
        return current=>{
            const now=visible.map(i=>sensor(current[i])),Q=anchorPlane(now,focal),q=now.map(Q.to);
            if(q.some(x=>!x))return null;
            const u=[q[1][0]-q[0][0],q[1][1]-q[0][1]],v=visible.length===2?[-u[1],u[0]]:[q[2][0]-q[0][0],q[2][1]-q[0][1]];
            const result=current.slice();
            for(const {i,s,t} of hidden){
                const h=Q.from([q[0][0]+s*u[0]+t*v[0],q[0][1]+s*u[1]+t*v[1]]);if(!h)return null;
                result[i]=[2*h[0]/o.sensorX,2*h[1]/o.sensorY];
            }
            return result;
        };
    }
    // EN: A geometric index, not a measured pointing error. It combines reconstruction
    // bias and sensitivity at camera-output and sensor-pixel scales, in mean and worst poses.
    // No camera IDs, fitted per-model weights or calculations in the firmware loop.
    // IT: Indice geometrico, non errore di mira misurato. Combina errore di ricostruzione
    // e sensibilita alle scale di uscita camera e pixel fisico, in media e nelle pose peggiori.
    // Nessun ID camera, peso tarato per modello o calcolo nel ciclo del firmware.
    function precision(s,o,ratio,distance) {
        const {width,height}=dimensions(s.inches,s.aspect),points=leds(s,ratio),focal=o.sensorX*.5*o.focal;
        const resolution=[o.coordinateX||o.sensorX,o.coordinateY||o.sensorY];
        const invalid={fine:Infinity,coarse:Infinity,worstFine:Infinity,worstCoarse:Infinity,full:0};
        let fine=0,coarse=0,worstFine=0,worstCoarse=0,count=0,full=0;
        // EN: Radial optics and centred LEDs are mirror symmetric: the same average
        // needs one height and two lateral positions (the mirrored side counts twice).
        // IT: Ottica radiale e LED centrati sono simmetrici: stessa media con una sola
        // altezza e due posizioni laterali, contando due volte il lato speculare.
        for(const depth of[distance,distance*STEP_FORWARD])for(const side of[0,.2])for(const up of[.1]){
            const weight=side?2:1;
            const eye=[side*width,up*height,-depth],f=norm([-eye[0],-eye[1],depth]),right=norm([f[2],0,-f[0]]),down=cross(f,right);
            const centre=points.map(p=>project(o,{eye,f,right,down},p));
            for(const [ax,ay] of[[0,0],[-.5,-.5],[.5,-.5],[-.5,.5],[.5,.5]])for(const roll of[-15,0,15]){
                const view=t=>{
                    const f=norm([t*ax*width-eye[0],t*ay*height-eye[1],depth]),base=norm([f[2],0,-f[0]]),vertical=cross(f,base),c=Math.cos(t*roll*R),sn=Math.sin(t*roll*R);
                    const pose={eye,f,right:base.map((v,i)=>v*c+vertical[i]*sn),down:vertical.map((v,i)=>v*c-base[i]*sn)};
                    return points.map(p=>project(o,pose,p));
                };
                const raw=view(1),visible=[0,1,2,3].filter(i=>Math.max(...raw[i].map(Math.abs))<=1);
                let ref=centre;
                // EN/IT: A nearby complete view along a centre-to-target sweep /
                // Riferimento completo vicino lungo una mira dal centro al bersaglio.
                if(visible.length<4&&centre.every(p=>Math.max(...p.map(Math.abs))<=1)){
                    let lo=0,hi=1;
                    for(let n=0;n<Math.ceil(Math.log2(Math.max(...resolution)));++n){const mid=(lo+hi)/2,q=view(mid);if(q.every(p=>Math.max(...p.map(Math.abs))<=1)){lo=mid;ref=q;}else hi=mid;}
                }
                const reconstruct=reconstruction(o,ref,visible,focal);
                const evaluate=raw=>{const quad=reconstruct&&reconstruct(raw);return quad&&pointer(s,o,ratio,quad);};
                const q=evaluate(raw);if(!q)return invalid;
                let viewFine=(q[0]-ax*s.aspect)**2+(q[1]-ay)**2,viewCoarse=viewFine;
                if(visible.length===4)full+=weight;
                for(const i of visible)for(let axis=0;axis<2;++axis){
                    const value=raw[i][axis],epsilon=1/resolution[axis];
                    raw[i][axis]=value+epsilon;const a=evaluate(raw);
                    raw[i][axis]=value-epsilon;const b=evaluate(raw);raw[i][axis]=value;
                    if(!a||!b)return invalid;
                    // EN/IT: Uniform quantisation variance: step squared / 12 /
                    // Varianza di quantizzazione uniforme: passo al quadrato / 12.
                    const delta=((a[0]-b[0])**2+(a[1]-b[1])**2)/12;
                    viewFine+=delta;viewCoarse+=delta*(resolution[axis]/(axis?o.sensorY:o.sensorX))**2;
                }
                fine+=weight*viewFine;coarse+=weight*viewCoarse;
                worstFine=Math.max(worstFine,viewFine);worstCoarse=Math.max(worstCoarse,viewCoarse);
                count+=weight;
            }
        }
        return {fine:fine/count,coarse:coarse/count,worstFine,worstCoarse,full:full/count};
    }

    function candidates(s,distance=0) {
        if(s.diamond) return [s.aspect];
        const far=distance*FAR_STEP>=FAR_HEIGHTS*dimensions(s.inches,s.aspect).height;
        const list=[DEFAULT_RATIO,.8,.85,1.3,1.4,1.5,1.6,1.68,s.aspect];
        return list
            .filter((r,i,a)=>r<=s.aspect && (r<=.85 || r>=1.5 || (r===s.aspect && r>=1.3) || (far && (r===1.3 || r===1.4))) && a.indexOf(r)===i);
    }
    // EN: Cache the whole progression once per setup. Start with the first usable
    // placement, even if the chosen playing distance is too short. Precision may
    // override the visibility-based reference only if both spatial scales improve.
    // IT: La progressione si calcola una volta per configurazione. Si parte dalla
    // prima disposizione utile, anche se la distanza scelta e troppo breve. La stima
    // di precisione cambia il riferimento solo se migliora a entrambe le scale.
    let planKey='',cachedPlan;
    function distancePlan(s,o) {
        const key=JSON.stringify([s.inches,s.aspect,s.diamond,o]);
        if(key===planKey)return cachedPlan;
        const height=dimensions(s.inches,s.aspect).height;
        const immediate=candidates(s,0);
        // EN/IT: Reuse the identical search endpoints for every rectangle /
        // Riutilizza gli stessi estremi della ricerca per tutti i rettangoli.
        const endpoints=new Map([[.3,poses(s,.3)],[MAX_DISTANCE,poses(s,MAX_DISTANCE)]]);
        const rows=candidates(s,MAX_DISTANCE).map(ratio=>{
            // EN/IT: Physical limits first; apply the FAR rule separately / Prima i limiti
            // ottici; la regola FAR si applica a parte, senza dividerla per STEP_FORWARD.
            const limit=border=>{
                const accepts=d=>fits(s,o,ratio,endpoints.get(d)||poses(s,d),border);
                if(!accepts(MAX_DISTANCE))return Infinity;
                let lo=.3,hi=MAX_DISTANCE;
                if(accepts(lo))return lo;
                for(let i=0;i<16;++i){const mid=(lo+hi)/2;if(accepts(mid))hi=mid;else lo=mid;}
                return hi;
            };
            const gate=immediate.includes(ratio)?0:FAR_HEIGHTS*height/FAR_STEP;
            const usable=Math.max(limit(0),gate),good=Math.max(limit(.05)/STEP_FORWARD,gate);
            return {ratio,usable,good,at:ratio<=DEFAULT_RATIO?usable:good};
        });
        const reference=rows.reduce((a,b)=>b.usable<a.usable||(b.usable===a.usable&&b.ratio>a.ratio)?b:a);
        let firstGood=reference.good;
        for(const row of rows)if(row.ratio>reference.ratio){
            if(row.at>firstGood)row.at=row.good;
            firstGood=Math.min(firstGood,row.good);
        }
        const scoreCache=new Map();
        const score=(row,distance)=>{
            const key=row.ratio+':'+distance;
            if(!scoreCache.has(key)){
                const q=precision(s,o,row.ratio,distance);
                scoreCache.set(key,q);
            }
            return scoreCache.get(key);
        };
        // EN: Preserve the reference when the model is ambiguous. A preferred layout
        // must improve both mean and worst poses and keep the first green control step;
        // visibility breaks precision ties.
        // IT: Se il modello e ambiguo si conserva il riferimento. La disposizione preferita
        // deve migliorare media e pose peggiori e conservare il primo passo verde;
        // a pari precisione decide la visibilita.
        const choose=(baseline,eligible,distance)=>{
            const alternatives=eligible.filter(r=>r!==baseline&&distanceCeil(r.good)<=distanceCeil(baseline.good));
            if(!alternatives.length)return baseline;
            const ref=score(baseline,distance);let best=baseline,bestRegret=1;
            const bounds=['fine','coarse','worstFine','worstCoarse'];
            for(const row of alternatives){
                const q=score(row,distance);
                if(!Number.isFinite(q.coarse)||!bounds.every(k=>q[k]<=ref[k]))continue;
                const regret=Math.max(...bounds.map(k=>q[k]/(Number.isFinite(ref[k])&&ref[k]>0?ref[k]:1)));
                if(regret<bestRegret||(regret===bestRegret&&q.full>score(best,distance).full)){best=row;bestRegret=regret;}
            }
            return best;
        };
        const firstDistance=Math.min(MAX_DISTANCE,distanceCeil(reference.usable));
        let baseline=reference;
        for(const row of rows)if(row.ratio>baseline.ratio&&row.at<=firstDistance)baseline=row;
        const initial=rows.filter(r=>r.usable<=firstDistance);
        const first=choose(baseline,initial.length?initial:rows,firstDistance),steps=[];
        let selected=first;
        // EN/IT: Only wider successors; no dependence on the slider history /
        // Solo successori piu larghi, senza dipendere dai movimenti dello slider.
        for(let n=Math.round(firstDistance/DISTANCE_STEP)+1;n<=MAX_DISTANCE/DISTANCE_STEP;++n){
            const distance=n/(1/DISTANCE_STEP);
            let baseline=selected;
            for(const row of rows)if(row.ratio>baseline.ratio&&row.at<=distance)baseline=row;
            const eligible=rows.filter(r=>r.ratio>=selected.ratio&&(r===selected||r.at<=distance));
            const next=choose(baseline,eligible,distance);
            if(next!==selected){steps.push({at:distance,row:next});selected=next;}
        }
        // EN/IT: First green distance of this same progression (the ideal-distance button
        // finds the same value) / Prima distanza verde della stessa progressione.
        let green=null;
        for(let n=Math.round(.3/DISTANCE_STEP);green===null&&n<=MAX_DISTANCE/DISTANCE_STEP;++n){
            const distance=n/(1/DISTANCE_STEP);let row=first;
            for(const step of steps)if(step.at<=distance)row=step.row;
            if(distance>=row.good)green=distance;
        }
        planKey=key;
        return cachedPlan={rows,first,steps,green};

    }
    function calculate(input, record, factor=1, custom=null, withDistance=true) {
        const s=settings(input),o=optics(record,factor,custom),dim=dimensions(s.inches,s.aspect);
        const plan=distancePlan(s,o);
        // EN: A red result keeps the next useful placement. The colour and distance
        // labels still refer to the actual input and to this same displayed rectangle.
        // IT: In rosso resta la prima disposizione utile successiva. Colore e distanze
        // si riferiscono sempre al valore impostato e al rettangolo effettivamente mostrato.
        let chosen=plan.first;
        for(const step of plan.steps)if(step.at<=s.distance)chosen=step.row;
        const status=s.distance<chosen.usable?'close':s.distance<chosen.good?'limited':'good';
        return {ratio:chosen.ratio,status,...dim,settings:s,optics:o,
            base:s.diamond?dim.width:chosen.ratio*dim.height,
            // EN: Green: the displayed layout's own limit. Otherwise the first green distance:
            //     farther back the progression may already show wider LEDs.
            // IT: In verde: il limite della disposizione mostrata. Altrimenti la prima distanza
            //     verde: piu' indietro la progressione puo' gia' mostrare LED piu' distanziati.
            minDistance:!withDistance?null:status==='good'?distanceCeil(chosen.good):plan.green,
            moveTo:status==='close'&&chosen.usable<=MAX_DISTANCE?distanceCeil(chosen.usable):null};
    }
    OF.EmitterLayout={precision,calculate,candidates,dimensions,optics,settings,poses,leds,project,coverage,fits,DEFAULT_RATIO,MAX_DISTANCE,DISTANCE_STEP,STEP_FORWARD,FAR_STEP,FAR_HEIGHTS};
})(typeof globalThis!=='undefined'?globalThis:this);
