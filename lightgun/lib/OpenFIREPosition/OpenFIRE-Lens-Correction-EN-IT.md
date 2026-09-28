# OpenFIRE — Lens Correction / Correzione ottica

[English](#english) · [Italiano](#italiano)

---

<a id="english"></a>

## English

### Purpose and scope

This document explains how to estimate initial values for `lensRadialK1` and `lensRadialK2` from nominal optical specifications, and how to choose `lensCorrectionMin` and `lensCorrectionMax` without clipping the expected correction within the sensor field.

The derivation specifically covers the **PAJ7025R2 and PAJ7025R3**, which have square sensors and equal horizontal and vertical fields of view.

**These are theoretical starting values, not a measured lens calibration.** Focal length and field-of-view specifications at the boundaries do not describe the complete distortion curve. Two coefficients can fit two nominal points exactly, but this does not establish accuracy between those points or on an individual camera. Real-camera validation remains necessary.

Numerical inputs come from the supplied PixArt specifications. The formulas below are a derivation for OpenFIRE's correction model, not official PixArt calibration coefficients.

### 1. Required specifications

| Datasheet item | Symbol | PAJ7025R2 | PAJ7025R3 | Purpose |
|---|---|---:|---:|---|
| Effective Focal Length | `f` | 1.484307 mm | 0.378 mm | Reference rectilinear image radius |
| Sensor Pixel Resolution | `N` | 98 × 98 | 98 × 98 | Physical active-area dimensions |
| Sensor Pixel Size | `p` | 11 × 11 µm | 11 × 11 µm | Physical pixel-to-mm conversion |
| Angle Field of View — Horizontal/Vertical | `FOV_s` | 38.3° | 111.3° | Constraint at a side midpoint |
| Angle Field of View — Diagonal | `FOV_d` | 52.2° | 140° | Constraint at a corner |
| Boundary notes, Y | `a`, `d` | 0.539 / 0.762 mm | 0.539 / 0.762 mm | Cross-check of physical half-dimensions |
| Distortion | `D` | < 2.8% | −30% | Indicative cross-check, not a complete curve |

The FOV values are **full angles**. Use half of each angle in the formulas, and convert degrees to radians when required by the trigonometric function.

The virtual 4096 × 4096 coordinate resolution is not the physical pixel array. Do not multiply 4096 by 11 µm: use the **98 physical pixels**.

F-number, gain, exposure, frame rate, SPI interface and power consumption do not enter this geometric estimate. Some affect detection quality, but they do not determine K1/K2. Image Circle, Back Focal Length and lens-element count do not replace the effective focal length or active sensor dimensions.

### 2. Assumptions

The R2/R3 derivation assumes:

- a square sensor, square pixels, and equal horizontal/vertical FOV;
- an optical centre coincident with the centre used by the firmware;
- reported coordinates proportional to position in the active sensor area, without an already-applied DSP rectification;
- approximately radial distortion, without tangential or decentring terms;
- nominal effective focal length used as the paraxial focal length of the reference rectilinear projection;
- nominal FOV and dimension specifications referring to the same active area.

OpenFIRE normalises each axis by its own half-size. The simple derivation below is appropriate for these square PAJ sensors; it must not be transferred automatically to rectangular sensors, cropped images or a different coordinate normalisation.

### 3. Correction model used by the firmware

Let x and y be the observed coordinates, and cx/cy the centre:

```text
cx = mouseResX / 2
cy = mouseResY / 2

nx = (x - cx) / cx
ny = (y - cy) / cy
q  = nx*nx + ny*ny                  // q = r²

C(q) = 1 + K1*q + K2*q*q           // correction factor before clamping
C_limited = clamp(C(q), Min, Max)

x_corrected = cx + (x - cx)*C_limited
y_corrected = cy + (y - cy)*C_limited
```

The fields correspond to:

| Firmware field | Meaning |
|---|---|
| `lensRadialK1` | Coefficient of q, or r² |
| `lensRadialK2` | Coefficient of q², or r⁴ |
| `lensCorrectionMax` | Upper limit of the multiplicative correction factor |
| `lensCorrectionMin` | Lower limit of the multiplicative correction factor |

For the square PAJ sensors, at the idealised field boundaries:

| Position | r | q |
|---|---:|---:|
| Centre | 0 | 0 |
| Side midpoint | 1 | 1 |
| Corner | √2 | 2 |

The correction acts on **both axes**, about the optical centre. LEDs and the aim point must use the same corrected coordinate space; working on local copies prevents repeated correction from accumulating in stored calibration values.

This is an **observed → corrected** map normalised by sensor half-size. It is not automatically compatible with coefficients from a calibration package. For example, the standard OpenCV model applies distortion to ideal coordinates normalised by focal length: direction and normalisation must be converted before reusing its coefficients. [OpenCV documentation](https://docs.opencv.org/4.13.0/d9/d0c/group__calib3d.html).

Simply reversing the signs of coefficients from another model is not sufficient, especially with strong distortion.

### 4. Physical sensor dimensions

Using millimetres:

```text
p_mm = 11 / 1000             = 0.011 mm
a    = N*p_mm/2 = 98*0.011/2 = 0.539 mm
d    = a*sqrt(2)             = 0.7622611101 mm
```

Here, a is the half-side and d is the half-diagonal. The datasheet's 0.762 mm value is consistent with the rounded half-diagonal. The proposed coefficients use d = a√2, avoiding that additional rounding.

### 5. Required correction at the two boundary points

For the reference rectilinear projection, a ray at half-angle θ has an ideal image radius:

$$
R_{ideal}(\theta)=f\tan\theta.
$$

The required multiplicative correction is:

$$
C=\frac{R_{ideal}}{R_{observed}}.
$$

At the side midpoint and corner:

$$
C_s=\frac{f\tan(FOV_s/2)}{a},\qquad
C_d=\frac{f\tan(FOV_d/2)}{a\sqrt{2}}.
$$

Explicit form when the input FOV values are in degrees and `tan()` expects radians:

```text
theta_s = FOV_s * pi / 360
theta_d = FOV_d * pi / 360
Cs = f * tan(theta_s) / a
Cd = f * tan(theta_d) / (a * sqrt(2))
```

These estimates assume that the specified FOV boundaries correspond to the stated physical radii. They are not measurements of the correction throughout the image.

### 6. Deriving K1 and K2

Evaluate the correction model at q=1 and q=2:

$$
1+K_1+K_2=C_s,
$$

$$
1+2K_1+4K_2=C_d.
$$

Solving the two equations gives:

$$
\boxed{K_2=\frac{C_d-2C_s+1}{2}},\qquad
\boxed{K_1=C_s-1-K_2}.
$$

Equivalently:

```text
K2 = (Cd - 2*Cs + 1) / 2
K1 = Cs - 1 - K2
// Equivalent: K1 = (4*Cs - Cd - 3) / 2
```

#### PAJ7025R3 example

```text
f     = 0.378 mm
a     = 0.539 mm
FOV_s = 111.3°  → theta_s = 55.65°
FOV_d = 140.0°  → theta_d = 70.00°

Cs = 1.0261407311
Cd = 1.3624550049

K2 =  0.1550867713
K1 = -0.1289460402
```

Rounded profile values:

| Field | Initial R3 value | Basis |
|---|---:|---|
| `lensRadialK1` | `-0.128946f` | Nominal two-point fit |
| `lensRadialK2` | `0.155087f` | Nominal two-point fit |
| `lensCorrectionMax` | `1.5f` | Chosen upper safety limit |
| `lensCorrectionMin` | `0.8f` | Retained lower safety limit |

The last four arguments of the current `MakeProfile()` interface are ordered **K1, K2, Max, Min**:

```cpp
// Theoretical inverse radial fit; not a measured lens calibration.
-0.128946f, // lensRadialK1
 0.155087f, // lensRadialK2
 1.5f,      // lensCorrectionMax
 0.8f       // lensCorrectionMin
```

A negative K1 does not, by itself, mean the correction is wrong: the complete polynomial matters. In this interpolant, the factor falls slightly below 1 in the inner region, then rises towards the corners. This does not prove that the real lens follows exactly that curve.

### 7. Why “Distortion −30%” is not enough

For radial geometric distortion defined as a fraction:

$$
D=\frac{R_{observed}-R_{ideal}}{R_{ideal}},
$$

the inverse correction factor at the same point is:

$$
C=\frac{1}{1+D}.
$$

Thus −30% means D = −0.30 and would imply C = 1/0.70 ≈ 1.42857, **not** 1.30 and **not** K1 = 0.30. However, the location and definition of that distortion value must be known. A single value cannot determine two coefficients or the complete curve. Geometric distortion and TV distortion are also different quantities. [Edmund Optics — Distortion](https://www.edmundoptics.com/knowledge-center/application-notes/imaging/distortion/).

The R3 diagonal factor obtained from the nominal focal length and FOV, Cd ≈ 1.362455, corresponds to D ≈ −26.60%. This is of the same order as −30%, but it is not identical. The specifications have not been treated as simultaneous exact constraints: focal length and FOV have tolerances, and the excerpt does not provide the detailed distortion definition.

### 8. Choosing Min and Max

`lensCorrectionMin` and `lensCorrectionMax` **cannot be uniquely calculated as optical properties from the datasheet**. They are limits on the multiplicative factor, particularly useful when reconstructed LEDs fall outside the sensor field.

To avoid clipping the nominal model over q ∈ [0,2], find the minimum and maximum of:

$$
C(q)=1+K_1q+K_2q^2.
$$

Evaluate C at:

- q=0;
- q=2;
- q*=−K1/(2K2), if K2 is nonzero and q* lies in [0,2].

The smallest and largest values are C_min and C_max. If K2=0, evaluating the interval endpoints is sufficient.

Limits that leave the model unchanged within the nominal field must satisfy:

$$
0<\text{lensCorrectionMin}\le C_{min},\qquad
\text{lensCorrectionMax}\ge C_{max}.
$$

Any additional margin is an engineering choice, not a unique result of the optical calculation.

For the R3 candidate:

```text
q*    ≈ 0.415722
C_min ≈ 0.973197
C_max ≈ 1.362455
```

Therefore:

- **Min=0.8** retains the previous lower bound and is below the nominal minimum;
- **Max=1.5** leaves approximately 10% margin above the nominal maximum;
- neither bound clips the R3 candidate inside the nominal sensor field;
- the previous Max=1.2 would clip the correction before reaching the corners.

The value 0.8 was not derived using a 10% margin: it is a retained design choice. Likewise, 1.5 is not “the lens distortion”; it is a factor limit. These bounds neither make off-sensor extrapolation exact nor guarantee valid quadrilateral geometry.

### 9. Checking the radial map

Before clamping, the corrected normalised radius is:

$$
r'=r(1+K_1r^2+K_2r^4).
$$

Within the intended domain, verify that the factor is positive and that:

$$
\frac{dr'}{dr}=1+3K_1r^2+5K_2r^4>0.
$$

This prevents the radial map from reversing radial order. For the R3 candidate over r ∈ [0,√2], the minimum derivative is approximately 0.95175, so this check passes. It does not establish agreement with the real lens.

### 10. Why retain K1=K2=0 for the R2

Mechanically applying the nominal R2 numbers produces:

```text
Cs ≈ 0.9562865381
Cd ≈ 0.9539441098
K1 ≈ -0.0643989788
K2 ≈  0.0206855169
```

These numbers are **not a proposed R2 profile change**. Under the geometric convention above, they imply approximately +4.57% distortion at a side midpoint and +4.83% at a corner, whereas the datasheet excerpt specifies <2.8% without further detail.

This discrepancy shows why nominal focal length/FOV values with tolerances must not be treated as an exact calibration. It proves neither that the datasheet is wrong nor that the lens has zero distortion. Given the satisfactory real-world behaviour already reported for the R2, retain:

```cpp
0.0f, // lensRadialK1
0.0f, // lensRadialK2
1.2f, // lensCorrectionMax
0.8f  // lensCorrectionMin
```

With K1=K2=0, the correction function returns immediately: the bounds are not used. Nonzero R2 coefficients should be justified by measurements, not solely by this nominal interpolation.

### 11. Refinement using measurements

If more measured pairs of observed radius and known angle become available, compute for each point i:

```text
Ci = f * tan(theta_i) / R_observed_i
qi = squared normalised radius, using the firmware's normalisation
```

Use nonzero observed radii; the centre has a limiting factor of 1 in this model. Fit:

$$
C_i-1\simeq K_1q_i+K_2q_i^2
$$

using least squares, preferably solved with QR/SVD. Check residual errors, optical centre and radial monotonicity. Measurements distributed across the field provide more information than two FOV boundary values.

Validation should include the centre and edges, several working distances and camera rotations, and the complete tracking/calibration pipeline. A numerically correct correction function does not by itself establish the accuracy of the lens model or missing-LED reconstruction.

All trigonometric and optical calculations in this document are **offline calculations used to choose profile constants**. They add no calculations, waits or workload to the 209 Hz firmware loop. The R3 coefficients remain theoretical candidates until confirmed with the actual lens.

---

<a id="italiano"></a>

## Italiano

### Scopo e ambito

Questo documento spiega come stimare valori iniziali per `lensRadialK1` e `lensRadialK2` dai dati ottici nominali e come scegliere `lensCorrectionMin` e `lensCorrectionMax` senza troncare la correzione prevista nel campo del sensore.

La derivazione riguarda specificamente **PAJ7025R2 e PAJ7025R3**, che hanno sensori quadrati e campi visivi orizzontale e verticale uguali.

**Sono valori teorici di partenza, non una calibrazione misurata della lente.** Focale e campi visivi ai bordi non descrivono l'intera curva di distorsione. Due coefficienti possono interpolare esattamente due punti nominali, ma ciò non garantisce la precisione fra quei punti o sul singolo esemplare. Resta necessaria la verifica sulla camera reale.

I dati numerici provengono dalle specifiche PixArt fornite. Le formule riportate sono una derivazione per il modello di correzione OpenFIRE, non coefficienti ufficiali di calibrazione PixArt.

### 1. Dati necessari

| Voce nelle specifiche | Simbolo | PAJ7025R2 | PAJ7025R3 | Impiego |
|---|---|---:|---:|---|
| Effective Focal Length | `f` | 1,484307 mm | 0,378 mm | Raggio dell'immagine rettilineare di riferimento |
| Sensor Pixel Resolution | `N` | 98 × 98 | 98 × 98 | Dimensioni fisiche dell'area sensibile |
| Sensor Pixel Size | `p` | 11 × 11 µm | 11 × 11 µm | Conversione pixel fisici → mm |
| Angle Field of View — Horizontal/Vertical | `FOV_s` | 38,3° | 111,3° | Vincolo al centro di un bordo |
| Angle Field of View — Diagonal | `FOV_d` | 52,2° | 140° | Vincolo all'angolo |
| Note Y ai bordi | `a`, `d` | 0,539 / 0,762 mm | 0,539 / 0,762 mm | Verifica delle semidimensioni fisiche |
| Distortion | `D` | < 2,8% | −30% | Controllo indicativo, non una curva completa |

I FOV sono angoli **totali**. Nelle formule si usa metà di ogni angolo, convertendo i gradi in radianti quando richiesto dalla funzione trigonometrica.

La risoluzione virtuale 4096 × 4096 non è la matrice fisica del sensore. Non moltiplicare 4096 per 11 µm: usare i **98 pixel fisici**.

F-number, guadagno, esposizione, frame rate, interfaccia SPI e consumo non entrano in questa stima geometrica. Alcuni influenzano la qualità del rilevamento, ma non determinano K1/K2. Image Circle, Back Focal Length e numero di elementi della lente non sostituiscono la focale efficace o le dimensioni dell'area sensibile.

### 2. Ipotesi

La derivazione R2/R3 assume:

- sensore quadrato, pixel quadrati e stesso FOV orizzontale/verticale;
- centro ottico coincidente con il centro utilizzato dal firmware;
- coordinate restituite proporzionali alla posizione nell'area sensibile, senza una rettificazione DSP già applicata;
- distorsione approssimabile come radiale, senza termini tangenziali o di decentramento;
- focale efficace nominale usata come focale paraassiale della proiezione rettilineare di riferimento;
- valori nominali dei FOV e delle dimensioni riferiti alla stessa area attiva.

OpenFIRE normalizza ciascun asse rispetto alla propria semidimensione. La derivazione semplice qui sotto è appropriata per queste PAJ quadrate; non va trasferita automaticamente a sensori rettangolari, immagini ritagliate o una diversa normalizzazione delle coordinate.

### 3. Modello di correzione usato dal firmware

Siano x e y le coordinate osservate e cx/cy il centro:

```text
cx = mouseResX / 2
cy = mouseResY / 2

nx = (x - cx) / cx
ny = (y - cy) / cy
q  = nx*nx + ny*ny                  // q = r²

C(q) = 1 + K1*q + K2*q*q           // fattore prima del clamp
C_limited = clamp(C(q), Min, Max)

x_corrected = cx + (x - cx)*C_limited
y_corrected = cy + (y - cy)*C_limited
```

I campi corrispondono a:

| Campo del firmware | Significato |
|---|---|
| `lensRadialK1` | Coefficiente di q, cioè r² |
| `lensRadialK2` | Coefficiente di q², cioè r⁴ |
| `lensCorrectionMax` | Limite superiore del fattore moltiplicativo di correzione |
| `lensCorrectionMin` | Limite inferiore del fattore moltiplicativo di correzione |

Per le PAJ quadrate, ai bordi idealizzati del campo:

| Posizione | r | q |
|---|---:|---:|
| Centro | 0 | 0 |
| Centro di un bordo | 1 | 1 |
| Angolo | √2 | 2 |

La correzione agisce su **entrambi gli assi**, rispetto al centro ottico. LED e punto di mira devono usare lo stesso spazio di coordinate corretto; lavorare su copie locali evita che correzioni ripetute si accumulino nei valori di calibrazione memorizzati.

È una mappa **osservato → corretto**, normalizzata alla semidimensione del sensore. Non è automaticamente compatibile con i coefficienti di un programma di calibrazione. Per esempio, il modello standard OpenCV applica la distorsione alle coordinate ideali normalizzate rispetto alla focale: direzione e normalizzazione vanno convertite prima di riusare i suoi coefficienti. [Documentazione OpenCV](https://docs.opencv.org/4.13.0/d9/d0c/group__calib3d.html).

Non basta quindi invertire il segno di coefficienti ottenuti con un altro modello, specialmente per distorsioni ampie.

### 4. Dimensioni fisiche del sensore

Usando millimetri:

```text
p_mm = 11 / 1000             = 0.011 mm
a    = N*p_mm/2 = 98*0.011/2 = 0.539 mm
d    = a*sqrt(2)             = 0.7622611101 mm
```

a è il semilato e d è la semidiagonale. Il dato 0,762 mm della scheda è coerente con la semidiagonale arrotondata. I coefficienti proposti usano d = a√2, evitando quell'arrotondamento aggiuntivo.

### 5. Correzione richiesta nei due punti di bordo

Per la proiezione rettilineare di riferimento, un raggio a semiangolo θ ha un raggio immagine ideale:

$$
R_{ideale}(\theta)=f\tan\theta.
$$

Il fattore moltiplicativo richiesto è:

$$
C=\frac{R_{ideale}}{R_{osservato}}.
$$

Al centro del bordo e all'angolo:

$$
C_s=\frac{f\tan(FOV_s/2)}{a},\qquad
C_d=\frac{f\tan(FOV_d/2)}{a\sqrt{2}}.
$$

Forma esplicita con FOV in ingresso espressi in gradi e `tan()` che richiede radianti:

```text
theta_s = FOV_s * pi / 360
theta_d = FOV_d * pi / 360
Cs = f * tan(theta_s) / a
Cd = f * tan(theta_d) / (a * sqrt(2))
```

Queste stime assumono che i bordi del FOV dichiarato corrispondano ai raggi fisici indicati. Non sono misure della correzione in ogni posizione dell'immagine.

### 6. Ricavare K1 e K2

Valutando il modello di correzione nei punti q=1 e q=2:

$$
1+K_1+K_2=C_s,
$$

$$
1+2K_1+4K_2=C_d.
$$

Risolvendo le due equazioni:

$$
\boxed{K_2=\frac{C_d-2C_s+1}{2}},\qquad
\boxed{K_1=C_s-1-K_2}.
$$

Equivalentemente:

```text
K2 = (Cd - 2*Cs + 1) / 2
K1 = Cs - 1 - K2
// Equivalent: K1 = (4*Cs - Cd - 3) / 2
```

#### Esempio PAJ7025R3

```text
f     = 0.378 mm
a     = 0.539 mm
FOV_s = 111.3°  → theta_s = 55.65°
FOV_d = 140.0°  → theta_d = 70.00°

Cs = 1.0261407311
Cd = 1.3624550049

K2 =  0.1550867713
K1 = -0.1289460402
```

Valori arrotondati del profilo:

| Campo | Valore iniziale R3 | Origine |
|---|---:|---|
| `lensRadialK1` | `-0.128946f` | Interpolazione nominale a due punti |
| `lensRadialK2` | `0.155087f` | Interpolazione nominale a due punti |
| `lensCorrectionMax` | `1.5f` | Limite superiore di sicurezza scelto |
| `lensCorrectionMin` | `0.8f` | Limite inferiore di sicurezza conservato |

Gli ultimi quattro argomenti dell'interfaccia attuale di `MakeProfile()` sono ordinati **K1, K2, Max, Min**:

```cpp
// Theoretical inverse radial fit; not a measured lens calibration.
-0.128946f, // lensRadialK1
 0.155087f, // lensRadialK2
 1.5f,      // lensCorrectionMax
 0.8f       // lensCorrectionMin
```

Il K1 negativo non indica da solo che la correzione sia sbagliata: conta il polinomio completo. In questo interpolante il fattore scende leggermente sotto 1 nella zona interna e poi aumenta verso gli angoli. Non è una prova che la lente reale segua esattamente quella curva.

### 7. Perché non basta “Distortion −30%”

Per la distorsione geometrica radiale definita come frazione:

$$
D=\frac{R_{osservato}-R_{ideale}}{R_{ideale}},
$$

il fattore inverso nello stesso punto è:

$$
C=\frac{1}{1+D}.
$$

Quindi −30% significa D = −0,30 e darebbe C = 1/0,70 ≈ 1,42857, **non** 1,30 e **non** K1 = 0,30. Bisogna però conoscere la posizione e la definizione di quel valore di distorsione. Un singolo valore non determina due coefficienti o l'intera curva. Distorsione geometrica e distorsione TV sono inoltre grandezze differenti. [Edmund Optics — Distortion](https://www.edmundoptics.com/knowledge-center/application-notes/imaging/distortion/).

Il fattore diagonale R3 ottenuto da focale e FOV nominali, Cd ≈ 1,362455, corrisponde a D ≈ −26,60%. È dello stesso ordine del −30%, ma non coincide. Non abbiamo trattato tutte le specifiche come vincoli esatti simultanei: focale e FOV hanno tolleranze e l'estratto non fornisce la definizione dettagliata della distorsione.

### 8. Scegliere Min e Max

`lensCorrectionMin` e `lensCorrectionMax` **non sono caratteristiche ottiche ricavabili univocamente dalla scheda**. Sono limiti del fattore moltiplicativo, particolarmente utili quando i LED ricostruiti si trovano fuori dal campo del sensore.

Per non troncare il modello nominale nell'intervallo q ∈ [0,2], trovare minimo e massimo di:

$$
C(q)=1+K_1q+K_2q^2.
$$

Valutare C in:

- q=0;
- q=2;
- q*=−K1/(2K2), se K2 è diverso da zero e q* appartiene a [0,2].

Il valore più piccolo e quello più grande sono C_min e C_max. Se K2=0 è sufficiente valutare gli estremi dell'intervallo.

Limiti che lasciano invariato il modello nel campo nominale devono soddisfare:

$$
0<\text{lensCorrectionMin}\le C_{min},\qquad
\text{lensCorrectionMax}\ge C_{max}.
$$

Un margine aggiuntivo è una scelta progettuale, non un risultato univoco del calcolo ottico.

Per la candidata R3:

```text
q*    ≈ 0.415722
C_min ≈ 0.973197
C_max ≈ 1.362455
```

Di conseguenza:

- **Min=0,8** mantiene il limite inferiore precedente ed è sotto il minimo nominale;
- **Max=1,5** lascia circa il 10% di margine sopra il massimo nominale;
- nessuno dei due limiti tronca la candidata R3 nel campo nominale del sensore;
- il precedente Max=1,2 taglierebbe la correzione prima di raggiungere gli angoli.

Il valore 0,8 non è stato ricavato applicando un margine del 10%: è una scelta progettuale conservata. Allo stesso modo, 1,5 non è “la distorsione della lente”: è un limite al fattore. Questi limiti non rendono esatta l'extrapolazione fuori sensore e non garantiscono la validità geometrica dei quadrilateri.

### 9. Controllare la mappa radiale

Prima del clamp, il raggio corretto normalizzato è:

$$
r'=r(1+K_1r^2+K_2r^4).
$$

Nel dominio di utilizzo verificare che il fattore sia positivo e che:

$$
\frac{dr'}{dr}=1+3K_1r^2+5K_2r^4>0.
$$

Questo evita che la mappa radiale inverta l'ordine dei raggi. Per la candidata R3, nell'intervallo r ∈ [0,√2], il minimo della derivata è circa 0,95175: il controllo è soddisfatto. Non dimostra la corrispondenza con la lente reale.

### 10. Perché mantenere K1=K2=0 sulla R2

Applicando meccanicamente i numeri nominali R2 si ottiene:

```text
Cs ≈ 0.9562865381
Cd ≈ 0.9539441098
K1 ≈ -0.0643989788
K2 ≈  0.0206855169
```

Questi numeri **non sono una proposta di modifica del profilo R2**. Con la convenzione geometrica precedente implicherebbero circa +4,57% di distorsione al centro di un bordo e +4,83% all'angolo, mentre l'estratto della scheda riporta <2,8% senza ulteriori dettagli.

Questa discrepanza mostra perché focale/FOV nominali con tolleranze non vadano trattati come una calibrazione esatta. Non dimostra né che la scheda sia sbagliata né che la lente abbia distorsione nulla. Dato il funzionamento reale già soddisfacente riportato per la R2, mantenere:

```cpp
0.0f, // lensRadialK1
0.0f, // lensRadialK2
1.2f, // lensCorrectionMax
0.8f  // lensCorrectionMin
```

Con K1=K2=0 la funzione di correzione esce subito: i limiti non vengono utilizzati. Coefficienti R2 non nulli andrebbero giustificati da misure, non soltanto da questa interpolazione nominale.

### 11. Affinamento tramite misure

Se diventano disponibili più coppie misurate di raggio osservato e angolo noto, per ogni punto i calcolare:

```text
Ci = f * tan(theta_i) / R_observed_i
qi = raggio normalizzato al quadrato, con la normalizzazione del firmware
```

Usare raggi osservati non nulli; al centro il fattore limite di questo modello è 1. Adattare:

$$
C_i-1\simeq K_1q_i+K_2q_i^2
$$

con un fit ai minimi quadrati, preferibilmente risolto con QR/SVD. Verificare errori residui, centro ottico e monotonicità radiale. Misure distribuite nel campo forniscono più informazioni di due soli valori di FOV ai bordi.

La verifica dovrebbe includere centro e bordi, più distanze operative e rotazioni della camera, nonché l'intera catena di tracking e calibrazione. Una funzione di correzione numericamente corretta non dimostra da sola la precisione del modello ottico o della ricostruzione dei LED mancanti.

Tutti i calcoli trigonometrici e ottici di questo documento sono **calcoli fuori linea per scegliere le costanti del profilo**. Non aggiungono calcoli, attese o carico al ciclo firmware a 209 Hz. I coefficienti R3 rimangono candidati teorici fino alla conferma con la lente reale.
