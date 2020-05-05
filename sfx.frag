#version 130
const float PI = radians(180.);
const float TAU = 2.*PI;
float clip(float a) { return clamp(a,-1.,1.); }
float smstep(float a, float b, float x) {return smoothstep(a, b, clamp(x, a, b));}
float theta(float x) { return smstep(0.,1e-3,x); }
float _sin(float a) { return sin(TAU * mod(a,1.)); }
float _sin_(float a, float p) { return sin(TAU * mod(a,1.) + p); }
float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float freqC1(float note){ return 32.7 * exp2(note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }
float minus1hochNminus1halbe(int n) { return sin(.5*PI*float(n)); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }
float fhelp(float x) { return 1. + .333*x; } // 1. + .33333*x + .1*x*x + .02381*x*x*x + .00463*x*x*x*x;
float s_atan(float a) { return .636 * atan(a); }
float doubleslope(float t, float a, float d, float s) { return smstep(-.00001,a,t) - (1.-s) * smstep(0.,d,t-a); }
float s_moothmin(float a, float k) {
    float ha = max(1.-2.*abs(abs(a)-1.), 0.);
    return a >= 0. ? min(a, 1.) - .5/6.*ha*ha*ha : max(a, -1.) + .5/6.*ha*ha*ha;
}
float s_moothmin(float a) { return s_moothmin(a,.5); }

#define SONGLENGTH 33.5667
#define NTIME 2
const float pos_B[2] = float[2](0.,20.);
const float pos_t[2] = float[2](0.,33.3333);
const float pos_BPS[1] = float[1](.6);
const float pos_SPB[1] = float[1](1.6667);
float BPS, SPB, BT;

float Tsample;

#define filterthreshold 1.e-3

//TEXCODE

float drop_phase(float time, float t1, float f0, float f1)
{
    float t = min(time, t1);
    float phi = f0*t + .5*(f1-f0)/t1*t*t;

    if(time > t1)
    {
        phi += f1 * (time - t1);
    }
    return phi;
}

float lpnoise(float t, float fq)
{
    t *= fq;
    float tt = fract(t);
    float tn = t - tt;
    return mix(pseudorandom(floor(tn) / fq), pseudorandom(floor(tn + 1.0) / fq), smstep(0.0, 1.0, tt));
}

float reverb_phase(float t, float amt)
{
    float r = lpnoise(t, 100.0) + 0.2*lpnoise(t, 550.0) + 0.1*lpnoise(t, 1050.0)*exp(-5.*t);
    return amt * r;
}

float env_AHDSR(float x, float L, float A, float H, float D, float S, float R)
{
    return (x<A ? x/A : x<A+H ? 1. : x<A+H+D ? (1. - (1.-S)*(x-H-A)/D) : x<=L-R ? S : x<=L ? S*(L-x)/R : 0.);
}

float sinshape(float x, float amt, float parts)
{
    return (1.-amt) * x + amt * sign(x) * 0.5 * (1. - cos(parts*PI*x));
}

float comp_SAW(int N, float inv_N, float PW) {return inv_N * (1. - _sin(float(N)*PW));}
float comp_TRI(int N, float inv_N, float PW) {return N % 2 == 0 ? .1 * inv_N * _sin(float(N)*PW) : inv_N * inv_N * (1. - _sin(float(N)*PW));}
float comp_SQU(int N, float inv_N, float PW) {return inv_N * (minus1hochN(N) * _sin(.5*float(N)*PW + .25) - 1.);}
float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}
float comp_OBO(int N, float inv_N, float PW) {return sqrt(inv_N) * (1. + _sin(float(N)*(1.5+PW)+.5*PI));}

float MADD(float t, float f, float p0, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, float LOWCUT, float keyF)
{
    float ret = 0.;
    float f_ = keyF > .99 ? 1. : (keyF < 1.e-3 ? f : pow(f, 1.-keyF));
    float INR = f_/CO;
    float IRESQ = 1./(RES_Q*f_);

    float p = f*t;
    float float_N, inv_N, comp_mix, filter_N;
    for(int N = 1 + int(LOWCUT/f - 1.e-3); N<=NMAX; N+=NINC)
    {
        float_N = float(N);
        inv_N = 1./float_N;
        comp_mix = MIX < -1. ? (MIX+2.) * comp_SAW(N,inv_N,PW)  - (MIX+1.) * comp_OBO(N,inv_N,PW)
                 : MIX <  0. ? (MIX+1.) * comp_TRI(N,inv_N,PW)  -     MIX  * comp_SAW(N,inv_N,PW)
                 : MIX < 1. ? (1.-MIX) * comp_TRI(N,inv_N,PW)  +     MIX  * comp_SQU(N,inv_N,PW)
                            : (MIX-1.) * comp_HAE(N,inv_N,PW)  + (2.-MIX) * comp_SQU(N,inv_N,PW);

        if(abs(comp_mix) < 1e-4) continue;

        filter_N = pow(1. + pow(float_N*INR,NDECAY),-.5) + RES * exp(-pow((float_N*f-CO)*IRESQ,2.));

        ret += comp_mix * filter_N * (_sin_(float_N * p, p0) + _sin_(float_N * p * (1.+DET), p0));
    }
    return s_moothmin(ret);
}

float MADD(float t, float f, float p0, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, int keyF)
{
    return MADD(t, f, p0, NMAX, NINC, MIX, CO, NDECAY, RES, RES_Q, DET, PW, 0., keyF);
}

float protokick(float t, float f_start, float f_end, float fdecay, float hold, float decay, float drive, float detune, float rev_amount, float rev_hold, float rev_decay, float rev_drive)
{
    float phi = drop_phase(t, fdecay, f_start, f_end);
    float rev_phi = phi + reverb_phase(t, 1.);
    return clamp(drive*.5*(_sin(phi)+_sin((1.-detune)*phi)),-1.,1.) * exp(-max(t-hold, 0.)/decay)
         + rev_amount*clamp(rev_drive*.5*(_sin(rev_phi)+_sin((1.-detune)*rev_phi)),-1.,1.) * exp(-max(t-rev_hold, 0.)/rev_decay);
}

float maceboss_vol(float _BEAT)
{
    return _BEAT<0 ? 0. : 1.;
}

uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iTexSize;
uniform sampler2D iSequence;
uniform float iSequenceWidth;

// Read short value from texture at index off
float rshort(in float off)
{
    float hilo = mod(off, 2.);
    off = .5*off;
    vec2 ind = vec2(mod(off, iSequenceWidth), floor(off/iSequenceWidth));
    vec4 block = texelFetch(iSequence, ivec2(ind), 0);
    vec2 data = mix(block.rg, block.ba, hilo);
    return round(dot(vec2(255., 65280.), data));
}

// Read float value from texture at index off
float rfloat(int off)
{
    float d = rshort(float(off));
    float sign = floor(d/32768.),
        exponent = floor(d*9.765625e-4 - sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    if(exponent == 0.)
         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
}

#define NTRK 3
#define NMOD 19
#define NPTN 3
#define NNOT 62
#define NDRM 52

int trk_sep(int index)      {return int(rfloat(index));}
int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);}
float trk_pre(int index)    {return     rfloat(index+1+4*NTRK);}
float trk_slide(int index)  {return     rfloat(index+1+5*NTRK);} // idea for future: change to individual note_slide_time
float mod_on(int index)     {return     rfloat(index+1+6*NTRK);}
float mod_off(int index)    {return     rfloat(index+1+6*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloat(index+1+6*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloat(index+1+6*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloat(index+1+6*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_pan(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+3*NNOT);}
float note_vel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+4*NNOT);}
float note_slide(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+5*NNOT);}
float note_aux(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+6*NNOT);}
float drum_rel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+7*NNOT);}

vec2 mainSynth(float time)
{
    float sL = 0.;
    float sR = 0.;
    float dL = 0.;
    float dR = 0.;

    if (time > SONGLENGTH) return vec2(0.);
    
    int _it;
    for(_it = 0; _it < NTIME - 2 && pos_t[_it + 1] < time; _it++);
    BPS = pos_BPS[_it];
    SPB = pos_SPB[_it];
    BT = pos_B[_it] + (time - pos_t[_it]) * BPS;

    float time2 = time - .0002;
    float sidechain = 1.;

    float amaysynL, amaysynR, amaydrumL, amaydrumR, B, Bon, Boff, Bprog, Bproc, L, tL, _t, _t2, vel, rel, pre, f, amtL, amtR, env, slide, aux;
    int tsep0, tsep1, _modU, _modL, ptn, psep0, psep1, _noteU, _noteL, syn, drum;

    for(int trk = 0; trk < NTRK; trk++)
    {
        tsep0 = trk_sep(trk);
        tsep1 = trk_sep(trk + 1);

        syn = trk_syn(trk);
        rel = trk_rel(trk) + 1.e-3;
        pre = trk_pre(trk);

        for(_modU = tsep0; (_modU < tsep1 - 1) && (BT > mod_on(_modU + 1) - pre); _modU++);
        for(_modL = tsep0; (_modL < tsep1 - 1) && (BT >= mod_off(_modL) + rel); _modL++);

        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            B = BT - mod_on(_mod) + pre;

            ptn   = mod_ptn(_mod);
            psep0 = ptn_sep(ptn);
            psep1 = ptn_sep(ptn + 1);

            for(_noteU = psep0; (_noteU < psep1 - 1) && (B > note_on(_noteU + 1)); _noteU++);
            for(_noteL = psep0; (_noteL < psep1 - 1) && (B >= note_off(_noteL) + rel); _noteL++);

            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                if(syn == 139)
                {
                    drum = int(note_pitch(_note));
                    rel = drum_rel(drum) + 1.e-3;
                }

                amaysynL  = 0.;
                amaysynR  = 0.;
                amaydrumL = 0.;
                amaydrumR = 0.;

                Bon   = note_on(_note);
                Boff  = note_off(_note) + rel;
                L     = Boff - Bon;
                tL    = L * SPB;
                Bprog = max(0., B - Bon); // I DO NOT GET THIS WEIRD FIX, but Revision is approaching
                Bproc = Bprog / L;
                _t    = Bprog * SPB; 
                _t2   = _t - .0002; // this is on purpose not max(0., _t - .0002), because I hope future-QM is clever
                vel   = note_vel(_note);
                amtL  = clamp(1. - note_pan(_note), 0., 1.);
                amtR  = clamp(1. + note_pan(_note), 0., 1.);
                slide = note_slide(_note);
                aux   = note_aux(_note);

                if(syn == 139)
                {
                    env = trk_norm(trk) * theta(Bprog) * theta(L - Bprog);
                    if(drum == 0) { sidechain = min(sidechain, 1. - vel * (clamp(1.e4 * Bprog,0.,1.) - pow(Bprog/(L-rel),8.)));}
                    else if(drum == 5){
                        amaydrumL = vel*.9*protokick(_t,242.,55.,.036,.03,.0666,1.42,.02,.25,.01,.1,.4)
      +.9*protokick(_t,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.)
      +.64*((clamp(2.27*_tri(drop_phase(_t,.03,241.,72.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t-.01))+.91*clamp(.9*_tri(drop_phase(_t,.03,241.,72.)+.91*lpnoise(_t,8164.)),-1.,1.)*exp(-20.76*_t)+.05*lpnoise(_t,10466.)*(1.-smstep(0.,.18,_t-.56))+.56*lpnoise(_t,7123.)*exp(-_t*5.45)+.11*lpnoise(_t,1134.)*exp(-_t*13.82))*smstep(0.,.004,_t));
                        amaydrumR = vel*.9*protokick(_t2,242.,55.,.036,.03,.0666,1.42,.02,.25,.01,.1,.4)
      +.9*protokick(_t2,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.)
      +.64*((clamp(2.27*_tri(drop_phase(_t2,.03,241.,72.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t2-.01))+.91*clamp(.9*_tri(drop_phase(_t2,.03,241.,72.)+.91*lpnoise(_t2,8164.)),-1.,1.)*exp(-20.76*_t2)+.05*lpnoise(_t2,10466.)*(1.-smstep(0.,.18,_t2-.56))+.56*lpnoise(_t2,7123.)*exp(-_t2*5.45)+.11*lpnoise(_t2,1134.)*exp(-_t2*13.82))*smstep(0.,.004,_t2));
                    }
                    else if(drum == 7){
                        amaydrumL = vel*.85*(clamp(1.15*_tri(drop_phase(_t,.13,157.,76.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t-.13))+.81*clamp(.24*_tri(drop_phase(_t,.13,157.,76.)+.81*lpnoise(_t,2401.)),-1.,1.)*exp(-14.8*_t)+.01*lpnoise(_t,4079.)*(1.-smstep(0.,.7,_t-.12))+.5*lpnoise(_t,5164.)*exp(-_t*19.79)+.76*lpnoise(_t,8446.)*exp(-_t*24.))*smstep(0.,.002,_t);
                        amaydrumR = vel*.85*(clamp(1.15*_tri(drop_phase(_t2,.13,157.,76.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t2-.13))+.81*clamp(.24*_tri(drop_phase(_t2,.13,157.,76.)+.81*lpnoise(_t2,2401.)),-1.,1.)*exp(-14.8*_t2)+.01*lpnoise(_t2,4079.)*(1.-smstep(0.,.7,_t2-.12))+.5*lpnoise(_t2,5164.)*exp(-_t2*19.79)+.76*lpnoise(_t2,8446.)*exp(-_t2*24.))*smstep(0.,.002,_t2);
                    }
                    else if(drum == 10){
                        amaydrumL = vel*(lpnoise(_t,10000.)*smstep(0.,.01,_t)*(1.-(1.-.13)*smstep(0.,.12,_t-.01))-.3*(1.00*lpnoise((_t-0.00),10000.)*smstep(0.,.01,(_t-0.00))*(1.-(1.-.13)*smstep(0.,.12,(_t-0.00)-.01))+.61*lpnoise((_t-1.20e-03),10000.)*smstep(0.,.01,(_t-1.20e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-1.20e-03)-.01))+.372*lpnoise((_t-2.40e-03),10000.)*smstep(0.,.01,(_t-2.40e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-2.40e-03)-.01))));
                        amaydrumR = vel*(lpnoise(_t,10000.)*smstep(0.,.01,_t)*(1.-(1.-.13)*smstep(0.,.12,_t-.01))-.3*(1.00*lpnoise((_t-0.00),10000.)*smstep(0.,.01,(_t-0.00))*(1.-(1.-.13)*smstep(0.,.12,(_t-0.00)-.01))+.61*lpnoise((_t-1.20e-03),10000.)*smstep(0.,.01,(_t-1.20e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-1.20e-03)-.01))+.372*lpnoise((_t-2.40e-03),10000.)*smstep(0.,.01,(_t-2.40e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-2.40e-03)-.01))));
                    }
                    else if(drum == 12){
                        amaydrumL = vel*1.2*fract(sin(_t*100.*.5)*50000.*.5)*doubleslope(_t,0.,.03,.1)*exp(-13.*Bprog);
                        amaydrumR = vel*1.2*fract(sin(_t2*100.*.5)*50000.*.5)*doubleslope(_t2,0.,.03,.1)*exp(-13.*Bprog);
                    }
                    else if(drum == 21){
                        amaydrumL = vel*((clamp(1.09*_tri(drop_phase(_t,.08,249.,77.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t-.04))+.97*clamp(.99*_tri(drop_phase(_t,.08,249.,77.)+.97*lpnoise(_t,9855.)),-1.,1.)*exp(-21.22*_t)+.03*lpnoise(_t,10655.)*(1.-smstep(0.,.58,_t-.81))+.71*lpnoise(_t,7520.)*exp(-_t*16.22)+.57*lpnoise(_t,4386.)*exp(-_t*29.48))*smstep(0.,.005,_t));
                        amaydrumR = vel*((clamp(1.09*_tri(drop_phase(_t2,.08,249.,77.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t2-.04))+.97*clamp(.99*_tri(drop_phase(_t2,.08,249.,77.)+.97*lpnoise(_t2,9855.)),-1.,1.)*exp(-21.22*_t2)+.03*lpnoise(_t2,10655.)*(1.-smstep(0.,.58,_t2-.81))+.71*lpnoise(_t2,7520.)*exp(-_t2*16.22)+.57*lpnoise(_t2,4386.)*exp(-_t2*29.48))*smstep(0.,.005,_t2));
                    }
                    
                    if(drum > 0)
                    {
                        dL += amtL * s_moothmin(env * amaydrumL);
                        dR += amtR * s_moothmin(env * amaydrumR);
                    }
                }
                else
                {
                    f = freqC1(note_pitch(_note) + mod_transp(_mod));

                    if(abs(slide) > 1e-3) // THIS IS SLIDEY BIZ
                    {
                        float Bslide = trk_slide(trk);
                        float fac = slide * log(2.)/12.;
                        if (Bprog <= Bslide)
                        {
                            float help = 1. - Bprog/Bslide;
                            f *= Bslide * (fhelp(fac) - help * fhelp(fac*help*help)) / Bprog;
                        }
                        else
                        {
                            f *= 1. + (Bslide * (fhelp(fac)-1.)) / Bprog;
                        }
                    }

                    env = theta(Bprog) * (1. - smstep(Boff-rel, Boff, B));
                    if(syn == 0){amaysynL = _sin(f*_t); amaysynR = _sin(f*_t2);}
                    else if(syn == 113){
                        time2 = time-0.007; _t2 = _t-0.007;
                        amaysynL = .9*vel*(.4*sinshape(s_atan(2.*(2.*fract(.2497*f*(_t-0.0*(1.+1.3*_sin(.9*_t)))+.02)-1.)+MADD((_t-0.0*(1.+1.3*_sin(.9*_t))),.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD((_t-0.0*(1.+1.3*_sin(.9*_t))),.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.)
      +.4*sinshape(s_atan(2.*(2.*fract(.2497*f*(_t-3.0e-03*(1.+1.3*_sin(.9*_t)))+.02)-1.)+MADD((_t-3.0e-03*(1.+1.3*_sin(.9*_t))),.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD((_t-3.0e-03*(1.+1.3*_sin(.9*_t))),.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.)
      +.4*sinshape(s_atan(2.*(2.*fract(.2497*f*(_t-6.0e-03*(1.+1.3*_sin(.9*_t)))+.02)-1.)+MADD((_t-6.0e-03*(1.+1.3*_sin(.9*_t))),.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD((_t-6.0e-03*(1.+1.3*_sin(.9*_t))),.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.))
      +.9*vel*(.8*(2.*fract(.2497*f*_t+.02)-1.)+MADD(_t,.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD(_t,.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.))*exp(-6.*max(_t-.115,0.))
      +.3*vel*sinshape(s_atan(2.*(2.*fract(.2497*f*_t+.02)-1.)+MADD(_t,.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD(_t,.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.)*exp(-6.*max(_t-.115,0.));
                        amaysynR = .9*vel*(.4*sinshape(s_atan(2.*(2.*fract(.2497*f*(_t2-0.0*(1.+1.3*_sin(.9*_t2)))+.02)-1.)+MADD((_t2-0.0*(1.+1.3*_sin(.9*_t2))),.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD((_t2-0.0*(1.+1.3*_sin(.9*_t2))),.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.)
      +.4*sinshape(s_atan(2.*(2.*fract(.2497*f*(_t2-3.0e-03*(1.+1.3*_sin(.9*_t2)))+.02)-1.)+MADD((_t2-3.0e-03*(1.+1.3*_sin(.9*_t2))),.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD((_t2-3.0e-03*(1.+1.3*_sin(.9*_t2))),.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.)
      +.4*sinshape(s_atan(2.*(2.*fract(.2497*f*(_t2-6.0e-03*(1.+1.3*_sin(.9*_t2)))+.02)-1.)+MADD((_t2-6.0e-03*(1.+1.3*_sin(.9*_t2))),.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD((_t2-6.0e-03*(1.+1.3*_sin(.9*_t2))),.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.))
      +.9*vel*(.8*(2.*fract(.2497*f*_t2+.02)-1.)+MADD(_t2,.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD(_t2,.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.))*exp(-6.*max(_t2-.115,0.))
      +.3*vel*sinshape(s_atan(2.*(2.*fract(.2497*f*_t2+.02)-1.)+MADD(_t2,.5009*f,0.,16,1,-1.,10.*f,1.,10.,1.,.015,1.,0.,0.)+.8*MADD(_t2,.25*f,0.,32,1,-.6,1000.,100.,0.,1.,.02,.3,0.,0.)),6.,3.)*exp(-6.*max(_t2-.115,0.));
                    }
                    else if(syn == 138){
                        
                        amaysynL = maceboss_vol(BT)*(env_AHDSR(Bprog,L,.007,0.,.01,1.,.01)*(sinshape(MADD(_t,.5*f,0.,256,2,1.+.9*(.55+(.4*clip((1.+1.)*_sin(4.*BT)))),(1088.+(863.*_sin_(2.*BT,.4))),10.,5.37,3.85,.005,.4*(.55+(.4*clip((1.+1.)*_sin(4.*BT)))),0.,0.),1.,3.)+.8*clip((1.+.5)*_sin(.499*f*_t))+.4*_sq_(1.01*f*_t,.95)));
                        amaysynR = maceboss_vol(BT)*(env_AHDSR(Bprog,L,.007,0.,.01,1.,.01)*(sinshape(MADD(_t2,.5*f,0.,256,2,1.+.9*(.55+(.4*clip((1.+1.)*_sin(4.*BT)))),(1088.+(863.*_sin_(2.*BT,.4))),10.,5.37,3.85,.005,.4*(.55+(.4*clip((1.+1.)*_sin(4.*BT)))),0.,0.),1.,3.)+.8*clip((1.+.5)*_sin(.499*f*_t2))+.4*_sq_(1.01*f*_t2,.95)));
                    }
                    
                    sL += amtL * trk_norm(trk) * s_moothmin(clamp(env,0.,1.) * amaysynL);
                    sR += amtR * trk_norm(trk) * s_moothmin(clamp(env,0.,1.) * amaysynR);
                }
            }
            
        }
    }
    float masterL = .5 * sidechain * sL + .67 * dL;
    float masterR = .5 * sidechain * sR + .67 * dR;
    return vec2(
        masterL,
        masterR);
}

void main()
{
    Tsample = 1./iSampleRate;
    float t = (iBlockOffset + gl_FragCoord.x + gl_FragCoord.y*iTexSize) * Tsample;
    vec2 s = mainSynth(t);
    vec2 v  = floor((0.5+0.5*s)*65535.0);
    vec2 vl = mod(v,256.0)/255.0;
    vec2 vh = floor(v/256.0)/255.0;
    gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}
