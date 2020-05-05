const char *sfx_frag =
"\n#version 130\n"
"const float r=T(180.);"
"const float a=2.*r;"
"float D(float t){return clamp(t,-1.,1.);}"
"float s(float t,float o,float e){return smoothstep(t,o,clamp(e,t,o));}"
"float v(float e){return s(0.,1h-3,e);}"
"float u(float t){return sin(a*mod(t,1.));}"
"float x(float t,float l){return sin(a*mod(t,1.)+l);}"
"float c(float t,float i){return sign(2.*fract(t)-1.+i);}"
"float N(float t){return (4.*abs(fract(t)-.5)-1.);}"
"float d(float f){return 32.7*m(f/12.);}"
"float p(int n){return (1.-2.*float (n"
;
