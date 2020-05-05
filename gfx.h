const char *gfx_frag =
""
""
"\n#version 130\n"
""
"out vec4 gl_FragColor;"
""
"uniform float iTime;"
""
"const vec2 L=vec2 (1920,1080);"
"const float i1=9.;"
""
"const float N1=acos(-1.);"
"const vec3 c=vec3 (1.,0.,-1.);"
"const float v=.45,"
"u1=7.,"
"W=-.01;"
""
"void e1(in vec3 r,out float i)"
"{"
"r=fract(r*.1031);"
"r+=dot(r,r.yzx+33.33);"
"i=fract((r.x+r.y)*r.z);"
"}"
""
"void C1(in vec2 z1,out float i)"
"{"
"vec3 r=fract(vec3 (z1.xyx)*.1031);"
"r+=dot(r,r.yzx+33.33);"
"i=fract((r.x+r.y)*r.z);"
"}"
""
"void F(in vec3 o,in vec3 N,out float i)"
"{"
"vec3 b=abs(o)-N;"
"i=length(max(b,0.0))"
"+min(max(b.x,max(b.y,b.z)),0.0);"
"}"
""
"void o1(in vec3 o,in vec3 N,in float D,out float i)"
"{"
"F(o,N,i);"
""
"float b;"
"F(o,N+c.zzx*D,b);"
"i=max(i,-b);"
"F(o,N+c.xzz*D,b);"
"i=max(i,-b);"
"F(o,N+c.zxz*D,b);"
"i=max(i,-b);"
"}"
""
"void g1(in float v1,in float y1,in float b1,out float i)"
"{"
"vec2 M=vec2 (y1,abs(v1)-0.5*b1);"
"i=min(max(M.x,M.y),0.0)+length(max(M,0.0));"
"}"
""
"void P(in vec2 k,in vec2 S,out vec2 m)"
"{"
"m=(k.x<S.x)?k:S;"
"}"
""
"void w1(in vec2 k,in vec2 S,out vec2 m)"
"{"
"m=(k.x>S.x)?abs(k):abs(S)*c.zx;"
"}"
""
"float _1(in float i)"
"{"
"return smoothstep(1.5/L.y,-1.5/L.y,i);"
"}"
""
"void C(in vec3 o,out vec2 m)"
"{"
"float i,"
"f=v,"
"s;"
""
"m=vec2 (o.z+v-W,0.);"
""
"for(float l=0.;l<u1;l+=1.)"
"{"
"vec3 _=mod(o,f*c.xxx)-.5*f,"
"g=(o-_)/f*2.;"
"if(max(abs(g.x),abs(g.y))>pow(2.,l+1.))"
"break;"
"e1(1.e2*g-100.*l,s);"
"if(s>.4)"
"{"
"o1(_,(.5*s-.05*l)*f*c.xxx,(.2-.01*l)*(.5*s-.05*l)*f,i);"
"s=2.*(s-.5);"
"if(s<.5)P(m,vec2 (i,1.),m);"
"else P(m,vec2 (i,-1.),m);"
"}else break;"
"if(s>.6)"
"{"
"F(_,(.45*s-.05*l)*f*c.xxx,i);"
"s=2.5*(s-.6);"
"P(m,vec2 (i,2.+floor(3.*s)),m);"
"}"
"f/=2.;"
"}"
"vec3 g=vec3 (-1.,-1.,-1.);"
"e1(1.e2*g,s);"
"o1(o-.5*g*v,.525*s*v*c.xxx,.28*.5*s*v,i);"
"P(m,vec2 (i,1.),m);"
"}"
""
"\n#define m1(n ,p )void n(in vec3 o ,out vec3 d ,in float X ){vec2 a ,w ;p(o ,a );p(o +X *c .xyy,w );d .x=w .x;p(o +X *c .yxy,w );d .y=w .x;p(o +X *c .yyx,w );d .z=w .x;d =normalize(d -a .x);}\n"
"m1(q,C)"
""
"\n#define Q(f1 ,r1 ,t1 ,step)void f1(out vec3 o ,in vec3 n ,inout float i ,in vec3 x ,in int h ,out int l ,out vec2 a ){for(l =0;l <h ;++l ){o =n +i *x ;r1(o ,a );if(a .x<1.e-4)return ;if(t1 ){l =h ;}i +=step;}}\n"
"Q(V,C,max(max(abs(o.x),abs(o.y)),abs(o.z))>1.1*v,min(a.x,5.e-4))"
"Q(p1,C,o.z>.5,min(a.x,1.e-2))"
""
"void B(in vec3 n,in vec3 x,in vec3 f,out float i)"
"{"
"vec3 u=min((f-n)/x,(-f-n)/x);vec2 s1=abs(n.yz+u.x*x.yz),"
"n1=abs(n.xz+u.y*x.xz),"
"l1=abs(n.xy+u.z*x.xy);"
"vec4 y=100.*c.xyyy;"
""
"y=mix(y,vec4 (u.x,c.xyy),float (all(lessThan(s1,f.yz)))*step(u.x,y.x));"
"y=mix(y,vec4 (u.y,c.yxy),float (all(lessThan(n1,f.xz)))*step(u.y,y.x));"
"y=mix(y,vec4 (u.z,c.yyx),float (all(lessThan(l1,f.xy)))*step(u.z,y.x));"
""
"i=y.r;"
"}"
""
"void G(in vec3 o,in vec3 d,in vec3 x,in vec3 t,inout vec3 e,in vec2 a)"
"{"
"if(a.y==-1.){"
"e=vec3 (0.09,0.09,0.09);"
"e=.6*e"
"+.8*e*max(dot(t-o,d),0.)"
"+1.5*e*pow(max(dot(reflect(t-o,d),x),0.),2.);"
"}"
"else if(a.y==0.){"
"e=.1*c.xxx;"
"e=.1*e"
"+.8*e*max(dot(t-o,d),0.)"
"+.5*e*pow(max(dot(reflect(t-o,d),x),0.),2.);"
"}"
"else if(a.y==1.){"
"e=.4*c.xxx;"
"e=.6*e"
"+.8*e*max(dot(t-o,d),0.)"
"+1.5*e*pow(max(dot(reflect(t-o,d),x),0.),2.);"
"}"
"else if(a.y==2.){"
"e=vec3 (1.00,0.29,0.24);"
"e=.7*e"
"+.8*e*max(dot(t-o,d),0.)"
"+.5*e*pow(max(dot(reflect(t-o,d),x),0.),2.);"
"}"
"else if(a.y==3.){"
"e=vec3 (0.33,0.86,0.10);"
"e=.7*e"
"+.8*e*max(dot(t-o,d),0.)"
"+.5*e*pow(max(dot(reflect(t-o,d),x),0.),2.);"
"}"
"else if(a.y==4.){"
"e=vec3 (1.00,0.59,0.12);"
"e=.7*e"
"+.8*e*max(dot(t-o,d),0.)"
"+.5*e*pow(max(dot(reflect(t-o,d),x),0.),2.);"
"}"
"}"
""
"void x1(out vec4 c1,in vec2 a1)"
"{"
"vec2 J=(a1-.5*L.xy)/L.y,"
"a,h1,"
"O;"
"vec3 e=c.yyy,"
"n=2.*c.zzx-.4*c.yyx+.05*c.zzy,"
"Y=n,"
"s=c.xzy,"
"p=c.yyy,"
"d1=cross(normalize(p-n),-s),"
"x,"
"d,"
"o,"
"z=c.yyy,"
"t,"
"j,"
"R,"
"T,"
"K,"
"E;"
"int h=1500,"
"l;"
"float i=0.,U;"
""
"p+=J.x*s+J.y*d1;"
"x=normalize(p-n);"
""
"B(n,x,v*c.xxx,i);"
"if(i>1.e1)"
"{"
"i=-(n.z+v-W)/x.z;"
"o=n+i*x;"
"C(o,a);"
"}"
"V(o,n,i,x,h,l,a);"
""
"t=c.xzx;"
"j=c.zxx;"
""
"if(l<h)"
"{"
"q(o,d,5.e-5);"
"d=round(d);"
"}"
"else "
"{"
"i=-(n.z+v-W)/x.z;"
"o=n+i*x;"
"C(o,a);"
"q(o,d,5.e-5);"
"d=round(d);"
"}"
""
"G(o,d,x,t,e,a);"
"G(o,d,x,j,z,a);"
"e=mix(e,z,.5);"
""
"if(a.y==0.)"
"{"
"Y=o;"
"R=reflect(x,d);"
"U=.01;"
"B(Y,R,v*c.xxx,U);"
"V(T,Y,U,R,h,l,O);"
""
"if(l<h)"
"{"
"q(T,E,5.e-4);"
"d=round(d);"
"G(T,E,R,t,z,O);"
"G(T,E,R,j,K,O);"
"z=mix(z,K,.5);"
"e=mix(e,z,.1);"
"}"
"}"
""
""
""
""
""
"T=o;"
""
"n=o;"
"x=normalize(t-o);"
"i=1.e-2;"
""
"{"
"float I=1.0;"
"float Z=1.e20;"
"for(int l=0;l<h;++l)"
"{"
"o=n+i*x;"
"C(o,a);"
"if(a.x<1.e-4)"
"{"
"I=0.;"
"break;"
"}"
"if(o.z>.5)break;"
"float _=a.x*a.x/(2.0*Z);"
"float b=sqrt(a.x*a.x-_*_);"
"I=min(I,100.0*b/max(0.0,i-_));"
"Z=a.x;"
"i+=min(a.x,5.e-3);"
"}"
"e=mix(.3*e,e,I);"
"}"
""
""
"e=e+1.*e*e*e;"
""
"c1=vec4 (clamp(e,0.,1.),1.);"
"}"
""
"void main()"
"{"
"vec4 e=vec4 (0.);"
"float A=sqrt(i1)-1.;"
"for(float l=-.5*A;l<=.5*A;l+=1.)"
"for(float H=-.5*A;H<=.5*A;H+=1.)"
"{"
"vec4 z;"
"x1(z,gl_FragCoord.xy+vec2 (l,H)*1./max(A,1.));"
"e+=z;"
"}"
"e/=i1;"
"gl_FragColor=e;"
"}"
""
;
