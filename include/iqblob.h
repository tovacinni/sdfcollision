// Created by inigo quilez - iq/2019
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//
//
// An animation test - a happy and blobby creature jumping and
// looking around. It gets off-model very often, but it looks
// good enough I think.
//
// Making-of and related math/shader/art explanations (6 hours
// long): https://www[1]outube.com/watch?v=Cfe5UQ-1L9Q
//
// Video capture: https://www[1]outube.com/watch?v=s_UOFo2IULQ
//
// You can buy a metal print of this shader here:
// https://www.redbubble.com/i/metal-print/Happy-Jumping-by-InigoQuilez/43594745.0JXQP

#include <algorithm>
#include <iostream>
#include <Eigen/Core>

//------------------------------------------------------------------
typedef Eigen::Vector2d vec2;
typedef Eigen::Vector3d vec3;
typedef Eigen::Vector4d vec4;
typedef Eigen::Matrix2d mat2x2;

template <typename T>
T clamp(const T& n, const T& lower, const T& upper) {
    return std::max(lower, std::min(n, upper));
}

template <typename T>
T mix(const T& a, const T& b, const T& alpha) {
    return alpha * a + (1.0-alpha) * b;
}

template <typename T>
vec2 mix(const vec2& a, const vec2& b, const T& alpha) {
    return alpha * a + (1.0-alpha) * b;
}


template <typename T>
T smoothstep(const T edge0, const T edge1, const T x) {
    T t;  /* Or genDType t; */
    t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

template <typename T>
T fract(const T& a) {
    return a - floor(a);
}

template <typename T>
T sign(const T& a) {
    return (a > 0.0) ? 1.0 : -1.0;
}

// http://iquilezles.org/www/articles/smin/smin.htm
double smin( double a, double b, double k )
{
    double h = std::max(k-abs(a-b),0.0);
    return std::min(a, b) - h*h*0.25/k;
}

// http://iquilezles.org/www/articles/smin/smin.htm
vec2 smin( vec2 a, vec2 b, double k )
{
    double h = clamp( 0.5+0.5*(b[0]-a[0])/k, 0.0, 1.0 );
    return mix( b, a, h ).array() - k*h*(1.0-h);
}

// http://iquilezles.org/www/articles/smin/smin.htm
double smax( double a, double b, double k )
{
    double h = std::max(k-abs(a-b),0.0);
    return std::max(a, b) + h*h*0.25/k;
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
double sdSphere( vec3 p, double s )
{
    return p.norm()-s;
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
double sdEllipsoid( vec3 p, vec3 r ) // approximated
{
    vec3 r_ = 1.0/r.array();
    vec3 rr_ = 1.0/(r.cwiseProduct(r)).array();
    double k0 = (p.cwiseProduct(r_)).norm();
    double k1 = (p.cwiseProduct(rr_)).norm();
    return k0*(k0-1.0)/k1;
}


vec2 sdStick(vec3 p, vec3 a, vec3 b, double r1, double r2) // approximated
{
    vec3 pa = p-a;
    vec3 ba = b-a;
    double h = clamp( pa.dot(ba)/ba.dot(ba), 0.0, 1.0 );
    vec2 temp;
    temp <<  ( pa - ba*h ).norm() - mix(r1, r2, h*h*(3.0-2.0*h)), h;
    return temp;
}

// http://iquilezles.org/www/articles/smin/smin.htm
vec4 opU( vec4 d1, vec4 d2 )
{
    return (d1[0]<d2[0]) ? d1 : d2;
}

double href;
double hsha;

double sdBlob( vec3 pos, double atime )
{
    hsha = 1.0;
    
    double t1 = fract(atime);
    double t4 = abs(fract(atime*0.5)-0.5)/0.5;

    double p = 4.0*t1*(1.0-t1);
    double pp = 4.0*(1.0-2.0*t1); // derivative of p

    vec3 cen;
    cen << 0.5*(-1.0 + 2.0*t4),
           pow(p,2.0-p) + 0.1,
           floor(atime) + pow(t1,0.7) -1.0;

    // body
    vec2 uu;
    uu << 1.0, -pp;
    uu = uu.normalized();
    vec2 vv;
    //uu << -uu[1], uu[0];
    vv << -uu[1], uu[0];
    
    double sy = 0.5 + 0.5*p;
    double compress = 1.0-smoothstep(0.0,0.4,p);
    sy = sy*(1.0-compress) + compress;
    double sz = 1.0/sy;

    vec3 q = pos - cen;
    double rot = -0.25*(-1.0 + 2.0*t4);
    double rc = cos(rot);
    double rs = sin(rot);
    
    mat2x2 _rcs;
    _rcs << rc, rs, -rs, rc; // If there is a bug, check row major vs col major here
    q.topRows<2>() = _rcs*q.topRows<2>();
    vec3 r = q;
    href = q[1];
    q[1] = uu.dot(q.segment<2>(1));
    q[2] = vv.dot(q.segment<2>(1));
    
    vec4 res;
    res << sdEllipsoid( q, vec3(0.25, 0.25*sy, 0.25*sz) ), 2.0, 0.0, 1.0;
    
    return res[0];
    
    if( res[0]-1.0 < pos[1] ) // bounding volume
    {
    double t2 = fract(atime+0.8);
    double p2 = 0.5-0.5*cos(6.2831*t2);
    r[2] += 0.05-0.2*p2;
    r[1] += 0.2*sy-0.2;
    vec3 sq;
    sq << abs(r[0]), r[1], r[2];

    // head
    vec3 h = r;
    double hr = sin(0.791*atime);
    hr = 0.7*sign(hr)*smoothstep(0.5,0.7,abs(hr));
    
    h[0] = cos(hr) * h[0] - sin(hr) * h[2]; // if bug, swizzle the row / col major
    h[2] = sin(hr) * h[0] + cos(hr) * h[2];
    //h[0]z = mat2x2(cos(hr),sin(hr),-sin(hr),cos(hr))*h[0]z;
    
    vec3 hq;
    hq << abs(h[0]), h[1], h[2];
    vec3 mystery0;
    mystery0 << 0.0,0.20,0.02;
    vec3 mystery1;
    mystery1 << 0.08,0.2,0.15;
    vec3 mystery2;
    mystery2 << 0.0,0.21,-0.1;
    vec3 mystery3;
    mystery3 << 0.20,0.2,0.20;
       double d  = sdEllipsoid( h-mystery0, mystery1 );
    double d2 = sdEllipsoid( h-mystery2, mystery3 );
    d = smin( d, d2, 0.1 );
    res[0] = smin( res[0], d, 0.1 );
    

    // belly wrinkles
    {
    double yy = r[1]-0.02-2.5*r[0]*r[0];
    res[0] += 0.001*sin(yy*120.0)*(1.0-smoothstep(0.0,0.1,abs(yy)));
    }
        
    vec2 _resxz; // rhs resxz
    _resxz << res[0], res[2];

    // arms
    {
    vec3 _temp0;
    _temp0 << 0.18-0.06*hr*sign(r[0]),0.2,-0.05;
    vec3 _temp1;
    _temp1 << 0.3+0.1*p2,-0.2+0.3*p2,-0.15;

    vec2 arms = sdStick( sq, _temp0, _temp1, 0.03, 0.06 );

    vec2 temp = smin( _resxz, arms, 0.01+0.04*(1.0-arms[1])*(1.0-arms[1])*(1.0-arms[1]) );
    res[0] = temp[0];
    res[2] = temp[1];
    }
        
    // ears
    {
    
    double t3 = fract(atime+0.9);
    double p3 = 4.0*t3*(1.0-t3);
    vec3 _temp0;
    _temp0 << 0.15,0.32,-0.05;
    vec3 _temp1;
    _temp1 << 0.2+0.05*p3,0.2+0.2*p3,-0.07;
    vec2 ear = sdStick( hq, _temp0, _temp1, 0.01, 0.04 );
    vec2 temp = smin( _resxz, ear, 0.01 );
    res[0] = temp[0];
    res[2] = temp[1];
    }
    
    // mouth
    {
    vec3 _temp0;
    _temp0 << 0.0,0.15+4.0*hq[0]*hq[0],0.15;
    vec3 _temp1;
    _temp1 << 0.1,0.04,0.2;
    d = sdEllipsoid( h-_temp0, _temp1 );
    res[3] = 0.3+0.7*clamp( d*150.0,0.0,1.0); //res.w
    res[0] = smax( res[0], -d, 0.03 );
    }

    // legs
    {
    double t6 = cos(6.2831*(atime*0.5+0.25));
    double ccc = cos(1.57*t6*sign(r[0]));
    double sss = sin(1.57*t6*sign(r[0]));
    vec3 base;
    base << 0.12,-0.07,-0.1;
    base[1] -= 0.1/sy;
    vec2 legs = sdStick( sq, base, base + vec3(0.2,-ccc,sss)*0.2, 0.04, 0.07 );
    vec2 temp = smin( _resxz, legs, 0.07 );
    res[0] = temp[0];
    res[2] = temp[1];
    }
        
    // eye
    {
    vec3 _temp0 = vec3(0.08,0.27,0.06);
    double blink = pow(0.5+0.5*sin(2.1*atime),20.0);
    double eyeball = sdSphere(hq-_temp0,0.065+0.02*blink);
    res[0] = smin( res[0], eyeball, 0.03 );
    
    
    vec3 cq;
    cq << hq[0] - 0.1,hq[1]-0.34,hq[2]-0.08;
    mat2x2 temp;
    temp << 0.8, 0.6, -0.6, 0.8; // beware rm vs cm
    cq.segment<2>(0) = temp*cq.segment<2>(0);
    vec3 temp2;
    temp2 << 0.06, 0.03, 0.03;
    d = sdEllipsoid( cq, temp2 );
    res[0] = smin( res[0], d, 0.03 );

    vec2 temp3;
    temp3 << 0.095,0.285;
    vec2 temp4;
    temp4 << 1.0,1.1;
    double eonorm = ((hq.segment<2>(0)-temp3).cwiseProduct(temp4)).norm();
    double eo = 1.0-0.5*smoothstep(0.01,0.04,eonorm);
    vec3 temp5;
    temp5 << 0.08,0.28,0.08;
    vec3 temp6;
    temp6 << 0.075,0.28,0.102;
    vec4 temp7;
    temp7 << sdSphere(hq-temp5,0.060),3.0,0.0,eo;
    vec4 temp8;
    temp8 << vec4(sdSphere(hq-temp6,0.0395),4.0,0.0,1.0);

    res = opU( res, temp7);
    res = opU( res, temp8);
    }
    }
    /*
    
    // ground
    double fh = -0.1 - 0.05*(sin(pos[0]*2.0)+sin(pos[2]*2.0));
    double t5f = fract(atime+0.05);
    double t5i = floor(atime+0.05);
    double bt4 = abs(fract(t5i*0.5)-0.5)/0.5;
    vec2  bcen = vec2( 0.5*(-1.0+2.0*bt4),t5i+pow(t5f,0.7)-1.0 );
    
    double k = length(pos[0]z-bcen);
    double tt = t5f*15.0-6.2831 - k*3.0;
    fh -= 0.1*exp(-k*k)*sin(tt)*exp(-max(tt,0.0)/2.0)*smoothstep(0.0,0.01,t5f);
    double d = pos[1] - fh;
    
    // bubbles
    {
    vec3 vp = vec3( mod(abs(pos[0]),3.0)-1.5,pos[1],mod(pos[2]+1.5,3.0)-1.5);
    vec2 id = vec2( floor(pos[0]/3.0), floor((pos[2]+1.5)/3.0) );
    double fid = id[0]*11.1 + id[1]*31.7;
    double fy = fract(fid*1.312+atime*0.1);
    double y = -1.0+4.0*fy;
    vec3  rad = vec3(0.7,1.0+0.5*sin(fid),0.7);
    rad -= 0.1*(sin(pos[0]*3.0)+sin(pos[1]*4.0)+sin(pos[2]*5.0));
    double siz = 4.0*fy*(1.0-fy);
    double d2 = sdEllipsoid( vp-vec3(0.5,y,0.0), siz*rad );
    
    d2 -= 0.03*smoothstep(-1.0,1.0,sin(18.0*pos[0])+sin(18.0*pos[1])+sin(18.0*pos[2]));
    d2 *= 0.6;
    d2 = min(d2,2.0);
    d = smin( d, d2, 0.32 );
    if( d<res[0] ) { res = vec4(d,1.0,0.0,1.0); hsha=sqrt(siz); }
    }

    // candy
    {
    double fs = 5.0;
    vec3 qos = fs*vec3(pos[0], pos[1]-fh, pos[2] );
    vec2 id = vec2( floor(qos[0]+0.5), floor(qos[2]+0.5) );
    vec3 vp = vec3( fract(qos[0]+0.5)-0.5,qos[1],fract(qos[2]+0.5)-0.5);
    vp[0]z += 0.1*cos( id[0]*130.143 + id[1]*120.372 + vec2(0.0,2.0) );
    double den = sin(id[0]*0.1+sin(id[1]*0.091))+sin(id[1]*0.1);
    double fid = id[0]*0.143 + id[1]*0.372;
    double ra = smoothstep(0.0,0.1,den*0.1+fract(fid)-0.95);
    d = sdSphere( vp, 0.35*ra )/fs;
    if( d<res[0] ) res = vec4(d,5.0,qos[1],1.0);
    }
    
    return res;
    */
    return res(0);
}
