#include	<md.h>
#include	<system.h>
#include	<math.h>

sysForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
int		i,j;
double		dedr[3],r[3];
double h11ph22, h11mh22, beta, ct, DELTA,dgpdr[3],dgmdr[3], eb1, eb2;
double  gasp(), gasm(), sM(), s1,s2,ds1,ds2, sqD, sqrt(), exp(), getCClKE();
tripdouble	r0,r1,r2,grad[2];
if (nsolute != 3) return;
i = natoms - nsolute;

r0.fx = pos[i+1].fx - pos[i].fx;
r0.fy = pos[i+1].fy - pos[i].fy;
r0.fz = pos[i+1].fz - pos[i].fz;
mvimage(&r0);
r1.fx = pos[i+2].fx - pos[i].fx;
r1.fy = pos[i+2].fy - pos[i].fy;
r1.fz = pos[i+2].fz - pos[i].fz;
mvimage(&r1);
//RC1 is (natoms-3)-(natoms-2) distance
RC1 = r[0] = sqrt(r0.fx*r0.fx + r0.fy*r0.fy + r0.fz*r0.fz);
RC2 = r[1] = sqrt(r1.fx*r1.fx + r1.fy*r1.fy + r1.fz*r1.fz);
r[2] = r0.fx*r1.fx+r0.fy*r1.fy+r0.fz*r1.fz;

cosTheta = ct = r[2]/(r[0]*r[1]);/* cos angle*/
h11ph22 = gasp(r[0],r[1],ct,dgpdr);
h11mh22 = gasm(r[0],r[1],ct,dgmdr);
s1 = sM(r[0],&ds1);/*ds1 = ds(r1)/dr1*/
s2 = sM(r[1],&ds2);/*ds2 = ds(r2)/dr2*/
beta = Qbeta*s1*s2;
DELTA = h11mh22+(VINT_WI1+VINT_WD1)+(VINT_BI1+VINT_BD1)+(VINT_CDI1+VINT_CDD1)- 
                ((VINT_WI2+VINT_WD2)+(VINT_BI2+VINT_BD2)+(VINT_CDI2+VINT_CDD2)); 
s_coord = h11mh22 - DELTA;
sqD = sqrt(DELTA*DELTA+4*beta*beta);
VSYS = (h11ph22-sqD)/2.;
rVIB = 0.0;
EKVIB += getCClKE(&rVIB);
if (beta < 1.0e-6){
	if (DELTA < 0.0) C1sq = 1.0;
	if (DELTA > 0.0) C1sq = 0.0;
}
else{ 
	C1sq = 2*beta*beta/(DELTA*DELTA+4*beta*beta+DELTA*sqrt(DELTA*DELTA+4*beta*beta));
}
/*update forces on solute atoms due to gas phase potential and coupling*/

/* First calculate derivatives with respect to r1, r2 (each diveded by r1 *
 * or r2 respectively) and cos theta (devided by minus r1*r2 )            */

dedr[0] = (0.5*dgpdr[0] - 2*(beta/sqD)*Qbeta*s2*ds1- 0.5*(DELTA/sqD)*dgmdr[0])/r[0];
dedr[1] = (0.5*dgpdr[1] - 2*(beta/sqD)*Qbeta*s1*ds2- 0.5*(DELTA/sqD)*dgmdr[1])/r[1];
dedr[2] = -(0.5*dgpdr[2]- 0.5*(DELTA/sqD)*dgmdr[2])/(r[0]*r[1]);
/*add biasing potential*/
eb1 = Bbias*exp(-beta_bias*fabs(r[0]-r[1]));
eb2 = Abias*exp(-alpha_bias*(r[0]-r[1])*(r[0]-r[1]));
Vbias = eb2-eb1;
dedr[0] += (eb2*alpha_bias*2*(r[0]-r[1])-eb1*beta_bias*sgn(r[0]-r[1]))/r[0];
dedr[1] -= (eb2*alpha_bias*2*(r[0]-r[1])-eb1*beta_bias*sgn(r[0]-r[1]))/r[1];
VSYS -= Vbias;
/*add reaction coordinate windowing constraints*/
VRC = 0.;
if (x1RC > x2RC){
	fprintf(stderr,"sysForce: reaction window not set properly\n");
	exit(1);
}
if (r[0]-r[1] > x2RC){
	VRC = 0.5*kRC*sq(r[0]-r[1]-x2RC)*(r[0]-r[1]-x2RC);
	dedr[0] += 1.5*kRC*sq(r[0]-r[1]-x2RC)/r[0];
	dedr[1] -= 1.5*kRC*sq(r[0]-r[1]-x2RC)/r[1];
}
if (r[0]-r[1] < x1RC){
	VRC = -0.5*kRC*sq(r[0]-r[1]-x1RC)*(r[0]-r[1]-x1RC);
	dedr[0] -= 1.5*kRC*sq(r[0]-r[1]-x1RC)/r[0];
	dedr[1] += 1.5*kRC*sq(r[0]-r[1]-x1RC)/r[1];
}
VSYS += VRC;
/*update forces*/
force[i+1].fx -= dedr[0]*r0.fx;
force[i+1].fy -= dedr[0]*r0.fy;
force[i+1].fz -= dedr[0]*r0.fz;
force[i+2].fx -= dedr[1]*r1.fx;
force[i+2].fy -= dedr[1]*r1.fy;
force[i+2].fz -= dedr[1]*r1.fz;
force[i].fx += dedr[0]*r0.fx + dedr[1]*r1.fx;
force[i].fy += dedr[0]*r0.fy + dedr[1]*r1.fy;
force[i].fz += dedr[0]*r0.fz + dedr[1]*r1.fz;
force[i+1].fx += (grad[0].fx = dedr[2]*(r1.fx-r[2]*r0.fx/(r[0]*r[0])));
force[i+1].fy += (grad[0].fy = dedr[2]*(r1.fy-r[2]*r0.fy/(r[0]*r[0])));
force[i+1].fz += (grad[0].fz = dedr[2]*(r1.fz-r[2]*r0.fz/(r[0]*r[0])));
force[i+2].fx += (grad[1].fx = dedr[2]*(r0.fx-r[2]*r1.fx/(r[1]*r[1])));
force[i+2].fy += (grad[1].fy = dedr[2]*(r0.fy-r[2]*r1.fy/(r[1]*r[1])));
force[i+2].fz += (grad[1].fz = dedr[2]*(r0.fz-r[2]*r1.fz/(r[1]*r[1])));
force[i].fx -= grad[1].fx + grad[0].fx;
force[i].fy -= grad[1].fy + grad[0].fy;
force[i].fz -= grad[1].fz + grad[0].fz;

/*update forces on solute atoms due to EVB interactions with solvents*/
force[i].fx -= (DELTA/sqD)*fevb[i].fx;
force[i].fy -= (DELTA/sqD)*fevb[i].fy;
force[i].fz -= (DELTA/sqD)*fevb[i].fz;
force[i+1].fx -= (DELTA/sqD)*fevb[i+1].fx;
force[i+1].fy -= (DELTA/sqD)*fevb[i+1].fy;
force[i+1].fz -= (DELTA/sqD)*fevb[i+1].fz;
force[i+2].fx -= (DELTA/sqD)*fevb[i+2].fx;
force[i+2].fy -= (DELTA/sqD)*fevb[i+2].fy;
force[i+2].fz -= (DELTA/sqD)*fevb[i+2].fz;
/*update forces on solvent atoms due to EVB interactions with solute*/
for (j=0;j<i;j++){
	force[j].fx -= (DELTA/sqD)*fevb[j].fx;
	force[j].fy -= (DELTA/sqD)*fevb[j].fy;
	force[j].fz -= (DELTA/sqD)*fevb[j].fz;
}

/*projection of total force along reaction coordinate*/
if (r[0]-r[1] > x1RC && r[0]-r[1] < x2RC){
NW += 1.0;
FRC += ((force[i+1].fx-force[i].fx)*r0.fx + (force[i+1].fy-force[i].fy)*r0.fy +
       (force[i+1].fz-force[i].fz)*r0.fz)/r[0] - ((force[i+2].fx-force[i].fx)*r1.fx +
       (force[i+2].fy-force[i].fy)*r1.fy + (force[i+2].fz-force[i].fz)*r1.fz)/r[1];
}
/*reactive flux calculations*/
if (fluxQ){
	if (etime == 0){
		flux_t = 0;
		init_RC = r[0]-r[1];
	}
	if (etime > 0){
		flux_t++;
		if (flux_t == 1) init_dir = r[0]-r[1]-init_RC;
		if ((r[0]-r[1] - init_RC)*init_dir > 0)
				corflux[flux_t] += 1;
	}
	totc1sq[flux_t] += C1sq;
	cprod[flux_t] = C1sq; 
}
}

/*
 *	**** function to calculate EVB potential -forces / r ****
 *	**** with respect to the internal coordinates r.      ****
 */

double gasp(r1,r2,c,dgpdr) // H11º+H22º
double r1, r2,c,*dgpdr;
{
double vm1, vm2, vion1, vion2, vid, vbend, dtheta, vbend_dt, vid_dr;;
double z, e, ex1, ex2, exp(), pow(), acos();
if (c+1 < 1.0e-8) c = -1;
dtheta = acos(c)-PI;
vbend = Qd*exp(-Qa*(r1+r2))*dtheta*dtheta;
z = sigid/(r1+r2);
ex1 = exp(-amors*(r1-reqmors));
ex2 = exp(-amors*(r2-reqmors));
vm1 = Dmors*ex1*(ex1-2);
vm2 = Dmors*ex2*(ex2-2);
vid =  exp(-Qb*(c+1))*bid*epsid*(pow(z,nid)-pow(z,2));
vion1 = exp(2*Qc*(c+1))*Zion*exp(-2*amors*(r1-rstar))-EA;
vion2 = exp(2*Qc*(c+1))*Zion*exp(-2*amors*(r2-rstar))-EA;
if(RC1<RC2)
   VVIB+=vm1;
else
   VVIB+=vm2;
e = vm1+vm2+vion1+vion2+2*vid+2*EA+2*Dmors+2*vbend;
/*derivative of gas wrs to r1 and r2*/
vid_dr =  - exp(-Qb*(c+1))*bid*epsid*(nid*pow(z,nid+1)-2*pow(z,3))/sigid; 
dgpdr[0] = 2*vid_dr-2*amors*(vion1+EA)-2*amors*Dmors*ex1*(ex1-1)-2*Qa*vbend;
dgpdr[1] = 2*vid_dr-2*amors*(vion2+EA)-2*amors*Dmors*ex2*(ex2-1)-2*Qa*vbend;
/*derivative of gas wrs to cos theta*/
if (c+1 < 1.0e-8)
	vbend_dt = -2*Qd*exp(-Qa*(r1+r2));
else
	vbend_dt = 2*Qd*exp(-Qa*(r1+r2))*dtheta/sqrt(1.-c*c);

dgpdr[2] = 2*Qc*(vion1+EA)+2*Qc*(vion2+EA) - 2*Qb*vid - 2*vbend_dt;
return(e);
}

double gasm(r1,r2,c,dgmdr) // H11º-H22º
double r1, r2,c,*dgmdr;
{
double vm1, vm2, vion1, vion2;
double z, e, ex1, ex2, exp();
ex1 = exp(-amors*(r1-reqmors));
ex2 = exp(-amors*(r2-reqmors));
vm1 = Dmors*ex1*(ex1-2);
vm2 = Dmors*ex2*(ex2-2);
vion1 = exp(2*Qc*(c+1))*Zion*exp(-2*amors*(r1-rstar))-EA;
vion2 = exp(2*Qc*(c+1))*Zion*exp(-2*amors*(r2-rstar))-EA;
/*e = vm1-vm2+vion1-vion2;*/
e = -vm1+vm2+vion1-vion2;
/*derivative of gas wrs to r1 and r2*/
dgmdr[0] = -2*amors*(vion1+EA)+2*amors*Dmors*ex1*(ex1-1);
dgmdr[1] = 2*amors*(vion2+EA)-2*amors*Dmors*ex2*(ex2-1);
/*derivative of gas wrs to cos theta*/

dgmdr[2] = 2*Qc*(vion1+EA)-2*Qc*(vion2+EA);
return(e);
}
double sM(r,ds) // overlap intergal S(r)
double r, *ds;
{
double p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,e1,e2,e3;
double overlap, A[6],B[6], AP[6], BP[6], p,t, exp(), pow();
p = 0.5*(mua+mub)*r/BOHR;
t = (mua-mub)/(mua+mub);
e1 = exp(-p);
e2 = -exp(-p*t);
e3 = -exp(p*t);
p1 = 1/p; p2 = p1*p1; p3 = p1*p2; p4 = p2*p2; p5 = p3*p2; p6 = p3*p3; p7=p6*p1; 
q1 = 1/(p*t); q2 = q1*q1; q3 = q1*q2; q4 = q2*q2; q5 = q3*q2; q6 = q3*q3; q7=q6*q1;
A[0] = e1*p1;
AP[0] = -A[0]-e1*p2;
A[1] = e1*(p1+p2);
AP[1] = -A[1]-e1*(p2+2*p3);
A[2] = e1*(p1+2*p2+2*p3);
AP[2] =  -A[2]-e1*(p2+4*p3+6*p4);
A[3] = e1*(p1+3*p2+6*p3+6*p4);
AP[3] =  -A[3]-e1*(p2+6*p3+18*p4+24*p5);
A[4] = e1*(p1+4*p2+12*p3+24*p4+24*p5);
AP[4] =  -A[4]-e1*(p2+8*p3+36*p4+96*p5+120*p6);
A[5] = e1*(p1+5*p2+20*p3+60*p4+120*p5+120*p6);
AP[5] =  -A[5]-e1*(p2+10*p3+60*p4+240*p5+600*p6+720*p7);
B[0] = e2*q1;
BP[0] = -B[0]-e2*q2;
B[1] = e2*(q1+q2);
BP[1] = -B[1]-e2*(q2+2*q3);
B[2] = e2*(q1+2*q2+2*q3);
BP[2] = -B[2]- e2*(q2+4*q3+6*q4);
B[3] = e2*(q1+3*q2+6*q3+6*q4);
BP[3] = -B[3]-e2*(q2+6*q3+18*q4+24*q5);
B[4] = e2*(q1+4*q2+12*q3+24*q4+24*q5);
BP[4] = -B[4]- e2*(q2+8*q3+36*q4+96*q5+120*q6);
B[5] = e2*(q1+5*q2+20*q3+60*q4+120*q5+120*q6);
BP[5] = -B[5]- e2*(q2+10*q3+60*q4+240*q5+600*q6+720*q7);
B[0] -= e3*q1;
BP[0] -=  e3*q1-e3*q2;
B[1] += e3*(q1-q2);
BP[1] += e3*(q1-q2)-e3*(q2-2*q3);
B[2] += e3*(-q1+2*q2-2*q3);
BP[2] += e3*(-q1+2*q2-2*q3)-
         e3*(-q2+4*q3-6*q4);
B[3] += e3*(q1-3*q2+6*q3-6*q4);
BP[3] += e3*(q1-3*q2+6*q3-6*q4)-
        e3*(q2-6*q3+18*q4-24*q5);
B[4] += e3*(-q1+4*q2-12*q3+24*q4-24*q5);
BP[4] += e3*(-q1+4*q2-12*q3+24*q4-24*q5)-
        e3*(-q2+8*q3-36*q4+96*q5-120*q6);
B[5] += e3*(q1-5*q2+20*q3-60*q4+120*q5-120*q6);
BP[5] += e3*(q1-5*q2+20*q3-60*q4+120*q5-120*q6)-
        e3*(q2-10*q3+60*q4-240*q5+600*q6-720*q7);
*ds = overlap = (1.0/16.0)*1.0/(sqrt(30.0))*(1.0/p6)*pow(1+t,2.5)*pow(1-t,3.5); 
overlap *= (A[2]*(B[1]+B[5]) - A[3]*(B[0]+B[4]) - B[3]*(A[0]+A[4]) + B[2]*(A[1]+A[5]));
*ds *=  (AP[2]*(B[1]+B[5]) - AP[3]*(B[0]+B[4]) - B[3]*(AP[0]+AP[4]) + B[2]*(AP[1]+AP[5])+ t*(A[2]*(BP[1]+BP[5]) - A[3]*(BP[0]+BP[4]) - BP[3]*(A[0]+A[4]) + BP[2]*(A[1]+A[5])));
*ds += overlap*6/p; 
*ds *= (0.5*(mua+mub)/BOHR);
return (overlap);
}

double getCClKE(double *q)
{
tripd d[3],normV;
double radV,sumEK;
int iA,iB;

iA = natoms-3;
if(RC1 < RC2)
   iB = iA+1;
else 
   iB = iA+2;
d[0].fx = pos[iA].fx -pos[iB].fx;
d[0].fy = pos[iA].fy -pos[iB].fy;
d[0].fz = pos[iA].fz -pos[iB].fz;
q[0] = sqrt(d[0].fx*d[0].fx+d[0].fy*d[0].fy+d[0].fz*d[0].fz);

radV = (d[0].fx * (tvel[iA].fx -tvel[iB].fx)+
        d[0].fy * (tvel[iA].fy -tvel[iB].fy)+
        d[0].fz * (tvel[iA].fz -tvel[iB].fz))/q[0];
// vibrational
return(0.5*mass[iA]*mass[iB]/(mass[iA]+mass[iB])*radV*radV);
/*CClKE[0] = 0.5*mass[iA]*mass[iB]/(mass[iA]+mass[iB])*radV*radV;
//translational

CClKE[2] = (sq(mass[iA] * tvel[iA].fx + mass[iB] * tvel[iB].fx)+
       sq(mass[iA] * tvel[iA].fy + mass[iB] * tvel[iB].fy)+
       sq(mass[iA] * tvel[iA].fz + mass[iB] * tvel[iB].fz))/
       (2*(mass[iA]+mass[iB]));
normV.fx = (tvel[iA].fx -tvel[iB].fx)-d[0].fx*radV/q[0];
normV.fy = (tvel[iA].fy -tvel[iB].fy)-d[0].fy*radV/q[0];
normV.fz = (tvel[iA].fz -tvel[iB].fz)-d[0].fz*radV/q[0];
//rotational
CClKE[1] = 0.5*mass[iA]*mass[iB]/(mass[iA]+mass[iB])*
        (normV.fx*normV.fx+normV.fy*normV.fy+normV.fz*normV.fz);

*/
}
