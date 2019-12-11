%Jake Westling
%PNGE420 Project 3
%Fall 2015 Semester
%Clear all old data
close all
clear
clc
%Ask for inputs
fprintf('Inputs:   \n')
qo = input('Enter the oil flow rate in STB/D:  ');
qg = input('Enter the gas flow rate in SCF/D:  ');
d = input('Enter the diameter in inches:  ');
z = input('Enter the distance in miles:  ');
p1 = input('Enter the inlet pressure in psia:  ');
theta = input('Enter the angle from horizontal in degrees:  ');
Tbar = input('Enter the average temperature in degrees F:  ');
oilgrav = input('Enter the oil gravity in degrees API:  ');
gamg = input('Enter the gas specific gravity:  ');
GOR = input('Enter the gas-oil ratio in SCF/STB:  ');
stream = input('Enter 1 for downstream or 2 for upstream:  ');
fprintf('\nOutputs:   \n')
deltap = 100;
fdeltap = 0;
zzz=0;
while abs(deltap-fdeltap)>1
    zzz = zzz+1;
    if zzz > 1
    deltap = fdeltap;
    end
if stream == 1
    pbar = p1+deltap/2;
elseif stream == 2
    pbar = p1-deltap/2;
end

%PVT Equations (given)
%Solving Rs and Pb
gamo = 141.5/(131.5 + oilgrav);
syms Pb Rs 
[solPb,solRs] = solve((1.09373*10^-3)*Rs^.5502*gamg^-1.71956*gamo^2.5486*(Tbar+460)^2.0967 == Pb, (0.0006*Pb^0.856)*(gamg^0.351)*(Tbar^1.829)*(oilgrav^1.462)*(GOR^-2.116)*(pbar^(3.867*Pb^-.306*gamg^-.083*Tbar^-.306*oilgrav^-.288*GOR^.525)) == Rs);
Pb=double(solPb);
Rs=double(solRs);
w = warning('query','last');
id = w.identifier;
warning('off','last');
Rs = Rs +8;
%Formation volume factor oil
F = Rs*(gamg/gamo)^0.5+1.25*(Tbar);
Bo = 0.972+0.000147*F^1.175;
%Oil viscosity
muo = exp(exp(1.8653-0.025086*oilgrav-0.5644*log10(Tbar)))-1+1.5;
%Z factor
Zg = 1.39*(((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))-0.92)^0.5-0.36*((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))-0.101+(1-(1.39*(((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))-0.92)^0.5-0.36*((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))-0.101))/exp(((0.62-0.23*((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0)))*(pbar/(678-50*(gamg-0.5)-206.7*0+440*0+606.7*0))+(0.066/(((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))-0.86)-0.037)*(pbar/(678-50*(gamg-0.5)-206.7*0+440*0+606.7*0))^2+0.32/10^(9*(((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))-1))*(pbar/(678-50*(gamg-0.5)-206.7*0+440*0+606.7*0))^6))+(0.132-0.32*log((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0)))*(pbar/(678-50*(gamg-0.5)-206.7*0+440*0+606.7*0))^(10^(0.3106-0.49*((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))+0.1824*((Tbar+460)/(326+315.7*(gamg-0.5)-240*0-83.3*0+133.3*0))^2));
%Gas viscosity
B5 = 0;
B6 = 0;
B7 = 0;
mug =(((1.709/100000-2.062/1000000*gamg)*Tbar+8.188/1000-6.15/1000*log10(gamg))+(B5*(8.48/1000*log10(gamg)+9.59/1000))+(B6*(9.08/1000*log10(gamg)+6.24/1000))+(B7*(8.49/1000*log10(gamg)+3.73/1000)))/((Tbar+460)/(326+315.7*(gamg-0.5)-240*B5-83.3*B6+133.3*B7))*exp((-2.462)+(2.97)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))+(-2.862/10)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^2+(8.054/1000)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^3+((Tbar+460)/(326+315.7*(gamg-0.5)-240*B5-83.3*B6+133.3*B7))*((2.808)+(-3.498)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))+(3.603/10)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^2+(-1.044/100)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^3)+((Tbar+460)/(326+315.7*(gamg-0.5)-240*B5-83.3*B6+133.3*B7))^2*((-7.933/10)+(1.396)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))+(-1.491/10)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^2+(4.41/1000)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^3)+((Tbar+460)/(326+315.7*(gamg-0.5)-240*B5-83.3*B6+133.3*B7))^3*((8.393/100)+(-1.864/10)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))+(2.033/100)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^2+(-6.095/10000)*(pbar/(678-50*(gamg-0.5)-206.7*B5+440*B6+606.7*B7))^3));
%other stuff
rhog = (0.0764*gamg*pbar*520)/(14.7*(Tbar+460)*Zg);
sigo = 22.3;
gamo = 141.5/(131.5 + oilgrav);
rhoo = (350*gamo + 0.0764*Rs*gamg)/(5.615*Bo);
qgs = 3.27*10^(-7)*Zg*qo*(GOR-Rs)*(Tbar+460)/pbar;
qLs = 6.49e-5*(qo*Bo);
vsL = qLs*4/(pi*(d/12)^2);
vsg = qgs*4/(pi*(d/12)^2);
vm = vsL + vsg;
GL = rhoo*vsL;
Gg = rhog*vsg;
Gm = GL + Gg;
lamda = qLs/(qLs + qgs);
NFR = vm^2/(32.17*d/12);
mum = 6.72e-4*(muo*lamda+mug*(1-lamda));
sigL = sigo;
NRens = Gm*d/12/mum;
NLv = 1.938*vsL*(rhoo/sigL)^0.25;
L1 = 316*lamda^0.302;
L2 = 0.0009252*lamda^(-2.4684);
L3 = 0.10*lamda^(-1.4516);
L4 = 0.5*lamda^(-6.738);
%Flow Type
if (lamda < 0.01 && NFR < L1)||(lamda >= 0.01 && NFR < L2)
    x = 1;
    if abs(deltap-fdeltap)<=1
    fprintf('Flow pattern is Segregated. \n')
    end
    a=0.98;
    b=0.4846;
    c=0.0868;
    dx=0.011;
    ex=-3.768;
    f=3.539;
    g=-1.614;
    HLO = a*lamda^b/NFR^c;
    C = (1-lamda)*log(dx*lamda^ex*NLv^f*NFR^g);
elseif ((lamda < 0.4)&&(lamda >=0.01) && (L3 < NFR)&&(NFR <= L1))||(lamda >= 0.4 && (L3 < NFR)&&(NFR <= L4))
    x = 2;
    if abs(deltap-fdeltap)<=1
    fprintf('Flow pattern is Intermittent. \n')
    end
    a=0.845;
    b=0.5351;
    c=0.0173;
    dx=2.96;
    ex=0.305;
    f=-0.4473;
    g=0.0978;
    HLO = a*lamda^b/NFR^c;
    C = (1-lamda)*log(dx*lamda^ex*NLv^f*NFR^g);
elseif (lamda >= 0.01)&&(L2 < NFR)&&(NFR < L3)
    x = 3;
    if abs(deltap-fdeltap)<=1
    fprintf('Flow pattern is Transition. \n')
    end
    A = (L3-NFR)/(L3-L2);
    B=(1-A);
    a=0.98;
    b=0.4846;
    c=0.0868;
    dx=0.011;
    ex=-3.768;
    f=3.539;
    g=-1.614;
    HLOseg = a*lamda^b/NFR^c;
    a=0.845;
    b=0.5351;
    c=0.0173;
    dx=2.96;
    ex=0.305;
    f=-0.4473;
    g=0.0978;
    HLOint = a*lamda^b/NFR^c;
    HLO = A*HLOseg+B*HLOint;
    a=0.98;
    b=0.4846;
    c=0.0868;
    dx=0.011;
    ex=-3.768;
    f=3.539;
    g=-1.614;
    C = (1-lamda)*log(dx*lamda^ex*NLv^f*NFR^g);
elseif (lamda < 0.4 && NFR >= L1)||(lamda >= 0.4 && NFR < L2)
    x = 4;
    if abs(deltap-fdeltap)<=1
    fprintf('Flow pattern is Distributed. \n')
    end
    a=1.065;
    b=0.5824;
    c=0.0609;
    HLO = a*lamda^b/NFR^c;
    C=0;
end
tri = 1 + C*(sind(1.8*theta)-0.333*(sind(1.8*theta)^3));
HL3 = HLO*tri;
rhotp = rhoo*HL3 + rhog*(1-HL3);
y = lamda/HL3^2;
if 1<y && y<1.2
    S=log(2.2*y-1.2);
else
S = log(y)/(-0.0523+3.182*log(y)-0.8725*(log(y))^2+0.01853*(log(y))^4);
end
ftpoverfns = exp(0.2487);
fns = 1/(2*log10(NRens/(4.5223*log10(NRens)-3.8215)))^2;
ftp = ftpoverfns*fns;
fdeltap = z*5280*(rhotp*sind(theta)+ftp*Gm*vm/(2*32.17*d/12))/(1-(rhotp*vm*vsg)/(32.17*pbar*144))/144;
end
poutlet = p1 - fdeltap;

fprintf('The average pressure is %0.0f psia. \n', pbar)
fprintf('Rs is %0.0f scf/stk bbl.\nBo is %0.3f RB/STB.\nOil viscosity is %0.1f cp.\nGas viscosity is %0.4f cp.\nOil surface tension is %0.1f dynes/cm. \nZ factor is %0.2f. \n', Rs, Bo, muo, mug, sigo, Zg)
fprintf('Specific oil gravity is %0.2f.\nOil density is %0.2f lbm/cu ft.\nGas density is %0.2f lbm/cu ft.\nIn Situ gas rate is %0.3f cu ft/sec.\nIn Situ oil rate is %0.3f cu ft/sec.\nIn Situ oil velocity is %0.3f ft/sec.\nIn Situ gas velocity is %0.3f ft/sec.\nIn Situ mixture velocity is %0.3f ft/sec.\nLiquid mass flux rate is %0.2f lbm/sec-sq ft.\nGas mass flux rate is %0.3f lbm/sec-sq ft.\nTotal mass flux rate is %0.2f lbm/sec-sq ft.\n', gamo, rhoo, rhog, qgs, qLs, vsL, vsg, vm, GL, Gg, Gm)
fprintf('No-slip holdup is %0.3f.\nFroude number is %0.4f.\nMixture viscosity is %0.5f lbm/ft sec.\nLiquid surface tension is %0.1f cu ft/sec.\nNo slip Reynold''s number is %0.0f.\nLiquid viscosity number is %0.3f.\nL1 is %0.0f.\nL2 is %0.4f.\nL3 is %0.3f.\nL4 is %0.0f.\n', lamda, NFR, mum, sigL, NRens, NLv, L1, L2, L3, L4)
fprintf('The horizontal holdup is %0.3f.\nThe inclination correction factor is %0.2f.\n', HLO, C)
fprintf('Liquid holdup inclination correction factor is %0.1f.\nLiquid holdup is %0.3f.\nDenity of two phase is %0.2f lbm/cu ft.\ny is %0.3f.\nS is %0.4f.\nFriction factor ratio is %0.3f.\n', tri, HL3, rhotp, y, S, ftpoverfns)
fprintf('No-slip friction factor is %0.3f.\nTwo-phase friction factor is %0.3f.\nChange in pressure is %0.0f psia.\nOutlet pressure is %0.0f psia.\n', fns, ftp, fdeltap, poutlet)