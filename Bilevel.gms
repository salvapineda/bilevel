*Linear bilevel programming

$OFFSYMLIST OFFSYMXREF OFFUELLIST OFFUELXREF offlisting

$Offdigit

option lp=cplex;
option mip=cplex;
option mpec=nlpec;
option mcp=nlpec
option nlp=conopt;
option minlp=sbb;
option iterlim = 1e9;
option limrow = 0;
option limcol = 0;
option solprint=off;
option sysout=off;
option reslim = 21600;

execseed = 1+gmillisec(jnow);

$include .filedata1.txt

$include .filedata2.txt

alias(m,mm);
alias(q,qq);
alias(r,rr);

SET iter/1*20/;
PARAMETER t, rho, M1(m),M2(m),M3(q),M4(q),M5(r),M6(r);

SETS sos /1,2/;

VARIABLES
OF_UL
OF_LL
;

POSITIVE VARIABLES
x(n)
y(m)
alpha(q)
beta(r)
;

BINARY VARIABLES
u1(m)
u2(q)
u3(r)
;

SOS1 VARIABLES
s1(m,sos)
s2(q,sos)
s3(r,sos);

EQUATIONS
UL_OF,UL_C1,LL_OF,LL_C1,LL_C2,LL_D1,
Com_reg,Com_pen,Com_mip_1,Com_mip_2,Com_mip_3,Com_mip_4,Com_mip_5,Com_mip_6,
SO1_1,SO1_2,SO2_1,SO2_2,SO3_1,SO3_2
;

UL_OF..		OF_UL =e= sum(n,c(n)*x(n))+sum(m,d(m)*y(m));
UL_C1(p)..	sum(n,A(p,n)*x(n)) =l= b(p);

LL_OF..		OF_LL =e= sum(m,e(m)*y(m));
LL_C1(q)..	h(q) - sum(n,F(q,n)*x(n)) - sum(m,G(q,m)*y(m)) =g= 0;
LL_C2(r)..	j(r) - sum(m,I(r,m)*y(m)) =g= 0;
LL_D1(m)..	e(m) + sum(q,G(q,m)*alpha(q)) + sum(r,I(r,m)*beta(r)) =g= 0;

Com_reg..		sum(m,y(m)*(e(m) + sum(q,G(q,m)*alpha(q)) + sum(r,I(r,m)*beta(r)))) + sum(q,alpha(q)*(h(q) - sum(n,F(q,n)*x(n)) - sum(m,G(q,m)*y(m)))) + 
			sum(r,beta(r)*(j(r) - sum(m,I(r,m)*y(m)))) =l= t;
Com_pen..		OF_UL =e= sum(n,c(n)*x(n))+sum(m,d(m)*y(m)) + rho*(sum(m,y(m)*(e(m) + sum(q,G(q,m)*alpha(q)) + sum(r,I(r,m)*beta(r)))) + sum(q,alpha(q)*(h(q) - 
			sum(n,F(q,n)*x(n)) - sum(m,G(q,m)*y(m)))) + sum(r,beta(r)*(j(r) - sum(m,I(r,m)*y(m)))));

Com_mip_1(m)..	y(m) =l= u1(m)*M1(m);
Com_mip_2(m)..	e(m) + sum(q,G(q,m)*alpha(q)) + sum(r,I(r,m)*beta(r)) =l= (1-u1(m))*M2(m);
Com_mip_3(q)..	alpha(q) =l= u2(q)*M3(q);
Com_mip_4(q)..	h(q) - sum(n,F(q,n)*x(n)) - sum(m,G(q,m)*y(m)) =l= (1-u2(q))*M4(q);
Com_mip_5(r)..	beta(r) =l= u3(r)*M5(r);
Com_mip_6(r)..	j(r) - sum(m,I(r,m)*y(m)) =l= (1-u3(r))*M6(r);

SO1_1(m)..      s1(m,'1') =e= e(m) + sum(q,G(q,m)*alpha(q)) + sum(r,I(r,m)*beta(r));
SO1_2(m)..      s1(m,'2') =e= y(m);
SO2_1(q)..      s2(q,'1') =e= h(q) - sum(n,F(q,n)*x(n)) - sum(m,G(q,m)*y(m));
SO2_2(q)..      s2(q,'2') =e= alpha(q);
SO3_1(r)..      s3(r,'1') =e= j(r) - sum(m,I(r,m)*y(m));
SO3_2(r)..      s3(r,'2') =e= beta(r);

MODEL BiLin /UL_OF,UL_C1,LL_C1,LL_C2,LL_D1/;

MODEL BiLin2 /UL_OF,LL_OF,UL_C1,LL_C1,LL_C2,LL_D1/;

MODEL LowLev /LL_OF,LL_C1,LL_C2/;

MODEL BiReg /UL_OF,UL_C1,LL_C1,LL_C2,LL_D1,Com_reg/;

MODEL BiPen /UL_C1,LL_C1,LL_C2,LL_D1,Com_pen/;

MODEL BiMIP /UL_OF,UL_C1,LL_C1,LL_C2,LL_D1,Com_mip_1,Com_mip_2,Com_mip_3,Com_mip_4,Com_mip_5,Com_mip_6/;

MODEL BiSOS /UL_OF,UL_C1,LL_C1,LL_C2,LL_D1,SO1_1,SO1_2,SO2_1,SO2_2,SO3_1,SO3_2/;

$onecho > cplex.op2	
	mipstart 1
	parallelmode 0
	threads 1
	nodesel 2
	varsel 3
	nodefileind 3
	workmem 56
$offecho

$onecho > cplex.opt
	mipstart 1
	parallelmode 0
	threads 1	
$offecho

BiMIP.optcr = 0;
BiMIP.optfile = 1;
BiSOS.optcr = 0;
BiSOS.optfile = 2;

SCALARS modelstat, solvestat, time, optca, optcr;

time = 0;

if (met eq 2,
  SOLVE BiLin minimizing OF_UL using lp;
	time = time + BiLin.resusd;
	x.up(n) = x.l(n);
	x.lo(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	time = time + LowLev.resusd;
	x.up(n) = inf;
	x.lo(n) = 0;	
	s1.l(m,'1') = e(m) + sum(q,G(q,m)*LL_C1.m(q)) + sum(r,I(r,m)*LL_C2.m(r));
  s1.l(m,'2') = y.l(m);
	s2.l(q,'1') = h(q) - sum(n,F(q,n)*x.l(n)) - sum(m,G(q,m)*y.l(m));
	s2.l(q,'2') = LL_C1.m(q);
	s3.l(r,'1') = j(r) - sum(m,I(r,m)*y.l(m));
	s3.l(r,'2') = LL_C2.m(r);		
	SOLVE BiSOS minimizing OF_UL using mip;
	time = BiSOS.resusd;
	modelstat = BiSOS.modelstat;
	solvestat = BiSOS.solvestat;
	x.fx(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	OF_LL.fx = OF_LL.l;
	SOLVE BiLin2 minimizing OF_UL using lp;
elseif (met ge 3) and (met le 13),
  if (met eq 3, M1(m) = 5; M2(m) = 5; M3(q) = 5; M4(q) = 5; M5(r) = 5; M6(r) = 5;
  elseif met eq 4, M1(m) = 5; M2(m) = 5; M3(q) = 5; M4(q) = 5; M5(r) = 5; M6(r) = 5;
  elseif met eq 5, M1(m) = 10; M2(m) = 10; M3(q) = 10; M4(q) = 10; M5(r) = 10; M6(r) = 10;
  elseif met eq 6, M1(m) = 20; M2(m) = 20; M3(q) = 20; M4(q) = 20; M5(r) = 20; M6(r) = 20;
  elseif met eq 7, M1(m) = 50; M2(m) = 50; M3(q) = 50; M4(q) = 50; M5(r) = 50; M6(r) = 50;
  elseif met eq 8, M1(m) = 100; M2(m) = 100; M3(q) = 100; M4(q) = 100; M5(r) = 100; M6(r) = 100;
  elseif met eq 9, M1(m) = 200; M2(m) = 200; M3(q) = 200; M4(q) = 200; M5(r) = 200; M6(r) = 200;
  elseif met eq 10, M1(m) = 500; M2(m) = 500; M3(q) = 500; M4(q) = 500; M5(r) = 500; M6(r) = 500;
  elseif met eq 11, M1(m) = 1000; M2(m) = 1000; M3(q) = 1000; M4(q) = 1000; M5(r) = 1000; M6(r) = 1000;
  elseif met eq 12, M1(m) = 10000; M2(m) = 10000; M3(q) = 10000; M4(q) = 10000; M5(r) = 10000; M6(r) = 10000;
  elseif met eq 13, M1(m) = 100000; M2(m) = 100000; M3(q) = 100000; M4(q) = 100000; M5(r) = 100000; M6(r) = 100000;
  );
	SOLVE BiMIP minimizing OF_UL using mip;
	time = BiMIP.resusd;
	modelstat = BiMIP.modelstat;
	solvestat = BiMIP.solvestat;
	x.fx(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	OF_LL.fx = OF_LL.l;
	SOLVE BiLin2 minimizing OF_UL using lp;
elseif met eq 14,
  SOLVE BiLin minimizing OF_UL using lp;
	time = time + BiLin.resusd;
	x.up(n) = x.l(n);
	x.lo(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	time = time + LowLev.resusd;
	x.up(n) = inf;
	x.lo(n) = 0;		
	LOOP(iter,
		t = 10000*power(10,-(ord(iter)-1));
		SOLVE BiReg minimizing OF_UL using nlp;
		time = time + BiReg.resusd;
	);
	modelstat = BiReg.modelstat;
	solvestat = BiReg.solvestat;	
	x.fx(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	display y.l, LL_C1.m, LL_C2.m;
	OF_LL.fx = OF_LL.l;
	SOLVE BiLin2 minimizing OF_UL using lp;
elseif met eq 15,
  SOLVE BiLin minimizing OF_UL using lp;
	time = time + BiLin.resusd;
	x.up(n) = x.l(n);
	x.lo(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	time = time + LowLev.resusd;
	x.up(n) = inf;
	x.lo(n) = 0;		
	LOOP(iter,
		rho = 1/(power(1.2,-(ord(iter)-1)));
		SOLVE BiPen minimizing OF_UL using nlp;
		time = time + BiPen.resusd;
	);	
	modelstat = BiPen.modelstat;
	solvestat = BiPen.solvestat;
	x.fx(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	OF_LL.fx = OF_LL.l;
	SOLVE BiLin2 minimizing OF_UL using lp;	
elseif (met ge 16) and (met le 18),
  SOLVE BiLin minimizing OF_UL using lp;
	time = time + BiLin.resusd;
	x.up(n) = x.l(n);
	x.lo(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	time = time + LowLev.resusd;
	x.up(n) = inf;
	x.lo(n) = 0;		
	LOOP(iter,
		t = 10000*power(10,-(ord(iter)-1));
		SOLVE BiReg minimizing OF_UL using nlp;
		time = time + BiReg.resusd;
	);
	x.up(n) = x.l(n);
	x.lo(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;
	time = time + LowLev.resusd;
	x.up(n) = inf;
	x.lo(n) = 0;		
	if (met eq 16,
	  M1(m) = 2*smax(mm,y.l(mm));
	  M2(m) = 2*smax(mm,e(mm) + sum(q,G(q,mm)*LL_C1.m(q)) + sum(r,I(r,mm)*LL_C2.m(r)));
	  M3(q) = 2*smax(qq,LL_C1.m(qq));
	  M4(q) = 2*smax(qq,h(qq) - sum(n,F(qq,n)*x.l(n)) - sum(m,G(qq,m)*y.l(m)));
	  M5(r) = 2*smax(rr,LL_C2.m(rr));
	  M6(r) = 2*smax(rr,j(rr) - sum(m,I(rr,m)*y.l(m)));
	elseif met eq 17, 
	  M1(m) = 5*smax(mm,y.l(mm));
	  M2(m) = 5*smax(mm,e(mm) + sum(q,G(q,mm)*LL_C1.m(q)) + sum(r,I(r,mm)*LL_C2.m(r)));
	  M3(q) = 5*smax(qq,LL_C1.m(qq));
	  M4(q) = 5*smax(qq,h(qq) - sum(n,F(qq,n)*x.l(n)) - sum(m,G(qq,m)*y.l(m)));
	  M5(r) = 5*smax(rr,LL_C2.m(rr));
	  M6(r) = 5*smax(rr,j(rr) - sum(m,I(rr,m)*y.l(m)));
	elseif met eq 18, 
	  M1(m) = 10*smax(mm,y.l(mm));
	  M2(m) = 10*smax(mm,e(mm) + sum(q,G(q,mm)*LL_C1.m(q)) + sum(r,I(r,mm)*LL_C2.m(r)));
	  M3(q) = 10*smax(qq,LL_C1.m(qq));
	  M4(q) = 10*smax(qq,h(qq) - sum(n,F(qq,n)*x.l(n)) - sum(m,G(qq,m)*y.l(m)));
	  M5(r) = 10*smax(rr,LL_C2.m(rr));
	  M6(r) = 10*smax(rr,j(rr) - sum(m,I(rr,m)*y.l(m)));  
	);
	u1.l(m)$(y.l(m) gt 1E-6) = 1;
	u1.l(m)$((e(m) + sum(q,G(q,m)*LL_C1.m(q)) + sum(r,I(r,m)*LL_C2.m(r))) gt 1E-6) = 0;
	u2.l(q)$(LL_C1.m(q) gt 1E-6) = 1;
	u2.l(q)$((h(q) - sum(n,F(q,n)*x.l(n)) - sum(m,G(q,m)*y.l(m))) gt 1E-6) = 0;
	u3.l(r)$(LL_C2.m(r) gt 1E-6) = 1;
	u3.l(r)$((j(r) - sum(m,I(r,m)*y.l(m))) gt 1E-6) = 0;
	SOLVE BiMIP minimizing OF_UL using mip;
	time = time + BiMIP.resusd;	
	x.fx(n) = x.l(n);
	SOLVE LowLev minimizing OF_LL using lp;	
	OF_LL.fx = OF_LL.l;
	SOLVE BiLin2 minimizing OF_UL using lp;	
	modelstat = BiMIP.modelstat;
	solvestat = BiMIP.solvestat;
);

$include .filesol2.txt
