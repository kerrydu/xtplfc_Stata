* Simulation codes for xtplfc
* 3June2020
* To run the simulation codes, two user-written packages outtable and graph2tex should be installed in advance.
* Christopher F Baum & Joao Pedro Azevedo, 2001. "OUTTABLE: Stata module to write matrix to LaTeX table," Statistical Software Components S419501.
* UCLA Statistical Consulting Group. graph2tex, downloadable from https://stats.idre.ucla.edu/stat/stata/ado/analysis

capture log close
sjlog using DGP1, replace

clear all
set matsize 2000
set obs 50

set seed 123456
matrix C = (1,0,0.42\0,1,0.85\0.42,0.85,1)
drawnorm x2f x3f d, corr(C)
gen id = _n
expand 40
bysort id : gen t = _n
xtset id t
gen y = 0
matrix D = (1,0.2,0.8\0.2,1,0\0.8,0,1)
drawnorm x1 x2e x3e, corr(D)
gen x2 = (x2f+x2e)/sqrt(2)
gen x3 = (x3f+x3e)/sqrt(2)
gen z=rnormal()
gen gf=1*x3+ 2*x3^2 - 0.25*(x3)^3

forv j=1/500{
	cap drop e
	cap drop *_sd
	qui drawnorm e
	qui replace y = x1 -x2 + z*gf+ d + e
	qui xtplfc y x1 x2, z(z) u(x3) maxnk(20) gen(fit`j') fast brep(500)
	matrix B=e(b)
	matrix B=B[1,1..2]
	matrix B1=(nullmat(B1)\B)
	matrix V=e(Vs)
	matrix V=vecdiag(V)
	matrix V=V[1,1..2]
	matrix V1=(nullmat(V1)\V)

	qui xtreg y x1 x2 x3 z,fe
	matrix B=e(b)
	matrix B=B[1,1..2]
	matrix B2=(nullmat(B2)\B)
	matrix V=e(V)
	matrix V=vecdiag(V)
	matrix V=V[1,1..2]
	matrix V2=(nullmat(V2)\V)
	
	
	qui xtreg y x1 x2 x3 z c.x3#c.z,fe
	matrix B=e(b)
	matrix B=B[1,1..2]
	matrix B3=(nullmat(B3)\B)
	matrix V=e(V)
	matrix V=vecdiag(V)
	matrix V=V[1,1..2]
	matrix V3=(nullmat(V3)\V)	
	

}


* Fig.1

egen av_fit = rowmean(fit*)

egen sd_fit = rowsd(fit*)


gen  c =  invnormal(1 - (100 - 95) / 200)

gen low = av_fit - c * sd_fit

gen up = av_fit + c * sd_fit

twoway (rarea low up x3, sort(x3) color(gs7)) ///
   (line av_fit x3, sort(x3) color(black) lpattern(solid)) ///
   (line gf x3, color(gs10) sort lpattern(longdash)), ///
   legend(label(1 confidence interval at 95%) label(3 DGP) label(2 average fit) ///
   cols(3) order(3 1 2)) xtitle("X3",height(5)) ytitle("g(X3)",height(5)) sch(sj) 
graph2tex , epsfile(fig1) caption(Average fit of g(X3)) label(fig1) 


* Table 1	
clear	
set obs 500
mat res1=J(3,8,.)

forv k=1/3{
	qui svmat B`k'
	su B`k'1,meanonly
	mat res1[`k',1]=r(mean)-1
	su B`k'2,meanonly
	mat res1[`k',2]=r(mean)+1
	
	qui svmat V`k'
	qui replace V`k'1=sqrt(V`k'1)
	qui gen B`k'1_lb=B`k'1-invnormal(0.975)*V`k'1
	qui gen B`k'1_ub=B`k'1+invnormal(0.975)*V`k'1
	
	qui gen CIlen`k'1=B`k'1_ub-B`k'1_lb
	qui su CIlen`k'1,d
	mat res1[`k',5]=r(p50)	
	
	qui replace V`k'2=sqrt(V`k'2)
	qui gen B`k'2_lb=B`k'2-invnormal(0.975)*V`k'2
	qui gen B`k'2_ub=B`k'2+invnormal(0.975)*V`k'2
	
	qui gen CIlen`k'2=B`k'2_ub-B`k'2_lb
	qui su CIlen`k'2,d
	mat res1[`k',6]=r(p50)	
	
	qui gen cov`k'1=(B`k'1_lb<=1 & B`k'1_ub>=1) 
	su cov`k'1,meanonly
	mat res1[`k',7]=r(mean)	
	
	qui gen cov`k'2=(B`k'2_lb<=-1 & B`k'2_ub>=-1) 
	su cov`k'2,meanonly
	mat res1[`k',8]=r(mean)			

	qui replace B`k'1=(B`k'1-1)^2
	su B`k'1,meanonly
	mat res1[`k',3]=r(mean)
	qui replace B`k'2=(B`k'2+1)^2
	su B`k'2,meanonly
	mat res1[`k',4]=r(mean)
	
	
}

mat list res1
outtable using res1, mat(res1)  format(%9.4f) replace

sjlog close, replace



///////////////////////////////////////////////////////////////////////////

capture log close
sjlog using DGP2, replace


clear all
set matsize 2000
set obs 50
set seed 789
gen a=rnormal()
gen id=_n
expand 80
bys id: gen year=_n
gen x=10*runiform()+0.5*a
gen u=-9+18*runiform()
gen gf=sin(_pi/3*u)
gen y=0
xtset id year
mata: gfmat=J(2000,500,.)
forv j=1/500{
    
    preserve
    qui replace y=a+0.3*x+0.3*L2.y+L1.y*gf+sqrt(2)*rnormal() if year>2
    qui drop if year<41
    qui gen L_y=L1.y
    qui gen L2_y=L2.y
    qui ivxtplfc y L2_y x, zvars(L_y) uvar(u) gen(g) endoz(1) ///
	                       maxnknots(20) ivz(L2_y,uflag(1)) fast brep(500)
    qui putmata gfhat=g_1,replace
    mata: gfmat[.,`j']=gfhat
    matrix B=e(b)
    matrix B=B[1,1..2]
    matrix B1=(nullmat(B1)\B)
	matrix V=e(Vs)
	matrix V=vecdiag(V)
	matrix V=V[1,1..2]
	matrix V1=(nullmat(V1)\V)	
	

    qui ivregress 2sls D.y D.L2.y D.x (D.L.y=L2.y), noconstant
    matrix B=e(b)
    matrix B=B[1,2..3]
    matrix B2=(nullmat(B2)\B)
	matrix V=e(V)
	matrix V=vecdiag(V)
	matrix V=V[1,1..2]
	matrix V2=(nullmat(V2)\V)	

    qui ivregress 2sls D.y D.L2.y D.x D.u (D.L.y=L2.y), noconstant
    matrix B=e(b)
    matrix B=B[1,2..3]
    matrix B3=(nullmat(B3)\B)
	matrix V=e(V)
	matrix V=vecdiag(V)
	matrix V=V[1,1..2]
	matrix V3=(nullmat(V3)\V)		
	
    restore

}


* Fig 2

qui drop if year<41
qui getmata (gfs*)=gfmat
egen av_fit = rowmean(gfs*)
egen sd_fit = rowsd(gfs*)
gen  c =  invnormal(1 - (100 - 95) / 200)
gen low = av_fit - c * sd_fit
gen up = av_fit + c * sd_fit

	twoway (rarea low up u, sort(u) color(gs7)) ///
	   (line av_fit u, sort(u) color(black) lpattern(solid)) ///
	   (line gf u, color(gs10) sort lpattern(longdash)), ///
	   legend(label(1 confidence interval at 95%) label(3 DGP) label(2 average fit) ///
	   cols(3) order(3 1 2)) xtitle("U", height(5)) ytitle("g(U)", height(5)) sch(sj) 
graph2tex , epsfile(fig2) caption(Average fit of g(U)) label(fig2) 
	
	
* Table 2
clear
set obs 500
mat res2=J(3,8,.)
forv k=1/3{
    qui svmat B`k'
    su B`k'1,meanonly
    mat res2[`k',1]=r(mean)-0.3
    su B`k'2,meanonly
    mat res2[`k',2]=r(mean)-0.3

	svmat V`k'
	qui replace V`k'1=sqrt(V`k'1)
	qui gen B`k'1_lb=B`k'1-invnormal(0.975)*V`k'1
	qui gen B`k'1_ub=B`k'1+invnormal(0.975)*V`k'1
	
	qui gen CIlen`k'1=B`k'1_ub-B`k'1_lb
	qui su CIlen`k'1,d
	mat res2[`k',5]=r(p50)	
	
	qui replace V`k'2=sqrt(V`k'2)
	qui gen B`k'2_lb=B`k'2-invnormal(0.975)*V`k'2
	qui gen B`k'2_ub=B`k'2+invnormal(0.975)*V`k'2
	
	qui gen CIlen`k'2=B`k'2_ub-B`k'2_lb
	qui su CIlen`k'2,d
	mat res2[`k',6]=r(p50)	
	
	qui gen cov`k'1=(B`k'1_lb<=0.3 & B`k'1_ub>=0.3) 
	su cov`k'1,meanonly
	mat res2[`k',7]=r(mean)	
	
	qui gen cov`k'2=(B`k'2_lb<=0.3 & B`k'2_ub>=0.3) 
	su cov`k'2,meanonly
	mat res2[`k',8]=r(mean)				
	
    qui replace B`k'1=(B`k'1-0.3)^2
    su B`k'1,meanonly
    mat res2[`k',3]=r(mean)
    qui replace B`k'2=(B`k'2-0.3)^2
    su B`k'2,meanonly
    mat res2[`k',4]=r(mean)

}
mat list res2
outtable using res2, mat(res2)  format(%9.4f) replace

sjlog close, replace




