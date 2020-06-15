{smcl}
{* *! version 1.0 21Aug2018}{...}
{cmd:help ivxtplfc} 
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:ivxtplfc} {hline 2}}Partially Linear Functional-Coefficient Panel Data Models with endogeneous variables{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmd:ivxtplcf} {varlist} {cmd:,} {cmdab:z:vars(}{varlist}{cmd:)} {cmdab:u:vars(}{varlist}{cmd:)} {cmdab:gen:erate(}{it:prefix}{cmd:)} 
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent :* {cmdab:z:vars(}{it:{help varlist:varlist}}{cmd:)}}specify variables that have functional coefficients.{p_end}
{p2coldent :* {cmdab:u:vars(}{it:{help varlist:varlist}}{cmd:)}}specify variables that enter into the functions.{p_end}
{p2coldent :* {cmdab:gen:erate(}{it:prefix}{cmd:)}}specify a prefix for the names to store fitted values of functional coefficients. {p_end}
{synopt : {cmd:endox(}{it:{help varlist:varlist}}{cmd:)}}specify endogeneous variables which enter into the model linearly. {p_end}
{synopt : {cmdab:endoz:flag(numlist)}}specify the orders of variables in {cmd:zvars(}{it:{help varlist:varlist}}{cmd:)} which are endogeneous variables. {p_end}
{synopt : {cmd:ivx(}{it:{help varlist:varlist}}{cmd:)}}specify instrumental variables which enter into the model linearly. {p_end}
{synopt : {cmd:ivz(}{it:{help varlist:varlist}}{cmd:[...])}}specify instrumental variables entering into the model nonlinearly. {p_end}
{synopt :{opt te}}specify including time fixed effects.{p_end}
{synopt :{cmd: power(numlist)}}specify the power (or degree) of the splines for the functions.{p_end}
{synopt:{cmdab:nk:nots(numlist)}}specify the number of knots used for the spline interpolation. {p_end}
{synopt:{cmdab:quan:tile}}specify creating knots based on empirical quantiles. {p_end}
{synopt:{cmdab:maxnk:nots(numlist)}}specify the maximun number of knots used for performing LSCV. {p_end}
{synopt:{cmdab:minnk:nots(numlist)}}specify the minimun number of knots used for performing LSCV. {p_end}
{synopt :{cmd: grid(string)}}specify the name for storing the grid points of the variable specified by {cmd:uvar(varname)}.{p_end}
{synopt :{opt pctile(#)}}specify the domain of the generating grid points. The default is {cmd:pctile(0)}. {p_end}
{synopt :{opt brep(#)}}specify the number of bootstrap replications. The default is {cmd:bootstrap(200)}.{p_end}
{synopt :{opt wild}}specify using the wild bootstrap. By default, residual bootstrap with cluster(panelvar) is performed.{p_end}
{synopt :{opt predict(prspec)}}store predicted values of
the conditional mean and fixed effects using variable names specified in {it:prspec}{p_end}
{synopt :{opt nodots}}suppress iteration dots.{p_end}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}.{p_end}
{synopt :{opt fast}}speed up using mata functions.{p_end}
{synopt :{opt tenfoldcv}}specify using ten-fold CV instead of LSCV.{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}* zvars(),uvars(), and generate() are required.{p_end}
{p 4 6 2}Note: Unbalanced panel data must be rectangularized in advance; use {help tsfill:tsfill},full.{p_end}


{title:Description}

{pstd}{opt ivxtplfc} estimates Zhang and Zhou's (2018) partially linear functional-coefficient panel data models with endogeneous variables. 

{pstd}The model can be expressed as:

{space 14}          Y_it=X_it'*\beta+Z_it'*g(u_it)+a_i+e_it

{pstd}where subscripts i and t present individual and time period, respectively; X_it 
and Z_it are vectors of covariates, respectively; G(U_it) is a vector of functional coefficients, 
and U_it is a vector of continuous variables (Specifically, Z_it'*G(U_it)=\sum(Z_jit*G_j(U_jit))); 
a_i represents the fixed effects; e_it is the idiosyncratic error. Some covariates in X_it and Z_it are allowed to be endogeous. 

{pstd}We approximate the functional coefficients by a linear combination of B-spline base functions (see {it:{help bspline:bspline}}). The fixed effects is removed by the first time difference. Then the transformed model is estimated through 
the two step least square (2SLS) technique (see {it:{help ivregress:ivregress}}). 


{title:Options}

{phang}{cmd:zvars(}{it:{help varlist:varlist}}{cmd:)} specifies the variables that have functional coefficients. {cmd:zvars()}is required.

{phang}{cmd:uvars(}{it:{help varlist:varlist}}{cmd:)} specifies (continuous) variables that enter into 
the functional coefficients interacted with variables specified by zvars() in order. {cmd:uvars()}is required.

{phang}{cmd:generate(}{it:prefix}{cmd:)} specifies a prefix for the variable names to store fitted values of functional coefficients. {cmd:generate()}is required.

{phang}{cmd:endox(}{it:{help varlist:varlist}}{cmd:)} specifies endogeneous variables which enter into the model linearly. 

{phang}{cmd:endozflag(}{it:{help numlist:numlist}}{cmd:)} specify the orders of variables in {cmd:zvars(}{it:{help varlist:varlist}}{cmd:)}
 which are endogeneous variables. For example, {cmd:endozflag(1 3)} indicates that the first and third variables
  specified in {cmd:zvars(}{it:{help varlist:varlist}}{cmd:)} are endogeneous.

{phang}{cmd:ivx(}{it:{help varlist:varlist}}{cmd:)}specifies instrumental variables which enter into the model linearly. 

{phang}{cmd:ivz(}{it:{help varlist:varlist}} {cmd:,} {cmdab:u:flag(numlist)} {cmd:[ivtype(numlist)])}specify 
instrumental variables entering into the model nonlinearly which interacts with the functions specified 
by the orders in {cmd:uflag(numlist)}.
Optionally, one may specify the type of nonlinear
 instrumental variables to be constructed. ivtype(#) means using the #th lag of the basis functions and 
 the final IVs are formed from "ivz*L#.S(u)" 
 (where S(u) are basis functions of u). By default, the first lag of the basis functions are used.

{phang}{opt te} specifies including time fixed effects.

{phang}{opt power(numlist)} (nonnegative integers) specifies the power (or degree) of the splines in order specified by uvars(). If absent, 3 is assumed for all the functions.

{phang}{opt nknots(numlist)} specifies the number of knots used for the spline interpolation in order specified by uvars(). 
If absent, 2 is assumed for all the functions. 

{phang}{opt quantile} specifies creating knots based on empirical quantiles. By default, the knots are generated by the rule of equal space.

{phang}{opt maxnknots(numlist)} specifies the maximun number of knots used for performing LSCV. If present, 
LSCV is employed to determine the optimal number of knots. In our practice, we perform the 
Leave-One-Out CV across the panelvar. That is to say, we leave one individual (with all observations during the sample period) out each time. 

{phang}{opt minnknots(numlist)} specifies the minimun number of knots used for performing LSCV. If absent, 2 is assumed.

{phang}{opt grid(string)} specifies the name for storing the grid points of the variable specified by {cmd:uvar(varlist)}. If present, the functional coefficients are estimated over the grid points. By default, 
they are estimated over the observations.

{phang}{opt pctile(#)} specifies the domain of the generating grid points. It can be only used when {cmd: grid(string)} is specified. The default is pctile(0).

{phang}{opt brep(#)} specifies the number of bootstrap replications. The default is brep(200). We recommend that you select the number of replications.{p_end}

{phang}{opt wild} specifies using the wild bootstrap. By default, residual bootstrap with cluster(panelvar) is performed.

{phang}{opt predict(prspec)} stores predicted values of
the dependent variable and fixed effects using variable names specified in {it:prspec}. {it:prspec} is the following:

{phang2}
{cmd:predict(}{varlist}|{it:stub}{cmd:*} [{cmd:, replace noai}]{cmd:)}

{pmore}
The option takes a variable list or a {it:stub}.  The first variable name
corresponds to the predicted outcome mean. The second name corresponds to fixed effects. 

{pmore}
When {cmd:replace} is used, variables with the names in {it:varlist} or
{it:stub}{cmd:*} are replaced by those in the new computation.  If
{cmd:noai} is specified, only a variable for the mean
is created.  

{phang}{opt nodots} suppress iteration dots.

{phang} {opt level(#)} set confidence level; default is {cmd:level(95)}.

{phang}{opt fast} speeds up using mata functions.

{phang}{opt tenfoldcv} specifies using ten-fold CV instead of LSCV.  


{title:Requirement}

{pstd}{opt ivxtplfc} can only be used if data are declared to be panel
data through the {helpb xtset} or {helpb tsset} command.  Before using 
{opt ivxtplfc}, you must install {opt bspline} and {opt moremata}.


{title:Example 1: Endogeneous variable as a linear covariate}

{pstd}Setup{p_end}
{phang2}{cmd:.} {bf:{stata "clear"}}{p_end}
{phang2}{cmd:.} {bf:{stata "set obs 100"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen a=rnormal()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen id=_n"}}{p_end}
{phang2}{cmd:.} {bf:{stata "expand 200"}}{p_end}
{phang2}{cmd:.} {bf:{stata "bys id: gen year=_n"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen x=10*runiform()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen z=5+2*rnormal()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen u=-3+20*runiform()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "mat C=(1,0.4\ 0.4,1)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "drawnorm e1 e2, corr(C)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "xtset id year"}}{p_end}
{phang2}{cmd:.} {bf:{stata "bys id (year): replace x=0.8*x[_n-1]+0.4*e1 if _n>1"}}{p_end}
{phang2}{cmd:.} {bf:{stata "replace x=x+0.3*a"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen gf1=sin(_pi/3*u)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen y=a+5*x+z*gf1+sqrt(2)*e2"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen lag2x=L2.x"}}{p_end}
{phang2}{cmd:.} {bf:{stata "drop if year<151 "}}{p_end}

{pstd}Fixed-effects sieve 2SLS estimation{p_end}

{phang2}{cmd:.} {bf:{stata "ivxtplfc y x, zvars(z) uvars(u) gen(g) endox(x) ivx(lag2x) maxnknots(20)"}}{p_end}

{pstd}Compute the 95% CI for the functional coefficient {p_end}

{phang2}{cmd:.} {bf:{stata "gen lb=g_1-1.96*g_1_sd"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen ub=g_1+1.96*g_1_sd"}}{p_end}

{pstd}Plot the fitted values of the functional coefficient{p_end}

{phang2}{cmd:.} {bf:{stata "local plot1 line gf1 u, sort"}}{p_end}
{phang2}{cmd:.} {bf:{stata "local plot2 line g_1 u, sort"}}{p_end}
{phang2}{cmd:.} {bf:{stata "local plot3 rarea lb ub u, color(gs12) sort"}}{p_end}
{phang2}{cmd:.} {bf:{stata `"twoway (`plot3') (`plot1') (`plot2' legend(label(1 "95% CI") label(2 "Real values") label(3 "Fitted values")))"'}}{p_end}

{title:Example 2: Endogeneous variable as a nonlinear covariate}

{pstd}Setup{p_end}
{phang2}{cmd:.} {bf:{stata "clear"}}{p_end}
{phang2}{cmd:.} {bf:{stata "set obs 100"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen a=rnormal()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen id=_n"}}{p_end}
{phang2}{cmd:.} {bf:{stata "expand 200"}}{p_end}
{phang2}{cmd:.} {bf:{stata "bys id: gen year=_n"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen x=rnormal()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen z=rnormal()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen u=-3+20*runiform()"}}{p_end}
{phang2}{cmd:.} {bf:{stata "mat C=(1,0.4\ 0.4,1)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "drawnorm e1 e2, corr(C)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "xtset id year"}}{p_end}
{phang2}{cmd:.} {bf:{stata "bys id (year): replace x=0.8*x[_n-1]+e1 if _n>1"}}{p_end}
{phang2}{cmd:.} {bf:{stata "replace x=x+0.2*a"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen gf1=3*sin(_pi/3*u)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen y=a+2*z+x*gf1+sqrt(2)*e2"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen lag2x=L2.x"}}{p_end}
{phang2}{cmd:.} {bf:{stata "drop if year<151 "}}{p_end}


{pstd}Fixed-effects sieve 2SLS estimation{p_end}

{phang2}{cmd:.} {bf:{stata "ivxtplfc y z, zvars(x) uvars(u) gen(g) endozflag(1) maxnknots(20) ivz(lag2x,uflag(1)) brep(500) fast"}}{p_end}

{pstd}Compute the 95% CI for the functional coefficient {p_end}

{phang2}{cmd:.} {bf:{stata "gen lb=g_1-1.96*g_1_sd"}}{p_end}
{phang2}{cmd:.} {bf:{stata "gen ub=g_1+1.96*g_1_sd"}}{p_end}

{pstd}Plot the fitted values of the functional coefficient{p_end}

{phang2}{cmd:.} {bf:{stata "local plot1 line gf1 u, sort"}}{p_end}
{phang2}{cmd:.} {bf:{stata "local plot2 line g_1 u, sort"}}{p_end}
{phang2}{cmd:.} {bf:{stata "local plot3 rarea lb ub u, color(gs12) sort"}}{p_end}
{phang2}{cmd:.} {bf:{stata `"twoway (`plot3') (`plot1') (`plot2' legend(label(1 "95% CI") label(2 "Real values") label(3 "Fitted values")))"'}}{p_end}
 

{title:Saved results}

{pstd}
{cmd:ivxtplfc} saves the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(r2)}}within R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted within R-squared{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(mss)}}model sum of squares{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}


{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ivxtplfc}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(ivlist)}}instrumental variables used for estimation.{p_end}
{synopt:{cmd:e(vcetype)}}type of variance-covariance.{p_end}
{synopt:{cmd:e(estfun)}}variables storing the estimated functional coefficients{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(model)}}Fixed-effect sieve 2SLS Estimation{p_end}
{synopt:{cmd:e(k#)}}list of knots for the #th function{p_end}


{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(nknots)}}number of knots{p_end}
{synopt:{cmd:e(power)}}power (or degree) of splines{p_end}
{synopt:{cmd:e(b)}}coefficient vector in the linear part{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators in the linear part{p_end}
{synopt:{cmd:e(bs)}}coefficient vector in the approximating model{p_end}
{synopt:{cmd:e(Vs)}}variance-covariance matrix of the estimators in the approximating model{p_end}


{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
An, Y., H. Cheng, D. Li. 2016. Semiparametric estimation of partially linear varying coefficient panel data models. {it: Advances in Econometrics} 36: 47-65.

{phang}Baltagi, B. H., and D. Li.  2002.  Series estimation of partially
linear panel data models with fixed effects. 
{it:Annals of Economics and Finance} 3: 103-116.

{phang}Du, K., Zhang, Y. and Zhou, Q. 2019. {browse "https://github.com/kerrydu/xtplfc_Stata/blob/master/manuscript.pdf":Estimating partially linear functional-coefficient panel data models with Stata}. 
{it:Working Paper} 

{phang}
Libois, F. and V. Verardi. 2013. Semiparametric fixed-effects estimator. {it:Stata Journal} 13: 329-336.

{phang}
Newson, R.  2000.  B-splines and splines parameterized by their
values at reference points on the x-axis. 
{browse "http://www.stata.com/products/stb/journals/stb57.pdf":{it:Stata Technical Bulletin} 57}: 20-27.  Reprinted in 
{it:Stata Technical Bulletin Reprints}, vol. 10, pp. 221-230. 

{phang}
Zhang, Y., and Q. Zhou. 2018. Partially linear functional-coefficient panel data models: Sieve Estimation and Specification testing. {it:Working paper}.



{title:Author}

{pstd}
Kerry Du{break}
School of Management{break}
Xiamen University{break}
Xiamen, China{break}
{browse "mailto:kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}

{pstd}
Yonghui Zhang{break}
School of Economics{break}
Renmin University of China{break}
Beijing, China{break}
{browse "mailto:yonghui.zhang@hotmail.com":yonghui.zhang@hotmail.com}

{pstd}
Qiankun Zhou{break}
Department of Economics{break}
Renmin University of China{break}
Baton Rouge,USA{break}
{browse "mailto:qzhou@lsu.edu":qzhou@lsu.edu}


{title:Also see}

{p 7 14 2}Help:  {helpb xtsemipar}, {helpb bspline} (if installed){p_end}
