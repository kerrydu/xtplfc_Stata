*! version 7.0.1
* 2020-6-15
*! version 7.0
* 2019-8-26
* change cvtwo option to tenfoldcv
* bugs fix

*! version 6.0.5 
* 2019-01-25
* 2018-12-11
* 2018-11-16
*! version 6.0.4 add time fixed effect option
*  Kerry Du, kerrydu@xmu.edu.cn
cap program drop xtplfc
program define xtplfc, eclass prop(xt)
		
	version 14
	
if replay()& "`e(cmd)'"=="xtplfc" {
	ereturn  display
	exit
}
	
syntax varlist, Zvars(varlist) Uvars(varlist) GENerate(string) [ ///
				    NKnots(numlist integer >=2)  power(numlist integer >0)  ///
				   QUANtile grid(string) pctile(integer 0) brep(integer 200) ///
				   level(integer 95)  WILD  predict(string) TE ///
				  MINNKnots(numlist integer >=2) MAXNKnots(numlist integer >=2) NODOTS FAST TENfoldcv]


				  
*check the required packages
****************************************************************
	capture which lmoremata.mlib
	local rc=_rc
	if (`rc') {
		display as error "{bf:xtplfc} requires {bf:moremata}"
		exit 198
		}
		
	
    capture which bspline
	local rc = _rc
	if (`rc') {
		display as error "{bf:xtplfc} requires {bf:bspline}"
		exit 198
		}

  	qui mata mata mlib index		  
				  
*********************************************************************				  
if "`fast'"!=""{
	xtplfc_fast `varlist', zvars(`zvars') uvars(`uvars') generate(`generate') `te' ///
	                       nknots(`nknots') power(`power') grid(`grid') pctile(`pctile') `quantile'    ///
	                       brep(`brep') level(`level') `wild'  predict(`predict') ///
	                       minnknots(`minnknots') maxnknots(`maxnknots') `nodots' `tenfoldcv'
	exit
}				   


*********************************************************************

	
	local vcetype=cond("`wild'"=="","Clustered Bootstrap","Wild Bootstrap")
	
	
*************************************************************************	
	
	if  `"`predict'"'!=""{
		
		cap checkpredict `predict' 
		if _rc{
			disp as error "Errors in predict( ),"
			checkpredict `predict' 
		}
		local pname  `r(pname)'
		local replace   `r(replace)'
		local noai       `r(noai)'
		
		
	}	
	
	if `"`predict'"'!=""& "`replace'"=="" {
			confirm new variable `pname'
		}
		
	if "`grid'"=="" & `pctile'!=0 {
	
		display as error "{bf:pctile(#)} should be combined with {bf:grid(newvarname)}."
		exit
	}
		
	
****************************************************************************

	if "`maxnknots'"!="" & "`nknots'"!="" {
		di "{err}specify {bf:maxnknots(#)} or {bf:nknots(#)}, "	///
			"not both"
		exit		
	}
	
***************************************************************************	

	gettoken depvar xvars: varlist
	local xvars: list uniq xvars
	local comvars: list zvars & xvars
	if !(`"`comvars'"'==""){
		disp as red "`comvars' included in both linear and nonlinear variable lists."
		display as error "For identification, there should be no common variables between Xvars and Zvars." 
		exit
	}
	
	local nx: word count `xvars'
	local nz: word count `zvars'
	local nu: word count `uvars'
	
	if `nz'!=`nu'{
		di as error "The number of variables in zvars() and uvars() are not equal."
		exit
	}
	
********************************************************************************	
* powermat nknotsmat maxnknotsmat minnknotsmat nknots4cv nsplinesmat


	tempname powermat nknotsmat maxnknotsmat minnknotsmat nsplinesmat mxmi nknots4cv
	if "`power'"==""{
		//disp "By default, the number of power (degree) is set to 3 for all functions."
		mat `powermat'=J(1,`nu',3)
	}
	else{
		local npow: word count `power'
		mata: powermat=strtoreal(tokens(st_local("power")))
		mata: st_matrix("`powermat'",powermat)
		if `npow'!=`nu'{
			di as error "# of integers in power() != # of variables in uvars()."
			exit		
		}
		
	}
	
	
	if "`maxnknots'"=="" & "`nknots'"!="" {
		mata: nknotsmat=strtoreal(tokens(st_local("nknots")))
		mata: st_matrix("`nknotsmat'",nknotsmat)
		local nk0: word count `nknots'
		if `nk0'!=`nz'{
			di as error "# of integers in knots() != # of variables in uvars()."
			exit
		}

	}
	
	
	if "`maxnknots'"=="" & "`nknots'"=="" {
		mat `nknotsmat'=J(1,`nz',2)
	}
	
	if "`maxnknots'"!=""{
		local nmaxnk: word count `maxnknots'
		if `nmaxnk'!=`nz'{
			di as error "# of integers in maxnknots != # of variables in uvars()."
			exit
		}
		mata: maxnknotsmat=strtoreal(tokens(st_local("maxnknots")))
		mata: st_matrix("`maxnknotsmat'",maxnknotsmat)
	
	}
	
	
	if "`minnknots'"!=""&"`maxnknots'"==""{
		di as error "To use LSCV, maxnknots() should be specified."
		exit
	
	}
	
	
	if "`minnknots'"!=""{
		local nminnk: word count `minnknots'	
		if `nminnk'!=`nz'{
			di as error "# of integers in minnknots != # of variables in maxnknots()."
			exit
		}

		mata: minnknotsmat=strtoreal(tokens(st_local("minnknots")))
		mata: st_matrix("`minnknotsmat'",minnknotsmat)
		mata: st_numscalar("`mxmi'",minnknotsmat<=maxnknotsmat)
		if `mxmi'!=1{
			di as error "Minimum # of knots > Maximum # of knots."
			di as error "Check minnknots() and maxnknots()."
			exit
		}
	}
	else{
		mata: minnknotsmat=J(1,`nz',2)
		mata: st_matrix("`minnknotsmat'",minnknotsmat)
	}
	

	if "`maxnknots'"!=""{

		mat   `nsplinesmat'=`maxnknotsmat'+`powermat'
		mata: nknots4cv=combmat(maxnknotsmat,minnknotsmat)
		mata: st_matrix("`nknots4cv'",nknots4cv)
	}
	else{
		mat   `nsplinesmat'=`nknotsmat'+`powermat'

	}

		
local qflag=cond("`quantile'"=="",1,0)	
//disp `qflag'
	
*********************************************************************	
  local vlabel=cond("`grid'"=="","observed","grid")
  local gridname
  * check new generating vars
  if `"`grid'"'!=""{
	local k=1
	foreach uk of local uvars{
		confirm new variable `grid'_`k'
		qui gen `grid'_`k'=.
		if `pctile'==0{
			local upctile = 100 - `pctile'
			su `uk',meanonly
			local gmin=r(min)
			local gmax=r(max)
			mata:gengrid(`gmin',`gmax',`=_N',"`grid'_`k'")
			}
		else{
			local upctile = 100 - `pctile'
			_pctile `uk', p(`pctile' `upctile')
			local gmin=r(r1)
			local gmax=r(r2)
			mata:gengrid(`gmin',`gmax',`=_N',"`grid'_`k'")
		}
		label var `grid'_`k' "Grid points of `uk' over [`pctile'-`upctile'] percentiles"
		local gridname `gridname' `grid'_`k'	
		local k=`k'+1	
	}
  }
  
  
    local uvarscopy `uvars'
    local nn=1
  	foreach zk of local zvars{
  		gettoken uk uvarscopy:uvarscopy
		confirm new variable `generate'_`nn'
		qui gen `generate'_`nn'=.
		label var `generate'_`nn' `"Estimated g(`uk') at the `vlabel' points w.r.t `zk'*g(`uk')"'
		confirm new variable `generate'_`nn'_sd
		qui gen `generate'_`nn'_sd=.
		label var  `generate'_`nn'_sd  `"`vcetype' S.E. of g(`uk') w.r.t `zk'*g(`uk')"'
		local nn=`nn'+1
		
	}
  
  
**********************************************************************

	
	*******************************
	/*
	_xt, trequired 
	local id=r(ivar)
	local time=r(tvar)
	
	/*
	//fill the gap
	tempvar dataflag
	qui gen `dataflag'=1
	
	qui tsfill,full
    */
	//continuous id and time
	tempvar cid ctime
	qui egen `cid'=group(`id')
	qui egen `ctime'=group(`time')
	
	// Rectangularize the data
	// Note: For survey data, there might be year gaps which must not be filled in.
	// re-xtset before using tsfill
	qui xtset `cid' `ctime' 
	tempvar dataflag
	qui gen `dataflag'=1
	qui tsfill,full	
	*/
	//qui tab `cid', nofreq
	//local nid=r(r)
	//qui tab `ctime',nofreq
	//local ntime=r(r)

	 qui xtset
	 if r(balanced)!="strongly balanced" {
		noi display as error "{bf:xtplfc} needs strongly balanced panel, use {bf:tsfill, full} to rectangularize your data!}."
		exit	 	
	 }


	local id=r(panelvar)    
    local time=r(timevar) 
    local tdelta=r(tdelta)   
	//continuous id and time
	tempvar cid ctime
	qui egen `cid'=group(`id')
	qui egen `ctime'=group(`time')
	
	qui count if `ctime'==`ctime'[1]
	local nid=r(N)
	qui count if `cid'==`cid'[1]
	local ntime=r(N)
	
	local te=cond("`te'"=="",0,`ntime'-1)
	
	sort `id' `time'
	/*
	qui count if `time'-`time'[_n-1]!=`tdelta' & `cid'==`cid'[_n-1]
	local gaps=r(N)
	if `gaps'!=0{
		noi di as error `"`time' is regularly spaced, but does not have intervals of 1"'
		exit
	}
	*/
	qui count if `time'-`time'[_n-1]!=`tdelta' & `cid'==`cid'[_n-1]
	local gaps=r(N)
	if `gaps'!=0{
		noi di as error `"`time' is regularly spaced, but does not have intervals of 1"'
		exit
	}  	
	

	if `te'>0{
		mata: tematrix=gendtime(`nid',`te')
		mata: _sttempvar("tematrix",tematrix)
	}
	
**********************************
**********generating tempvar for regression

	
	local Dxvars
	
	foreach v of local xvars {

		local Dxvars `Dxvars' D.`v'
		
	}
	
	local Dxvars2 `Dxvars' `tematrix'  // adding time fixed effets
	
	local nx2=`nx'+`te'
   
   
   local labtematrix
   forv j=1/`te'{
		local temp=`time'[`j']
		local labtematrix `labtematrix' _I`time'_`temp'
   }
   
   local xvars2 `xvars' `labtematrix'
   
	


  forv k=1/`nz'{
	local nspline=`nsplinesmat'[1,`k']
	forv j=1/`nspline'{
	   tempvar  sp`k'_`j'
	}
  }
  
  
 /* 
  local k=1
  foreach uk of local uvars{
	su `uk', meanonly
	local min_`k'=r(min)
	local max_`k'=r(max)
	local k=`k'+1
  
  }
*/	

	if "`nodots'"!=""{
		local nodots qui
	}
	
qui putmata U=(`uvars'),replace	

	
***********Using CV to determine the optimal number of knots...******	
/////////////////////////////////////////////////////////////////////
if "`maxnknots'"!=""{
	
		if "`tenfoldcv'"!=""{
				// Faster computation for CV
				// The idea is originated from the npivcv command by Dongwoo Kim (University College London)
				local rngstate = c(rngstate) // record the ranom state
				set seed 1004
				tempvar splitdummy				
                qui splitsample if `ctime'==`ctime'[1], gen(`splitdummy') n(10)
				qui bys `cid' (`ctime'): replace `splitdummy'=`splitdummy'[_n-1] if _n>1

				set seed `rngstate' 	
			
		}		
	
	
		di
		tempvar ehat_sq  yhat 
		tempname mse 		
		disp as green "Using Cross-Validation to determine the optimal number of knots..."
		local nnknots=rowsof(`nknots4cv')
		mat `mse'=J(`nnknots',1,.)
		qui gen `ehat_sq'=.
	
	//qui xtset
	forvalues kk=1/`nnknots'{ 
	
	     **generating spline bases
		 local splines	
		 forv zk=1/`nz' {
			local uk: word `zk' of `uvars'
			local nknot`zk'=`nknots4cv'[`kk',`zk']
			local power_`zk'=`powermat'[1,`zk']
			//local udist=(`max_`zk''-`min_`zk'')/`nknot`zk''
			mata: _gknotlist("klist`zk'", `qflag', U[.,`zk'], `nknot`zk'',`power_`zk'')
			local sp`zk'
			local nspline=`nknot`zk''+`power_`zk''
				forv zj=1/`nspline'{
					local sp`zk' `sp`zk'' `sp`zk'_`zj''
				}

			//qui bspline `sp`zk'', xvar(`uk') knots(`min_`zk''(`udist')`max_`zk'') power(`power_`zk'')	
			qui bspline `sp`zk'', xvar(`uk') knots(`klist`zk'') power(`power_`zk'')	noexknot
			local splines `splines' `sp`zk''
		  }
		  
		  
		
		local HZ
		local lHZ
		local k=0
		foreach zk of local zvars {
			local j=0
			local k=`k'+1
			foreach uk of local sp`k' {
				local j=`j'+1
				tempvar `zk'`j'
				qui gen ``zk'`j''=`zk'*`uk' 
			    local HZ `HZ' D.``zk'`j''
				local lHZ `lHZ' `zk'`j'
			}
		}
		
		
		local allreg  `Dxvars2' `HZ' 
	
		if "`tenfoldcv'"!=""{
			qui xtset `id' `time'
			 forvalues jt=1/10 {
				 qui reg D.`depvar' `allreg'  if `splitdummy'~=`jt', noconstant
				 qui predict `yhat' if `splitdummy'==`jt',xb
				 qui replace `ehat_sq'=(D.`depvar'-`yhat')^2 if `splitdummy'==`jt'
				 qui cap drop `yhat'
			 }		
			
		}
		else{
             //Leave-One-Out CV
			 forvalues jt=1/`nid' {
				 qui reg D.`depvar' `allreg'  if `cid'~=`jt', noconstant
				 qui predict `yhat' if `cid'==`jt',xb
				 qui replace `ehat_sq'=(D.`depvar'-`yhat')^2 if `cid'==`jt'
				 qui cap drop `yhat'
			 }
				
		 }		

		qui su `ehat_sq' if ~missing(`ehat_sq')
		mat `mse'[`kk',1]=r(sum)/r(N)

		cap drop `splines'
		cap drop `lHZ'
		`nodots' displaydot, dot(`kk')

	}
		//mat list `mse'
		mata: minpos=_minpos(st_matrix("`mse'"))
		mata: nknotsmat=nknots4cv[minpos,.]
		mata: st_matrix("`nknotsmat'",nknotsmat)
		
		di
		local i=1
		foreach uj of local uvars{
			local ujk: word `i' of `zvars'
			disp as green "The optimal number of knots w.r.t. `ujk'*g(`uj') is  " `=`nknotsmat'[1,`i']' "."
			local i=`i'+1
		}
		


}
//////////////////////////////////////////////////////////////////	

		
	local N=_N
	local splines

  mat `nsplinesmat'	=`nknotsmat'+`powermat'
  mata: nsplinesmat=st_matrix("`nsplinesmat'")
  //mata: nsplinesmat=st_matrix("`nsplinesmat'")
  mata: bzindex=0,cumsum(nsplinesmat)
  tempname bzindex	
  mata: st_matrix("`bzindex'",bzindex)
 
 //local uvarscopy `uvars' 
 forv zk=1/`nz' {
 
	local nspline_`zk'=`nsplinesmat'[1,`zk']
	local nknots_`zk'=`nknotsmat'[1,`zk']
	local power_`zk'=`powermat'[1,`zk']
	local uk: word `zk' of `uvars'
	mata: _gknotlist("klist`zk'", `qflag', U[.,`zk'], `nknots_`zk'',`power_`zk'')
	local sp`zk'
	

		forv zj=1/`nspline_`zk''{
			local sp`zk' `sp`zk'' `sp`zk'_`zj''
		}
		qui bspline `sp`zk'', xvar(`uk') knots(`klist`zk'') power(`power_`zk'') noexknot
		local splines `splines' `sp`zk''
		local knotlist_`zk' `r(knots)'
	
}	
	
	
	
	//disp " r(nknot) =`r(nknot)'"

	qui sort `id' `time'
	qui putmata HH=(`splines'),replace
	
	local HZ
	local lHZ
	local labHZ
	 //qui xtset
	//qui xtset `cid' `ctime'
	
	local k=0
	foreach zk of local zvars {
		local k=`k'+1
		local j=0
		local uk: word `k' of `uvars'
		foreach uki of local sp`k' {
			local j=`j'+1
			tempvar `zk'`j' 
			qui gen ``zk'`j''=`zk'*`uki'
		    local HZ `HZ' D.``zk'`j'' 
		    local lHZ `lHZ' ``zk'`j''
		    local labHZ `labHZ' `zk'*S(`uk')_`j'
		}
	}
	


***********************************************************************************

	local allreg  `Dxvars2' `HZ' 
	mata: Xalls  =  tokens(st_local("Dxvars2"))
	mata: Zalls  =  tokens(st_local("HZ"))
	mata: allreg =  tokens(st_local("allreg"))
	mata: xpos   =  posof(Xalls,allreg)
	mata: zpos   =  posof(Zalls,allreg)
   
   	qui reg D.`depvar' `allreg' , noconstant 

	local dof    =   e(df_r)
	local r2     =   e(r2)
	local rmse   =   e(rmse)
	local mss    =   e(mss)
	local rss    =   e(rss)
	local r2_a   =   e(r2_a)
	local nobs   =   e(N)  	

    tempvar touse0 Dxb ehat aui Dehat xb
    qui gen `touse0'=e(sample)
    qui predict `Dehat', residuals
    qui predict `Dxb', xb
	tempname b0 b1 V B b_sp V1 V2
 	mat `b0'=e(b)
	mata: b0=st_matrix("`b0'")

////////////////////////////////
   sort `id' `time'
   mata: XHZ=st_data(.,"`xvars' `lHZ'")
   if `nx'>0{
    mata: b00=b0[1..`nx'],b0[(`nx2'+1)..length(b0)]  	
   }
   else{
   	mata: b00=b0[(`nx2'+1)..length(b0)] 
   }

   mata: xb=XHZ*b00'
   if `te'>0{
   	mata: timeeffects=b0[(`nx'+1)..`nx2'],-sum(b0[(`nx'+1)..`nx2'])
   	mata: xb=xb+J(`nid',1,timeeffects')
   }
   sort `id' `time'
   getmata `xb'=xb,replace

   qui gen `ehat'=`depvar'-`xb'
   qui bys `cid' (`ctime'): egen `aui'=mean(`ehat')
   qui replace `ehat'=`ehat'-`aui'
   qui replace `xb'=`xb'+`aui'
   
    if "`predict'"!=""{
		gettoken v pname: pname
		qui gen `v'=`xb'
		label var `v' "Fitted values of the dependent variable"
		if "`noai'"==""{
		   qui gen `pname'=`aui'
		   label var `pname' "Estimates of the fixed effects"				
		
		}
   }  
   
   mata: st_matrix("`b_sp'",b0[1,zpos])
   mat `B'=`b_sp'
	if "`Dxvars2'"!=""{
		mata: b1=b0[1,xpos]
		mata:st_matrix("`b1'",b1)
		mat `B'=(`b1',`b_sp')
	}


*******************************************************************	

/////////////////////////////////////////////////////////////	
	local vlabel=cond("`grid'"=="","observed","grid")	
    local estfun
    qui sort `id' `time'	
	mata: bz=b0[1,zpos]
	//mata: bz

		
   if "`grid'"!=""{
	    qui sort `id' `time'
		
	 forv zk=1/`nz' {
		local gds_`zk'	
		forv j=1/`nspline_`zk''{
			tempvar gd_`zk'`j'
			local gds_`zk' `gds_`zk'' `gd_`zk'`j''
		}

		quietly bspline `gds_`zk'', xvar(`grid'_`zk')  knots(`klist`zk'') power(`power_`zk'') noexknot
		local gsplist `gsplist' `r(splist)'		
	    qui sort `id' `time'
		qui putmata HH=(`gsplist'),replace
      }
   }
   
 	mata: gfun=estgfun(bz',HH,bzindex)	
	local j=1
	foreach zz of local zvars {

		mata: gtemp=gfun[.,`j']
		qui getmata `generate'_`j'=gtemp,replace		
		local estfun `estfun' `generate'_`j'
		local j=`j'+1
	}


*****************************************************
if `brep'>0 {
		di
		disp "Computing the bootstrap standard errors..."
		//local vcetype "Clustered Bootstrap"
		local bootcmd bs_xtresid, id(`cid') xb(`Dxb') residual(`Dehat') regvars(`allreg') 
		if "`wild'"!=""{
			//local vcetype "Wild Bootstrap"
			local bootcmd bs_xtwild, xb(`xb') residual(`ehat') regvars(`allreg') 
		}
		

	//mata: bsmat=J(`bootstrap',`nx',.)	
	mata: bsmat=J(`brep',length(b0),.)
	tempname bb

	forvalues z = 1/`nz'{
		mata: gf_`z'=J(`N',`brep',.)
	}
		
		qui sort `id' `time'
		//qui xtset
		//qui xtset `cid' `ctime'
		
		forv k=1/`brep' {
			
			`nodots' displaydot, dot(`k')
			qui `bootcmd'
			mat `bb'=e(bs)
			mata: bb=st_matrix("`bb'")
			mata: bs2=bb[1,zpos]
			mata: bsmat[`k',.]=bb[1,xpos],bs2
			
			//mata: bzindex=0,cumsum(nsplinesmat)	
			//mata: st_matrix("`bzindex'",bzindex)
			forvalues z=1/`nz'{
				mata: bz=bs2[1,(1+bzindex[1,`z'])..bzindex[1,`z'+1]]
				mata: gf_`z'[.,`k']=HH[.,(1+bzindex[1,`z'])..bzindex[1,`z'+1]]*bz'

			}
			
			

	   }

			mata: V=quadvariance(bsmat)
			mata: st_matrix("`V'",V)	   
		
			if "`Dxvars2'"!=""{
				mat `V1'=`V'[1..`nx2',1..`nx2']
			}

			* compute the s.e. for the smooth function
			*sort `id' `time'
		   local k=1
		   foreach zk of local zvars{
				mata: gf_m=rowsum(gf_`k')/`brep'
				mata: sd=(gf_`k':-(gf_m#J(1,`brep',1))):^2
				mata: sd=rowsum(sd)/(`brep'-1)
				mata: sd=sd:^0.5
				qui getmata `generate'_`k'_sd=sd,replace
				local k=`k'+1
		   
		   }

			

}
else{
	mat `V'=J(colsof(`b0'),colsof(`b0'),.)
	if "`Dxvars2'"!=""{
		mat `V1'=`V'[1..`nx2',1..`nx2']
	}	
}
/*

else{
   sort `id' `time'
   tempname omega 
   qui putmata flag=(`touse0'),replace

   local vcetype "Robust"
		   
   qui putmata cid=`cid' ehat=`Dehat' if `touse0'==1, replace
   mata: XHZ=st_data(.,"`Dxvars2' `HZ'","`touse0'")
		   

	   if (`"`hacopt'"'==""&"`weightvar'"==""){
			
			mata: omega=lrvar(cid,XHZ,ehat)
			
	   }
	   else{
		   local vcetype "HAC Robust"
		   if ("`weightvar'"!="") qui replace `Dehat'=`Dehat'*`weightvar'
		//qui tsset `cid' `ctime'
		qui xtlrcov `Dxvars2' `HZ', id(`cid') time(`ctime') wvar(`Dehat') `hacopt' 
		mat `omega'=r(omega)
		mata: omega=st_matrix("`omega'")

	   }
	   	//mata: colsum(XHZ)
		mata: invXX=invsym(quadcross(XHZ,XHZ))
	
		mata: V=invXX*omega*invXX // (nT*invXX)*(omega/nT)*(nT*invXX)/nT
		mata: st_matrix("`V'",V)		
		if `"`Dxvars2'"'!=""{
			mat `V1'=`V'[1..`nx2',1..`nx2']

		}
		mat `V2'=`V'[`=`nx2'+1'..rowsof(`V'),`=`nx2'+1'..rowsof(`V')]
		mata: V2=st_matrix("`V2'")
		mata: gsd=gcovar(V2,HH,bzindex)


		   sort `id' `time'
		   local k=1
		   foreach zk of local zvars{
				local uk: word `k' of `uvars'
				mata: tempsd=gsd[,`k']
				qui getmata `generate'_`zk'_sd=tempsd,replace
				//label var `generate'_`zk'_sd  "Robust S.E. of g(`uk') w.r.t `zk'*g(`uk')"
				local k=`k'+1
		   
		   }
		   

	   
	   
}
*/


	mat colnames `B'= `xvars2'  `labHZ'
	mat colnames `V'= `xvars2'  `labHZ'
	mat rownames `V'= `xvars2'  `labHZ'
	
	

	noi disp in green""
	noi disp in green "Partially linear functional-coefficient panel data model."
	noi disp in green "Fixed-effect series semiparametric estimation."

if (`nx2'>0){	
		ereturn local vcetype `vcetype'
		mat colnames `b1'= `xvars2' 
		mat colnames `V1'= `xvars2' 
		mat rownames `V1'= `xvars2' 

			if `brep'>0{
				ereturn post  `b1' `V1', dep(`depvar') obs(`nobs') esample(`flag') buildfvinfo
			}
			else{
				ereturn post  `b1', dep(`depvar') obs(`nobs') esample(`flag') buildfvinfo
			}
		*disp 
		disp "Estimation results: linear part" _c
		//noi di ""
		noi di in green "{col 48} Number of obs        =" in yellow %8.0f `nobs'
		noi di in green "{col 48} Within R-squared     =" in yellow %8.4f `r2'
		noi di in green "{col 48} Adj Within R-squared =" in yellow %8.4f  `r2_a'
		noi di in green "{col 48} Root MSE             =" in yellow %8.4f `rmse'		
		ereturn display, level(`level') 		
	}
	else{
	 ereturn post, dep(`depvar') obs(`nobs') esample(`touse0') buildfvinfo	
	}

		
		noi disp in green "Estimated functional coefficient(s) are saved in the data."
		//ereturn clear
		//ereturn post  `B' `V',esample(`flag') dep(`depvar') obs(`nobs') dof(`dof')
		//ereturn mat    knotv  = `knotv'
		ereturn mat bs=`B'
	    ereturn mat Vs=`V'
		ereturn mat nknots = `nknotsmat'
		ereturn mat power  = `powermat'
		ereturn scalar df_m   = `nobs'-`dof'
		ereturn scalar df_r   = `dof'
		ereturn scalar r2     = `r2'
		ereturn scalar rmse   = `rmse'
		ereturn scalar mss    = `mss'
		ereturn scalar rss    = `rss'
		ereturn scalar r2_a   = `r2_a'
		ereturn scalar N      = `nobs'
		
		forv j=1/`nu'{
			ereturn local k`j' `knotlist_`j''
		}
		
		ereturn local cmd "xtplfc" 
		ereturn local title "Partially linear functional-coefficient panel data model"
		ereturn local model "Fixed-effect Series Semiparametric Estimation" 
		ereturn local depvar "`depvar'" 
		ereturn local estfun `estfun'
		ereturn local basis "B-spline"
		ereturn local vcetype `vcetype'
		//qui keep if `dataflag'==1
		//qui xtset `id' `time'
		sort `id' `time'
	

end

   
   
 **************************************************
 

cap program drop  bs_xtresid
program define bs_xtresid,eclass
        version 14
        syntax, id(varname numeric) xb(varname numeric) RESidual(varname numeric) regvars(varlist ts fv) 

        tempvar idx ys
		* idx randomly selects the observations with replacement
        qui putmata id=`id',replace
        mata: mm_panels(id,idinfo=.)
       // mata: length(idinfo)
        mata: idx=mm_sample(length(idinfo),.,idinfo)
        qui getmata `idx'=idx, replace

        * the new dependent variable using resample residuals
        qui gen double `ys' = `xb' + `residual'[`idx']
		tempname b
		qui reg `ys' `regvars',noconstant
		
        mat `b'=e(b)
        ereturn matrix bs=`b'
end  



cap program drop  bs_xtwild
program define bs_xtwild,eclass
        version 14
        syntax,  xb(varname numeric) RESidual(varname numeric) regvars(varlist ts fv) 

		tempvar randu ys flag Dys
		
		qui gen `randu'=uniform()
		qui gen `flag'=(1-sqrt(5))/2
		qui replace `flag'=(sqrt(5)+1)/2 if `randu'>((1+sqrt(5))/(2*sqrt(5)))
		
		qui gen `ys'=`xb'+`residual'*`flag'
		
		qui reg  D.`ys' `regvars',noconstant
		
		//list D.`ys'
		tempname b
        mat `b'=e(b)
        
        ereturn matrix bs=`b'
end   



**********************************************************************




cap program drop checkhacopt

program define checkhacopt, rclass

	version 14
	syntax, [WVar(varname) *]
	local opt `options'
	local 0 `", `options'"'
	syntax, [noCENTer  ///
		dof(integer 0) ///
		vic(string) vlag(integer 0)  ///
		KERNel(string) BWIDth(real 0) bmeth(string) blag(real 0) bweig(numlist) bwmax(real 0)]

	return local weightvar `wvar'
	return local hacopt `opt'
	
	

end

cap program drop checkpredict

program define checkpredict, rclass

	version 14
	syntax namelist [, REPLACE NOAI]
	local nv: word count `namelist'
	if `nv'>1{
		return local pname `namelist'
	}
	else if "`noai'"==""{
		return local pname `namelist'_yhat `namelist'_ai
	}
	else{
		return local pname `namelist'
	}
	return local replace `replace'
	return local noai `noai'
end

 cap program drop displaydot
 program define displaydot
 version 14
 syntax, dot(integer)
 if mod(`dot',50) == 0 {
	di _c `dot'
	di 
}
else{
  disp _c "."
}
  end
 
/* 

*/

*!This program is copyed from the Stata Blog of David M. Drukker
* https://blog.stata.com/2016/03/17/programming-an-estimation-command-in-stata-making-predict-work/#e(predict)
cap program drop getcinfo
program define getcinfo, rclass
    syntax varlist(ts fv)
    _rmcoll `varlist' , noconstant expand
    local cnames `r(varlist)'
    local p : word count `cnames'
    
    tempname b mo
 
    matrix `b' = J(1, `p', 0)
    matrix colnames `b' = `cnames' 
    _ms_omit_info `b'
    matrix `mo' = r(omit)
 
    return local  cnames "`cnames'"
    return matrix mo = `mo'
end



*!This program is a minor modification of getcinfo written by David M. Drukker
cap program drop getcinfolist
program define getcinfolist, rclass
    syntax varlist(ts fv)
    _rmcoll `varlist', noconstant expand
    local cnames `r(varlist)'
    local p : word count `cnames'
    
    tempname b mo
 
    matrix `b' = J(1, `p', 0)
    matrix colnames `b' = `cnames' 
    _ms_omit_info `b'
    matrix `mo' = r(omit)
    /*
    mata: mflag=st_matrix("`mo'")
 	mata: templist=tokens(st_local("cnames"))
	mata: templist=select(templist,mflag:==0) 
	mata: st_local("cnames",concat(templist," "))
	*/
		
    *return local  cnames "`cnames'"
  rmcollv, vlist(`cnames') mat(`mo') 
  return local  cnames `r(cnames)'
end
***********************
cap program drop rmcollv
program define rmcollv, rclass
syntax, Vlist(string) Mat(string)

mata: rmcollv("`vlist'","`mat'","cnames")

return local cnames `cnames'
end

*********************************************

* 2018-11-14
*! version 5.0.4
* 16Sep2018 revised
*  Kerry Du, kerrydu@sdu.edu.cn
cap program drop xtplfc_fast
program define xtplfc_fast, eclass prop(xt)
		
	version 14
	
syntax varlist, Zvars(varlist) Uvars(varlist) GENerate(string) [ ///
				   NKnots(numlist integer >=2)  power(numlist integer >0)  ///
				   grid(string) pctile(integer 0) brep(integer 200) QUANtile TE ///
				   level(integer 95)  WILD  predict(string) ///
				  MINNKnots(numlist integer >=2) MAXNKnots(numlist integer >=2) NODOTS TENfoldcv]



	local vcetype=cond("`wild'"=="","Clustered Bootstrap","Wild Bootstrap")
	


/*
//check options for lrcov command
	if  `"`hacopt'"'!=""{
		local vcetype HAC robust
		cap checkhacopt, `hacopt' 
		if _rc{
			disp as error "Errors in  hacopt( ),"
			checkhacopt,  `hacopt'
		}
		local hacopt `r(hacopt)'
		local weightvar `r(weightvar)'
		
		
	}
*/
	
*************************************************************************	
	
	if  `"`predict'"'!=""{
		
		cap checkpredict `predict' 
		if _rc{
			disp as error "Errors in predict( ),"
			checkpredict `predict' 
		}
		local pname  `r(pname)'
		local replace   `r(replace)'
		local noai       `r(noai)'
		
		
	}	
	
	if `"`predict'"'!=""& "`replace'"=="" {
			confirm new variable `pname'
		}
		
	if "`grid'"=="" & `pctile'!=0 {
	
		display as error "{bf:pctile(#)} should be combined with {bf:grid(newvarname)}."
		exit
	}
		
	
****************************************************************************

	if "`maxnknots'"!="" & "`nknots'"!="" {
		di "{err}specify {bf:maxnknots(#)} or {bf:nknots(#)}, "	///
			"not both"
		exit		
	}
	
***************************************************************************	

	gettoken depvar xvars: varlist
	local xvars: list uniq xvars
	local comvars: list zvars & xvars
	if !(`"`comvars'"'==""){
		disp as red "`comvars' included in both linear and nonlinear variable lists."
		display as error "For identification, there should be no common variables between Xvars and Zvars." 
		exit
	}
	
	local nx: word count `xvars'
	local nz: word count `zvars'
	local nu: word count `uvars'
	
	if `nz'!=`nu'{
		di as error "The number of variables in zvars() and uvars() are not equal."
		exit
	}
	
********************************************************************************	
* powermat nknotsmat maxnknotsmat minnknotsmat nknots4cv nsplinesmat


	tempname powermat nknotsmat maxnknotsmat minnknotsmat nsplinesmat mxmi nknots4cv
	if "`power'"==""{
		//disp "By default, the number of power (degree) is set to 3 for all functions."
		mat `powermat'=J(1,`nu',3)
	}
	else{
		local npow: word count `power'
		mata: powermat=strtoreal(tokens(st_local("power")))
		mata: st_matrix("`powermat'",powermat)
		if `npow'!=`nu'{
			di as error "# of integers in power() != # of variables in uvars()."
			exit		
		}
		
	}
	
	
	if "`maxnknots'"=="" & "`nknots'"!="" {
		mata: nknotsmat=strtoreal(tokens(st_local("nknots")))
		mata: st_matrix("`nknotsmat'",nknotsmat)
		local nk0: word count `nknots'
		if `nk0'!=`nz'{
			di as error "# of integers in knots() != # of variables in uvars()."
			exit
		}

	}
	
	
	if "`maxnknots'"=="" & "`nknots'"=="" {
		mat `nknotsmat'=J(1,`nz',2)
	}
	
	if "`maxnknots'"!=""{
		local nmaxnk: word count `maxnknots'
		if `nmaxnk'!=`nz'{
			di as error "# of integers in maxnknots != # of variables in uvars()."
			exit
		}
		mata: maxnknotsmat=strtoreal(tokens(st_local("maxnknots")))
		mata: st_matrix("`maxnknotsmat'",maxnknotsmat)
		//mata: nknotsmat=0 // flag to use CV
	
	}
	
	
	if "`minnknots'"!=""&"`maxnknots'"==""{
		di as error "To use LSCV, maxnknots() should be specified."
		exit
	
	}
	
	
	if "`minnknots'"!=""{
		local nminnk: word count `minnknots'	
		if `nminnk'!=`nz'{
			di as error "# of integers in minnknots != # of variables in maxnknots()."
			exit
		}

		mata: minnknotsmat=strtoreal(tokens(st_local("minnknots")))
		mata: st_matrix("`minnknotsmat'",minnknotsmat)
		mata: st_numscalar("`mxmi'",minnknotsmat<=maxnknotsmat)
		if `mxmi'!=1{
			di as error "Minimum # of knots > Maximum # of knots."
			di as error "Check minnknots() and maxnknots()."
			exit
		}
	}
	else{
		mata: minnknotsmat=J(1,`nz',2)
		mata: st_matrix("`minnknotsmat'",minnknotsmat)
	}
	

		
	
	
*********************************************************************	
  local vlabel=cond("`grid'"=="","observed","grid")
  local gridname
  * check new generating vars
  if `"`grid'"'!=""{
	local k=1
	foreach uk of local uvars{
		confirm new variable `grid'_`k'
		qui gen `grid'_`k'=.
		if `pctile'==0{
			local upctile = 100 - `pctile'
			su `uk',meanonly
			local gmin=r(min)
			local gmax=r(max)
			mata:gengrid(`gmin',`gmax',`=_N',"`grid'_`k'")
			}
		else{
			local upctile = 100 - `pctile'
			_pctile `uk', p(`pctile' `upctile')
			local gmin=r(r1)
			local gmax=r(r2)
			mata:gengrid(`gmin',`gmax',`=_N',"`grid'_`k'")
		}
		label var `grid'_`k' "Grid points of `uk' over [`pctile'-`upctile'] percentiles"
		local gridname `gridname' `grid'_`k'	
		local k=`k'+1	
	}
  }
  
  
    local uvarscopy `uvars'
    local nn=1
  	foreach zk of local zvars{
  		gettoken uk uvarscopy:uvarscopy
		confirm new variable `generate'_`nn'
		qui gen `generate'_`nn'=.
		label var `generate'_`nn' `"Estimated g(`uk') at the `vlabel' points w.r.t `zk'*g(`uk')"'
		confirm new variable `generate'_`nn'_sd
		qui gen `generate'_`nn'_sd=.
		label var  `generate'_`nn'_sd  `"`vcetype' S.E. of g(`uk') w.r.t `zk'*g(`uk')"'
		local nn=`nn'+1
		
	}
  
  
**********************************************************************
/*
	_xt, trequired 
	local id=r(ivar)
	local time=r(tvar)
	

	//continuous id and time
	tempvar cid ctime
	qui egen `cid'=group(`id')
	qui egen `ctime'=group(`time')
	
	// Rectangularize the data
	// Note: For survey data, there might be year gaps which must not be filled in.
	// re-xtset before using tsfill
	qui xtset `cid' `ctime' 
	tempvar dataflag
	qui gen `dataflag'=1
	qui tsfill,full	
*/

	 qui xtset
	 if r(balanced)!="strongly balanced" {
		noi display as error "{bf:xtplfc} needs strongly balanced panel, use {bf:tsfill, full} to rectangularize your data!}."
		exit	 	
	 }

	local id=r(panelvar)    
    local time=r(timevar) 
    local tdelta=r(tdelta)   
	//continuous id and time
	tempvar cid ctime
	qui egen `cid'=group(`id')
	qui egen `ctime'=group(`time')	 
	 
	
	qui count if `cid'==`cid'[1]
	local ntime=r(N)
	//disp "ntime=`ntime'"
	local te=cond("`te'"=="",0,`ntime'-1)
	
    sort `id' `time'
   
	qui count if `time'-`time'[_n-1]!=`tdelta' & `cid'==`cid'[_n-1]
	local gaps=r(N)
	if `gaps'!=0{
		noi di as error `"`time' is regularly spaced, but does not have intervals of 1"'
		exit
	}  
   
   local labtematrix
   forv j=1/`te'{
		local temp=`time'[`j']
		local labtematrix `labtematrix' _I`time'_`temp'
   }
   
   local xvars2 `xvars' `labtematrix'	
   local nx2: word count `xvars2'
	
******************************************************************
/*
if "`cvtwo'"!=""{

	// Faster computation for CV
	// The idea is originated from the npivcv command by Dongwoo Kim (University College London)
	local state = c(rngstate) // record the ranom state
	set seed 1004
	tempvar samplesplit splitdummy
	qui gen double `samplesplit' = rnormal(0, 1) if `time'==`time'[1]
	quietly summarize `samplesplit', detail
	local med = r(p50)
	qui gen byte `splitdummy' = (`samplesplit' > `med')	if `time'==`time'[1]
	qui bys `cid' (`ctime'): replace `splitdummy'=`splitdummy'[_n-1] if _n>1
	set rngstate `state'
	

}
*/	
	
		if "`tenfoldcv'"!=""{
				// Faster computation for CV
				// The idea is originated from the npivcv command by Dongwoo Kim (University College London)
				local rngstate = c(rngstate) // record the ranom state
				set seed 1004
				tempvar splitdummy				
                qui splitsample if `ctime'==`ctime'[1], gen(`splitdummy') n(10)
				qui bys `cid' (`ctime'): replace `splitdummy'=`splitdummy'[_n-1] if _n>1

				set seed `rngstate' 	
			
		}	

	
*******************************************************************

	local yX `depvar' `xvars'	
	tempname V b
	local didots=cond("`nodots'"=="",1,0)
	local qflag=cond("`quantile'"=="",1,0)	
    di

/*	
	mata: xtplfc("`cid'","`ctime'","`yX'","`zvars'","`uvars'","`nknotsmat'", ///
	           "`maxnknotsmat'","`minnknotsmat'","`powermat'",`qflag',`te',`brep',"`V'","`b'","`generate'", ///
			   "`gridname'",`didots',"`wild'")
*/	

	mata: xtplfc("`cid'","`ctime'","`yX'","`zvars'","`uvars'","`nknotsmat'", ///
	           "`maxnknotsmat'","`minnknotsmat'","`powermat'",`qflag',`te',`brep',"`V'","`b'","`generate'", ///
			   "`gridname'",`didots',"`wild'","`splitdummy'")
		   
	
	qui replace `touse'=0 if missing(`touse')
	tempname nsplinesmat
	//mat list `nknotsmat'
	mat `nsplinesmat'=`nknotsmat'+`powermat'
	local nobs= `nobs'
	local dof=`DFR'	   
	//local nx: word count `xvars'
	if `nx2'>0{
		tempname b1 V1
		mat `b1'=`b'[1,1..`nx2']
		mat `V1'=`V'[1..`nx2',1..`nx2']
		mat colnames `b1'= `xvars2'
		mat colnames `V1'= `xvars2'
		mat rownames `V1'= `xvars2'
	
	}
	local labz
	local kk=1
	foreach z of local zvars{
		local uk: word `kk' of `uvars'
		local nspline=`nsplinesmat'[1,`kk']
		forv j=1/`nspline'{
			local labz `labz' `z'*S(`uk')_`j'
		
		}
		local kk=`kk'+1
	
	}
	
	mat colnames `b'= `xvars2' `labz'
	mat colnames `V'= `xvars2' `labz'
	mat rownames `V'= `xvars2' `labz'
	//mat list `b'

	//qui keep if `dataflag'==1
	//qui xtset `id' `time'
	sort `id' `time'
	noi disp in green""
	noi disp in green "Partially linear functional-coefficient panel data model."
	noi disp in green "Fixed-effect series semiparametric estimation."

    //noi di in green "" 
if (`nx2'>0){	
		ereturn local vcetype `vcetype'

			if `brep'>0{
				ereturn post  `b1' `V1', dep(`depvar') obs(`nobs') esample(`flag') buildfvinfo
			}
			else{
				ereturn post  `b1', dep(`depvar') obs(`nobs') esample(`flag') buildfvinfo
			}
		
		disp "Estimation results: linear part" _c

		noi di in green "{col 48} Number of obs        =" in yellow %8.0f `nobs'
		noi di in green "{col 48} Within R-squared     =" in yellow %8.4f `r2'
		noi di in green "{col 48} Adj Within R-squared =" in yellow %8.4f  `r2_a'
		noi di in green "{col 48} Root MSE             =" in yellow %8.4f `RMSE'		
		ereturn display, level(`level') 		
	}
	else{
		ereturn post, dep(`depvar') obs(`nobs') esample(`touse') buildfvinfo
	
	}

		noi disp in green "Estimated functional coefficient(s) are saved in the data."
		
 
		mata: U=st_data(.,"`uvars'")
        forv j=1/`nu'{
		    local nk=`nknotsmat'[1,`j']
			local np=`powermat'[1,`j']
			mata: s=bsknots(`qflag',U[.,`j'],`nk',`np')
			mata: st_local("k`j'",concat(strofreal(s)," "))
			ereturn local k`j' `k`j''
		}
 
		ereturn mat bs=`b'
	    ereturn mat Vs=`V'
		ereturn mat nknots = `nknotsmat'
		ereturn mat power  = `powermat'
		ereturn scalar df_m   = `DFM'
		ereturn scalar df_r   = `DFR'
		ereturn scalar r2     = `r2'
		ereturn scalar rmse   = `RMSE'
		ereturn scalar mss    = `MSS'
		ereturn scalar rss    = `RSS'
		ereturn scalar r2_a   = `r2_a'
		ereturn scalar N      = `nobs'
		ereturn local cmd "xtplfc" 
		ereturn local title "Partially linear functional-coefficient panel data model"
		ereturn local model "Fixed-effect Series Semiparametric Estimation" 
		ereturn local depvar "`depvar'" 
		ereturn local estfun `estfun'
		ereturn local basis "B-spline"
		ereturn local vcetype `vcetype'


    if "`predict'"!=""{
		gettoken v pname: pname
		qui gen `v'=`yhat'
		label var `v' "Fitted values of the dependent variable"
		if "`noai'"==""{
		   qui gen `pname'=`ui'
		   label var `pname' "Estimates of the fixed effects"				
		
		}
   }  


end




****************************************************************************

/*
cap program drop splitsample

program define splitsample

version 14

syntax [if] [in], gen(string) [n(integer 2)]

marksample touse

*count if `touse'
confirm new var `gen'

tempvar rn
qui gen `rn'=uniform() if `touse'

*su `rn'

local n1=100/`n'

local n2=100*(`n'-1)/`n'

_pctile `rn', p(`n1'(`n1')`n2' )

		qui gen int `gen'=.


		forv j=1/`=`n'-1' {

		   qui replace `gen'=`j' if `rn'<=r(r`j') & missing(`gen') & `touse'
		}

		qui replace `gen'=`n' if missing(`gen') & `touse'

end
*/
* 2020-6-15
cap program drop splitsample

program define splitsample

version 14

syntax [if] [in], gen(string) [n(integer 2)]

marksample touse

*count if `touse'
confirm new var `gen'

tempvar rn
qui gen `rn'=uniform() if `touse'

qui xtile `gen' = `rn' if `touse', n(`n')

end
