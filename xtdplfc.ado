*! version 7.0
* 2019-8-26
* change cvtwo option to tenfoldcv
* bugs fix
*! version 6.0.5 
* 2019-01-25
*12Dec2018
*8Dec2018
*9Nov2018
*By Kerry Du
cap program drop xtdplfc
program define xtdplfc, eclass prop(xt)

	version 14
if replay()& "`e(cmd)'"=="xtdplfc" {
	ereturn  display
	exit
}

syntax varlist, GENerate(string)  Uvars(varlist) [Zvars(varlist) endox(varlist) ENDOZflag(numlist integer >0) ///
				   ivx(varlist)  lags(integer 1) lagyinz(numlist integer >0) ONLYivxz ///
				   Ivtype(numlist integer >=0  )  NKnots(numlist integer >=2) power(numlist integer >0)  ///
				   grid(string) pctile(integer 0) brep(integer 200) QUANtile ///
				   level(integer 95)  WILD  predict(string) ivz(string) TE ///
				   MINNKnots(numlist integer >=2) MAXNKnots(numlist integer >=2) NODOTS FAST  TENfoldcv]	

				   
				   
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
/*
	capture which lrcov
	local rc = _rc
	if (`rc') {
		display as error "{bf:xtplfc} requires {bf:lrcov}"
		exit 198
		}				   
*/
*************************************************************				   
				   
				   
				   
***************************************************************				   
if"`fast'"!=""{
xtdplfc_fast `varlist', generate(`generate')  uvars(`uvars') zvars(`zvars') endox(`endox') ///
				   endozflag(`endozflag') ivx(`ivx') lags(`lags') lagyinz(`lagyinz') `onlyivxz' ///
				   ivtype(`ivtype')  nknots(`nknots') power(`power')  ///
				   grid(`grid') pctile(`pctile') brep(`brep') `quantile' ///
				   level(`level')  `wild'  predict(`predict') `te' ///
				   minnknots(`minnknots') maxnknots(`maxnknots') `nodots' ivz(`ivz')  `tenfoldcv'	

		exit
}				   
				   
***************************************************************				   
				   

if `lags'<1{
		display as error "No lag dependent variables included."
		display as error "Use xtplfc or ivxtplfc instead."
		exit	
}

**************************************************************
	
	local nz0: word count `zvars'
	local nu: word count `uvars'	
	local comv: list zvars & uvars
	if "`comv'"!=""{
		di as error "`comv' included in both covariates and functions"
		di as error "Check zvars() and uvars()"
		exit 198
	}


	gettoken depvar indepvars : varlist

****************************************************************
****Analyze the time structure	
 	allocy `depvar',lags(`lags') lagyinz(`lagyinz')
    local Dy1inx   `r(Dy1inx)'
    local laby1inx `r(laby1inx)'
    local Doyinx 	`r(Doyinx)'
    local oyinx `r(oyinx)'
    local laboyinx  `r(laboyinx)'
    local laby1inz `r(laby1inz)'
    local oyinz `r(oyinz)'
    local laboyinz `r(laboyinz)' 
    local yinz `r(yinz)'
    local yinx `r(yinx)'

	local nyinz: word count `yinz'
	local nz=`nz0'+`nyinz'
	
	if `nz'!=`nu'{
		di as error "# of variables in uvars() not equal to that specified by zvars() and lagyinz() "
		exit
	}
**********************************************************************
***Add lagyinz in the front of zvars
    
	local zvars `laby1inz' `laboyinz' `zvars'

*************************************************************************

if `"`ivz'"'!=""{
	checkivz `ivz'
	local ivzlist `r(ivlist)'
	local ivuflag `r(uflag)'
	local ivztypelist `r(typelist)'	
} 


foreach v of local ivuflag{
	if `v'>`nz'{
		di as error "Error in ivz()"
		di as error "# specified by uflag() > # of variables in uvars()"
		exit
	}
}

************************************************************************
if "`endozflag'"!=""{
	numlist "`endozflag'",sort
	local endozflag `r(numlist)'
	local endozflag: list uniq endozflag	
}


foreach v of local endozflag{
	if `v'>`nz'{
		di as error "Error in endozflag()/lagyinz()"
		di as error "# specified by endozflag()/lagyinz() > # of variables in uvars()"
		exit
	}
	if `v'<`nyinz'{
		di as error "The first `nyinz' variables in uvars() interact with the lags of depvar specified by lagyinz()"
		di as error "The # in endozflag() should be > `nyinz'"
		exit
	}

}


if "`lagyinz'"!=""{
	local endozflag2 1 `endozflag'	// adding the first lag of depvar
}
else{
	local endozflag2  `endozflag'
}


//local endozfalg2: list uniq endozflag2

*********************************************************************
**********************************************************************

//index the uvars
qui numlist "1/`nu'"
local allflags=r(numlist)

local exgozflag: list allflags - endozflag2

local nendoz: word count `endozflag2'
if "`ivtype'"==""{
	forv i=1/`nendoz'{
		local ivtype `ivtype' 1   // by default
	}
	
}
else{
	local nivtype: word count `ivtype'
	if `nivtype'!=`nendoz'{
		di as error "# of integer in ivtype()!= # of variables specified by endozflag() and lagyiz()"
		exit
	}
}

//mata: ivtype=tokens(st_local("ivtype"))


***********************************************************************

	if `"`endox'"'!="" & "`ivx'"==""{
			disp as error "Instruments should be provided for variables in endox()." 
			disp as error "Using ivx()."
			exit 198
	}


*************************************************************************


*************************************************************************

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



********************************************************************************	
* powermat nknotsmat maxnknotsmat minnknotsmat nknots4cv nsplinesmat nknots4cv


	tempname powermat nknotsmat maxnknotsmat minnknotsmat mxmi nsplinesmat nknots4cv
	if "`power'"==""{
		//disp "By default, the number of power (degree) is set to 3 for all functions."
		mat `powermat'=J(1,`nu',3)
		mata: powermat=st_matrix("`powermat'")
	}
	else{
		local npow: word count `power'
		if `npow'!=`nu' {
			di as error "# of integers in power() != # of variables in uvars()."
			exit			
		}
		mata: powermat=strtoreal(tokens(st_local("power")))
		mata: st_matrix("`powermat'",powermat)	
		
	}	

	
	if "`maxnknots'"=="" & "`nknots'"!="" {

		local nk0: word count `nknots'
		if `nk0'!=`nu'{
			di as error "# of integers in knots() != # of variables in uvars()."
			exit
		}

		mata: nknotsmat=strtoreal(tokens(st_local("nknots")))
		mata: st_matrix("`nknotsmat'",nknotsmat)

	}
	
	
	if "`maxnknots'"=="" & "`nknots'"=="" {
		mat `nknotsmat'=J(1,`nz',2)
		mata: nknotsmat=st_matrix("`nknotsmat'")
		//disp "By default, # of knots is set to 2 for all functions."
	}
	
	if "`maxnknots'"!=""{
		local nmaxnk: word count `maxnknots'
		if `nmaxnk'!=`nu'{
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
			di as error "# of integers in minnknots != # of variables in uvars()."
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
	
////////////////////////////////////////////////////////////	
	if "`maxnknots'"!=""{
		//mata: nsplines=maxnknotsmat+powermat
		//mata: st_matrix("`nsplinesmat'",nsplines)
		mat   `nsplinesmat'=`maxnknotsmat'+`powermat'
		mata: nknots4cv=combmat(maxnknotsmat,minnknotsmat)
		mata: st_matrix("`nknots4cv'",nknots4cv)
	}
	else{
		mat `nsplinesmat'=`nknotsmat'+`powermat'

	}

*********************************************************************
*deal with collinearity
*********************************************************************
	local comv: list indepvars | endox
	local exgox: list comv - endox

	local comv: list comv & zvars

	if !(`"`comv'"'==""){
		disp as red "`comv' included in both linear and nonlinear variable lists." 
		exit 198
	}


	
	if `"`endox'"'!=""&`"`exgox'"'!=""{
	 		getcinfo2list, alist(`endox') blist(`exgox')
	 		local endox `r(cname1)'
	 		local exgox `r(cname2)'
	 	}
	 else if `"`endox'"'==""&`"`exgox'"'!=""{
	 		getcinfoalist, alist(`exgox') 
	 		local exgox `r(cnames)'
	 	}
	 else if `"`endox'"'!=""&`"`exgox'"'==""{
	  		getcinfoalist, alist(`endox') 
	 		local endox `r(cnames)'	
	 	}


 


	if "`ivx'"!=""{
	 		getcinfoalist, alist(`ivx') 
	 		local ivx `r(cnames)'
		}


	local xvars `endox' `exgox'

/////////////////////////////////////////////////////////////////////
* construct IVs

	if "`onlyivxz'"==""{
		local ivtypecopy `ivtype'
		if "`laby1inz'"!=""{
			local ivzlist L2.`depvar' `ivzlist'
			local ivuflag  1 `ivuflag'
			gettoken ivtemp ivtypecopy: ivtypecopy
			local ivztypelist `ivtemp' `ivztypelist' 
		}
		else{
			local ivx L2.`depvar' `ivx'
		}

		foreach v of local endox{
			local ivx `ivx' L2.`v'
		}

		foreach v of local endozflag{
			local zv: word `v' of `zvars'
			local ivzlist `ivzlist' L2.`zv'
		}

		local ivuflag `ivuflag' `endozflag'
		local ivztypelist `ivztypelist' `ivtypecopy'

	}

	local ivlist `ivx'
	local ivzlistcopy `ivzlist'
	local ivztypelistcopy `ivztypelist'
	foreach uk of local ivuflag{
		gettoken zk ivzlistcopy: ivzlistcopy
		gettoken ivk ivztypelistcopy: ivztypelistcopy
		local uk0: word `uk' of `uvars'
		local ivname `zk'*L`ivk'.[Splines(`uk0')]
		local ivname=subinstr(`"`ivname'"',"L0.[","L.[",.)
		local ivlist `ivlist' `ivname'
	}		
*********************************************************************	
  local vlabel=cond("`grid'"=="","observed","grid")
  local qflag=cond("`quantile'"=="",1,0)	

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
		local k=`k'+1	
	}
  }
  
  

 /*   
	local zvars2name	
	
	forv uk=1/`nu'{
		
		local zk: word `uk' of `zvars'
		local uk0: word `uk' of `uvars'
		if `"`zk'"'==`"L.`depvar'"'{
			local zkname L1_`depvar'
		}
		else{
			local zkname=strtoname(`"`zk'"',1)	
		}
		local zvars2name `zvars2name' `zkname'
		confirm new variable `generate'_`zkname'
		confirm new variable `generate'_`zkname'_sd
		qui gen `generate'_`zkname'=.
		qui gen `generate'_`zkname'_sd=.	
		label var `generate'_`zkname' `"Estimated g(`uk0') at the `vlabel' points w.r.t `zk'*g(`uk0')"'	
		label var  `generate'_`zkname'_sd  `"`vcetype' S.E. of g(`uk0') w.r.t `zk'*g(`uk0')"'
	
	}	
*/

	forv uk=1/`nu'{
		
		local zk: word `uk' of `zvars'
		local uk0: word `uk' of `uvars'
		/*
		if `"`zk'"'==`"L.`depvar'"'{
			local zkname L1_`depvar'
		}
		else{
			local zkname=strtoname(`"`zk'"',1)	
		}
		local zvars2name `zvars2name' `zkname'
		*/
		confirm new variable `generate'_`uk'
		confirm new variable `generate'_`uk'_sd
		qui gen `generate'_`uk'=.
		qui gen `generate'_`uk'_sd=.	
		label var `generate'_`uk' `"Estimated g(`uk0') at the `vlabel' points w.r.t `zk'*g(`uk0')"'	
		label var  `generate'_`uk'_sd  `"`vcetype' S.E. of g(`uk0') w.r.t `zk'*g(`uk0')"'
	
	}	
  
  
*********************************************************************	
/*
	_xt, trequired 
	local id=r(ivar)
	local time=r(tvar)
	
/*	
//fill the gap
	tempvar dataflag
	qui gen `dataflag'=1
	qui tsfill,full

//continuous id and time
	tempvar cid ctime
	qui egen `cid'=group(`id')
	qui egen `ctime'=group(`time')
	qui tab `cid', nofreq
	local nid=r(r)
	qui tab `ctime',nofreq
	local ntime=r(r)
	local te=cond("`te'"=="",0,`ntime'-1)
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
		noi display as error "{bf:xtplfc} needs strongly balanced panel, use {bf:tsfill, full} to rectangularize your data!."
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
	
	if `te'>0{
		mata: tematrix=gendtime(`nid',`te')
		mata: _sttempvar("tematrix",tematrix)
	}
		
*********************************************************************
**********define tempvar for regression
  forv k=1/`nz'{
	local nspline=`nsplinesmat'[1,`k']
	forv j=1/`nspline'{
	   tempvar  sp`k'_`j'
	}
/*	
  	local uk: word `k' of `uvars'
  	su `uk',meanonly
 	local min_`k'=r(min)
	local max_`k'=r(max) 	
*/	
	
  }
  
/*
  forv j=1/`nu'{
  	local uj: word `j' of `uvars'
  	su `uj',meanonly
 	local min_`j'=r(min)
	local max_`j'=r(max) 	
  }
*/	

	if "`nodots'"!=""{
		local nodots qui
	}
  
***************************************************************************
	local Dexgox
	
	foreach v of local exgox {
		local Dexgox `Dexgox' D.`v'
	}
	
	local Dendox
	
	foreach v of local endox {
		local Dendox `Dendox' D.`v'
	}
	
	local Dexgox2 `Dexgox' `tematrix'

//qui xtset
//////////////////////////////////////////////////////////////////////////
	if "`maxnknots'"!=""{
			qui putmata U=(`uvars'),replace	
			tempvar ehat_sq  yhat 
			tempname mse 		
			di
			disp as green "Using Cross-Validation to determine the optimal number of knots..."
			local nnknots=rowsof(`nknots4cv')
			mat `mse'=J(`nnknots',1,.)
			qui gen `ehat_sq'=.

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
			
			
		//qui xtset
		forvalues kk=1/`nnknots'{ 
		

			 **generating spline bases
			 local splines	
			 
			 forv ii=1/`nu' {
				
				local nknot`ii'=`nknots4cv'[`kk',`ii']
				local power_`ii'=`powermat'[1,`ii']
				mata: _gknotlist("klist`ii'", `qflag', U[.,`ii'], `nknot`ii'',`power_`ii'')
				local sp`ii'
				local nspline=`nknot`ii''+`power_`ii''
				

				forv zi=1/`nspline'{
					local sp`ii' `sp`ii'' `sp`ii'_`zi''
					}
				local uii: word `ii' of `uvars'
				qui bspline `sp`ii'', xvar(`uii') knots(`klist`ii'') power(`power_`ii'') noexknot	
				local splines `splines' `sp`ii''
			  }	  

		        qui xtset `id' `time'
				forv uk=1/`nu'{
					local HZ`uk'
					local lHZ`uk'
					local zk: word `uk' of `zvars'
					local uk0: word `uk' of `uvars'
					foreach uki of local sp`uk'{
						tempvar z`uki'
						qui gen `z`uki''=`zk'*`uki'
						local HZ`uk' `HZ`uk'' D.`z`uki''
						local lHZ`uk' `lHZ`uk'' `z`uki''
					}
				
				}
				
					
					local DEXGZ
					local EXGZ
					foreach uk of local exgozflag {	
						local DEXGZ `DEXGZ' `HZ`uk''
						local EXGZ `EXGZ' `lHZ`uk''
					}
					
					
					local DENDZ
					local ENDZ
					foreach uk of local endozflag2 {
						local DENDZ `DENDZ' `HZ`uk''
						local ENDZ `ENDZ' `lHZ`uk''					

					}
					

					local instz
					local IVZs `ivztypelist'
					local ivzlistcopy `ivzlist'
					foreach uk of local ivuflag{
						gettoken zk ivzlistcopy: ivzlistcopy
						local zn=subinstr(`"`zk'"',".","",.)
						gettoken zktype IVZs: IVZs
						local uk0: word `uk' of `uvars'				
						foreach uki of local sp`uk' {
							tempvar I`zn'`uki'
							qui gen `I`zn'`uki''=`zk'*L`zktype'.`uki'
							local instz `instz' `I`zn'`uki''	
		
						}					

					}



				  if "`tenfoldcv'"==""{

					   forvalues jt=1/`nid' {



						  qui ivregress 2sls D.`depvar' (`Dy1inx' `Dendox' `DENDZ' = `ivx' `instz' )  ///
									   `DEXGZ' `Doyinx' `Dexgox2' if `cid'~=`jt', noconstant 


							qui predict `yhat' if `cid'==`jt',xb
							qui replace `ehat_sq'=(D.`depvar'-`yhat')^2 if `cid'==`jt'
							qui cap drop `yhat'
						}
					}
					else{
					   qui xtset `id' `time'
					   forvalues jt=1/10 {



						  qui ivregress 2sls D.`depvar' (`Dy1inx' `Dendox' `DENDZ' = `ivx' `instz' )  ///
									   `DEXGZ' `Doyinx' `Dexgox2' if `splitdummy'!=`jt', noconstant 


							qui predict `yhat' if `splitdummy'==`jt',xb
							qui replace `ehat_sq'=(D.`depvar'-`yhat')^2 if `splitdummy'==`jt'
							qui cap drop `yhat'
						}					
					}

					qui su `ehat_sq' if ~missing(`ehat_sq')
					mat `mse'[`kk',1]=r(sum)/r(N)
					cap drop `ENDZ'
					cap drop `EXGZ'
					cap drop `splines'
					cap drop `instz'
					`nodots' displaydot, dot(`kk')

		   }


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
			//di 

	}	



*************************************************************************	

  qui putmata U=(`uvars'),replace	
  mat `nsplinesmat'	=`nknotsmat'+`powermat'
  mata: nsplinesmat=st_matrix("`nsplinesmat'")
  mata: bzindex=0,cumsum(nsplinesmat)
  tempname bzindex	
  mata: st_matrix("`bzindex'",bzindex)

	local N=_N
	local splines

		 local splines	
	
		 forv ii=1/`nu' {
			
			local nknot`ii'=`nknotsmat'[1,`ii']
			local power_`ii'=`powermat'[1,`ii']
			mata: _gknotlist("klist`ii'", `qflag', U[.,`ii'], `nknot`ii'',`power_`ii'')
			local sp`ii'
			local nspline=`nknot`ii''+`power_`ii''
			

			forv zi=1/`nspline'{
				local sp`ii' `sp`ii'' `sp`ii'_`zi''
				}
			local uii: word `ii' of `uvars'
			qui bspline `sp`ii'', xvar(`uii') knots(`klist`ii'') power(`power_`ii'') noexknot		
			local splines `splines' `sp`ii''
			local knotlist_`ii' `r(knots)'
		  }	  

			sort `id' `time'
			qui putmata HH=(`splines'),replace		
			
			
			//qui xtset
			qui xtset `id' `time' 
			local Zalls
			local lablHZ
			forv uk=1/`nu'{
				local HZ`uk'
				local lHZ`uk'
				local zk: word `uk' of `zvars'
				local uk0: word `uk' of `uvars'
				local j=1
				foreach uki of local sp`uk'{
					tempvar z`uki'
					qui gen `z`uki''=`zk'*`uki'
					local HZ`uk' `HZ`uk'' D.`z`uki''
					local lHZ`uk' `lHZ`uk'' `z`uki''
					local lablHZ `lablHZ' `zk'*S(`uk0')_`j'
					local j=`j'+1
				}
				local Zalls `Zalls' `HZ`uk''
			
			}
			
				
				local DEXGZ
				local EXGZ
				foreach uk of local exgozflag {	
					local DEXGZ `DEXGZ' `HZ`uk''
					local EXGZ `EXGZ' `lHZ`uk''
				}
				
				
				local DENDZ
				local ENDZ
				foreach uk of local endozflag2 {
					local DENDZ `DENDZ' `HZ`uk''
					local ENDZ `ENDZ' `lHZ`uk''					

				}
				

				local instz
				local IVZs `ivztypelist'
				local ivzlistcopy `ivzlist'
				foreach uk of local ivuflag{
					gettoken zk ivzlistcopy: ivzlistcopy
					local zn=subinstr(`"`zk'"',".","",.)
					gettoken zktype IVZs: IVZs
					local uk0: word `uk' of `uvars'				
					foreach uki of local sp`uk' {
						tempvar I`zn'`uki'
						qui gen `I`zn'`uki''=`zk'*L`zktype'.`uki'
						local instz `instz' `I`zn'`uki''
						
	
					}					

				}


	qui ivregress 2sls D.`depvar' (`Dy1inx' `Dendox' `DENDZ' = `ivx' `instz' )  ///
	 					   `DEXGZ' `Doyinx' `Dexgox2', noconstant 			


   	local Xalls `Dy1inx' `Doyinx' `Dendox' `Dexgox2'
	local nbxs: word count `Xalls'
	local allreg `Dy1inx' `Dendox' `DENDZ' `DEXGZ' `Doyinx' `Dexgox2'
	//su `allreg'

	local nobs   =   e(N) 
	
	mata: allreg=tokens(st_local("allreg"))
	mata: Xalls=tokens(st_local("Xalls"))
	mata: Zalls=tokens(st_local("Zalls"))

	mata: xpos=posof(Xalls,allreg)
	mata: zpos=posof(Zalls,allreg)
	
	
  tempvar flag xb ehat aui Dehat yhat res2 Dxb
    //tempvar flag xb ehat aui Dehat res2 ystar Dxb
   qui gen `flag'= e(sample)
   qui predict `Dehat', residuals
   qui predict `Dxb', xb
   local dfm: word count `allreg'
   local dof=`nobs'-`dfm'
   local dfr=`dof'
   tempvar e_sq yd_sq
   qui gen `e_sq'=`Dehat'^2
   su `e_sq', meanonly
   local rss=r(sum)
   local rmse=(`rss'/`nobs')^0.5
   su D.`depvar', meanonly
   local dymean=r(mean)
   qui gen `yd_sq'=(`Dxb'-r(mean))^2
   su `yd_sq',meanonly
   local mss=r(sum)
   qui replace `yd_sq'=(D.`depvar'-`dymean')^2
   su `yd_sq',meanonly
   local r2=1-`rss'/r(sum)
   local r2_a=1-(`rss'/`dfr')/(r(sum)/(`nobs'-1))
  //sum `Dehat'
   tempname b0 b1 b2 V1 V2 B V
    mat `b0'=e(b)
	mata: b0=st_matrix("`b0'")
	
	mata: b2=b0[1,zpos]
	mata: st_matrix("`b2'",b2)
	mat `B'=`b2'
	
	sort `id' `time'
	//qui xtset
	mata: XXX=st_data(.,"`laby1inx' `endox' `ENDZ' `EXGZ' `laboyinx' `exgox'")
	mata: b00=b0[1,1..(cols(b0)-`te')]
	mata: yhat=XXX*b00'
	if `te'>0{
		mata: timeeffects=b0[1,(cols(b0)-`te'+1)..cols(b0)],-sum(b0[1,(cols(b0)-`te'+1)..cols(b0)])
		mata: timeeffects=J(`nid',1,timeeffects')
		mata: yhat=yhat+timeeffects
		tempvar timeeffects
		qui getmata `timeeffects'=timeeffects,replace
		local sumlagys `timeeffects'
	}
	else{
        local sumlagys	0	
	}


	qui getmata `yhat'=yhat, replace
	qui gen `ehat'=`depvar'-`yhat'
	qui bys `cid' (`ctime'): egen `aui'=mean(`ehat') 
	qui replace `ehat'=`ehat'-`aui'
	qui replace `yhat'=`yhat'+`aui' if `ctime'>`lags'
	
    if "`predict'"!=""{
		gettoken v pname: pname
		qui gen `v'=`yhat'
		label var `v' "Fitted values of the dependent variable"
		if "`noai'"==""{
		   qui gen `pname'=`aui'
		   label var `pname' "Estimates of the fixed effects"				
		
		}
   }  	
	
	//qui replace `depvar'=`yhat'
	
	if "`Xalls'"!=""{
		mata: b1=b0[1,xpos]
		mata:st_matrix("`b1'",b1)
		mat `B'=(`b1',`b2')
	}




	//local sumlagys `sumlagys'+`aui'
	local q=1
	foreach z of local xvars{

		local _b`q'=_b[D.`z']
		local sumlagys `sumlagys'+`_b`q''*`z'
		local q=`q'+1
	
	}
	

	
	if "`Dy1inx'"!=""{
	
		local _bl1y=_b[D.L.`depvar']
		local sumlagys `sumlagys'+`_bl1y'*L.`yhat'
	}
	

	foreach j of local oyinx{
	
		local _bl`j'=_b[D.L`j'.`depvar']
		local sumlagys `sumlagys'+`_bl`j''*L`j'.`yhat'
	}
	
	
	mata: bz=b0[1,zpos]

/*
	local estfun
	mata: gfun=estgfun(bz',HH,bzindex)
	local k=1
	local ZZname `zvars2name'
	foreach j of local zvars {
		gettoken zj ZZname: ZZname
		mata: temvec=gfun[,`k']
		qui getmata `generate'_`zj'=temvec,replace
		local estfun `estfun' `generate'_`zj'
		local j=subinstr(`"`j'"',`".`depvar'"',".`yhat'",.)
		local sumlagys `sumlagys'+`generate'_`zj'*`j'
		//disp "`generate'_`zj'*`j'"
		local k=`k'+1
	
	}
*/

	local estfun
	mata: gfun=estgfun(bz',HH,bzindex)
	local k=1
	foreach j of local zvars {
		mata: temvec=gfun[,`k']
		qui getmata `generate'_`k'=temvec,replace
		local estfun `estfun' `generate'_`k'
		local j=subinstr(`"`j'"',`".`depvar'"',".`yhat'",.)
		local sumlagys `sumlagys'+`generate'_`k'*`j'
		//disp "`generate'_`zj'*`j'"
		local k=`k'+1
	}
	
	//disp "`sumlagys'"



  
//disp "`sumlagys'"

  **************************************

   if "`grid'"!=""{	
	 forv zk=1/`nz' {
		//local uk: word `zk' of `uvars'
		//local udist=(`max_`zk''-`min_`zk'')/`nknots_`zk''
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
	  mata: gfun0=estgfun(bz',HH,bzindex)
   }
   else{
	  mata:gfun0=gfun
   }
   
   
   
	
*************************************************************	


*************************************************************

if `brep'>0 {
	di
	disp "Computing the bootstrap standard errors..."
	
	local bootcmd bs_resid3, id(`cid') res(`ehat') replace(`res2')
	if "`wild'"!=""{
		local bootcmd bs_wild3, res(`ehat') replace(`res2')
	}
	
	
	mata: bsmat=J(`brep',length(b0),.)


	forvalues z = 1/`nz'{
		mata: gf_`z'=J(`N',`brep',.)
	}
	
	sort `id' `time'
	//qui xtset `cid' `ctime'
	//qui xtset
	qui gen `res2'=.
	
	local ivx2 `ivx'
	local Dy1inx2
	if "`Dy1inx'"!=""{
		local Dy1inx2 D.L.`yhat'
		gettoken ivfirstlag ivx2: ivx2
		local ivx2 L2.`yhat' `ivx2'
	}
	
	//disp "`ivx2'"
	
	local Doyinx2
	foreach lgt of local oyinx{
		local Doyinx2 `Doyinx2' D.L`lgt'.`yhat'
	}
	
	
	tempname bb 
	forv tt=1/`brep' {
			qui replace `yhat'=`depvar'
			`nodots' displaydot, dot(`tt')
			qui `bootcmd'
			qui replace `yhat'=`sumlagys'+`aui'+`res2' if `ctime'>`lags'

			local j=0
			foreach	t of local yinz{
				local j=`j'+1
				local uk: word `j' of `uvars'
				foreach uki of local sp`j'{
					qui replace `z`uki''=L`t'.`yhat'*`uki'

				}

			}
				

			if "`onlyivxz'"==""&"`laby1inz'"!=""{
					local zktype: word 1 of `ivztypelist'
					foreach j of local sp1{
						qui replace `IL2`depvar'`j''=L2.`yhat'*L`zktype'.`j'
						
					 }

			        
			     }


	//su `yhat'	
		qui ivregress 2sls D.`yhat' (`Dy1inx2' `Dendox' D.(`ENDZ') = `ivx2' `instz' )  ///
	 					                   D.(`EXGZ') `Doyinx2' `Dexgox2', noconstant 
		//su `depvar'
		mat `bb'=e(b)
		mata: bb=st_matrix("`bb'")
		mata: bb2=bb[1,zpos]
		if `nbxs'>0	{
			mata: bsmat[`tt',.]=bb[1,xpos],bb2
		}
		else{
			mata: bsmat[`tt',.]=bb2
		}

			
			forvalues z=1/`nz'{
				mata: bz=bb2[1,(1+bzindex[1,`z'])..bzindex[1,`z'+1]]
				mata: gf_`z'[.,`tt']=HH[.,(1+bzindex[1,`z'])..bzindex[1,`z'+1]]*bz'

			}
    }
	
	//qui replace `depvar'=`depvar0'

	mata: V=quadvariance(bsmat)
	mata: st_matrix("`V'",V)
	if "`Xalls'"!=""{
		mat `V1'=`V'[1..`nbxs',1..`nbxs']
	}
	
	sort `id' `time'
		* compute the s.e. for the smooth function
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
	if "`Xalls'"!=""{
		mat `V1'=`V'[1..`nbxs',1..`nbxs']
	}
}
/*

else{
	
		   tempname omega
		   qui sort `id' `time'
		   local WW `ivx' `instz' `DEXGZ' `Doyinx' `Dexgox2'
		   *su `WW'
		   
		   mata: W=st_data(.,"`WW'","`flag'")

		   mata: HZ=st_data(.,"`Zalls'","`flag'")

		   qui putmata cid=`cid' Dehat=`Dehat'  if `flag'==1, replace
		   if (`"`hacopt'"'==""& "`weightvar'"==""){

				mata: omega=lrvar(cid,W,Dehat)
				*mata: sum(omega)
		   }
		   else{
		   
			   if ("`weightvar'"!="") {
					qui replace `Dehat'=`Dehat'*`weightvar'
				}		   	
			   //qui xtlrcov `WW', id(`cid') wvar(`Dehat') `hacopt' touse(`flag')
			   qui xtlrcov `WW', id(`cid') time(`ctime') wvar(`Dehat') `hacopt' 
			   mat `omega'=r(omega)
			   mata: omega=st_matrix("`omega'")	
		   }


	    mata: omega=omega/rows(W)
	

		if "`Xalls'"!=""{

			sort `id' `time'
			mata: X=st_data(.,"`Xalls'","`flag'")
			mata: gsd=estcov(W,(X,HZ),omega,"`V'",`nbxs',bzindex,HH)
			mat `V1'=`V'[1..`nbxs',1..`nbxs']

			
		}
		else{
		
			mata: gsd=estcov(W,HZ,omega,"`V'",`nbxs',bzindex,HH)
		
		}
	
		
	  //// qui getmata `generate'_sd*=gsd,replace //////
	   local k=1

	   foreach v of local zvars2name{
	   		mata: tempgsd=gsd[.,`k']
			qui getmata `generate'_`v'_sd=tempgsd,replace
			local k=`k'+1
	   
	   }
	   	   


}
*/

	local k=1
	foreach j of local zvars{
		mata: temvec=gfun0[,`k']
		qui getmata `generate'_`k'=temvec,replace
		local k=`k'+1
	
	}



//	qui keep if `dataflag'==1
//	qui xtset `id' `time'
	sort `id' `time'
	mat colnames `B' = `laby1inx' `laboyinx' `endox' `exgox' `labtematrix' `lablHZ' 
	mat colnames `V' = `laby1inx' `laboyinx' `endox' `exgox' `labtematrix' `lablHZ' 
	mat rownames `V' = `laby1inx' `laboyinx' `endox' `exgox' `labtematrix' `lablHZ'  	

	noi disp in green""
	noi disp in green "Partially linear functional-coefficient dynamic panel data model."
	noi disp in green "Fixed-effect sieve 2SLS estimation."

	if ("`Xalls'"!=""){
			mat colnames `b1' = `laby1inx' `laboyinx' `endox' `exgox' `labtematrix'
			mat colnames `V1' = `laby1inx' `laboyinx' `endox' `exgox' `labtematrix'
			mat rownames `V1' = `laby1inx' `laboyinx' `endox' `exgox' `labtematrix'
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
			noi di in green "{col 48} Root MSE             =" in yellow %8.4f `rmse'		
			ereturn display, level(`level') 		
		}
		else{
			ereturn post, dep(`depvar') obs(`nobs') esample(`flag') buildfvinfo
		}

		noi disp in green "Estimated functional coefficient(s) are saved in the data."
		noi disp in green "Instruments for differenced equation: `ivlist'"

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
		ereturn local cmd "xtdplfc" 
		ereturn local title "Partially linear functional-coefficient dynamic panel data model"
		ereturn local model "Fixed-effect Sieve 2SLS Estimation" 
		ereturn local depvar "`depvar'" 
		ereturn local estfun `estfun'
		ereturn local basis "B-spline"
		ereturn local vcetype `vcetype'
		ereturn local ivlist `ivlist'
/*		
		mata: U=st_data(.,"`uvars'")
        forv j=1/`nu'{
		    local nk=`nknotsmat'[1,`j']
			local np=`powermat'[1,`j']
			mata: s=bsknots(U[.,`j'],`nk',`np')
			mata: st_local("k`j'",concat(strofreal(s)," "))
			ereturn local k`j' `k`j''
		}
	
	
    if "`predict'"!=""{
		gettoken v pname: pname
		qui gen `v'=`yhat'
		label var `v' "Fitted values of the dependent variable"
		if "`noai'"==""{
		   qui gen `pname'=`ui'
		   label var `pname' "Estimates of the fixed effects"				
		
		}
     }  			
*/	
		forv j=1/`nu'{
			ereturn local k`j' `knotlist_`j''
		}	
	

end



********************************************************************
cap program drop checkivz
program define checkivz, rclass
	version 14
	syntax varlist, Uflag(numlist) [ Ivtype(numlist integer >=0) ]
	
	
	numlist "`uflag'"
	local uflag `r(numlist)'

	local niv: word count `varlist'
	local nf: word count `uflag'
	

	if `niv'!=`nf'{
		di as error "Error in ivz()"
		di as error "# of instruments != # of integers specified in uflag()"
		exit
	}


	if "`ivtype'"==""{
		local ivtype 1
		local typelist
		foreach v of local varlist{
			local typelist `typelist' `ivtype'
		}		
	}
	else{
		local typelist `ivtype'
		local ntype: word count `ivtype'
		if `niv'!=`ntype'{
			di as error "Error in ivz()"
			di as error "# of instruments != # of integers specified in ivtype()"
			exit
		}		
	}

	return local typelist `typelist'
	return local ivlist `varlist'
	return local uflag `uflag'

end


***************************************************

cap program drop  bs_resid3
program define bs_resid3
        version 14
        syntax, id(varname numeric)  RESidual(varname numeric) replace(varname)

        qui putmata id=`id',replace
        mata: mm_panels(id,idinfo=.)
       // mata: length(idinfo)
        mata: idx=mm_sample(length(idinfo),.,idinfo)
        tempvar idx
        qui getmata `idx'=idx, replace 
        qui replace `replace'=`residual'[`idx']

end     



cap program drop  bs_wild3
program define bs_wild3
        version 14
        syntax, RESidual(varname numeric) replace(varname)

		tempvar randu ys flag 
		
		qui gen `randu'=uniform()
		qui gen `flag'=(1-sqrt(5))/2
		qui replace `flag'=(sqrt(5)+1)/2 if `randu'>((1+sqrt(5))/(2*sqrt(5)))
		
		qui replace `replace'=`residual'*`flag'

end 

   
****************************************************
****************************************************



cap program drop allocy
program define allocy, rclass
version 14
syntax varname(ts fv),lags(numlist integer max=1 >=1) [lagyinz(numlist integer)]

local depvar `varlist'
local lag2ys
forv j=2/`lags'{
	local lag2ys `lag2ys' `j'
}


local yinx `lag2ys'
local oyinx `lag2ys'
if "`lagyinz'"!=""{
	numlist "`lagyinz'", sort
	local yinz "`r(numlist)'"
	local yinz: list uniq yinz  
	local oyinx: list lag2ys - yinz 
	
}

local firstlag 1
local alllags `firstlag' `lag2ys'
if "`lagyinz'"!=""{
	*local alllags `firstlag' `lag2ys'
	local ifinlag: list yinz - alllags
	if "`ifinlag'"!="" {
		disp as error "{bf:lagyinz(numlist)}  has elements outside of allowed range specified by {bf:lags(#)}."
		exit
	}
}


local Dy1inx
local laby1inx
local laby1inz
local y1inz: list firstlag & yinz
if "`y1inz'"=="" {
	local Dy1inx D.L.`depvar'
	local laby1inx L.`depvar'
}
else{
	local laby1inz L.`depvar'
}



local Doyinx
local laboyinx
foreach j of local oyinx{
	
	local Doyinx    `Doyinx'   D.L`j'.`depvar'
	local laboyinx  `laboyinx' L`j'.`depvar'
}



local oyinz: list yinz - firstlag

local laboyinz
foreach j of local oyinz{
	local laboyinz `laboyinz'  L`j'.`depvar'
}

local yinx: list alllags - yinz

return local Dy1inx   `Dy1inx'
return local laby1inx `laby1inx'
return local oyinx `oyinx'
return local Doyinx `Doyinx'
return local laboyinx  `laboyinx'
return local laby1inz `laby1inz'
return local laboyinz `laboyinz'
return local oyinz `oyinz'
return local yinz `yinz'
return local yinx  `yinx'
end

*********************************************

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
 

cap program drop expandtsfv
program define expandtsfv,rclass
version 14
syntax ,[v(varlist ts fv)] Local(string)
foreach j of local v{
	cap _rmcoll `j', expand noconstant
	local vlist `vlist' `r(varlist)'
}
c_local `local' `vlist'

end


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
*!This program is a minor modification of getcinfo written by David M. Drukker
cap program drop getcinfo2list
program define getcinfo2list, rclass
    syntax, alist(varlist ts fv) blist(varlist ts fv)
    _rmcoll `alist' `blist' , noconstant expand
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
	rmcollv, vlist(`cnames') mat(`mo') 
	local cnames `r(cnames)'

	qui _rmcoll `alist', noconstant expand
	local cname1 `r(varlist)'
	local n1: word count `cname1'
	local n: word count `cnames'
	
	qui _rmcoll `blist', noconstant expand
	local cname2 `r(varlist)'

	local cname1: list cname1 & cnames
	local cname2: list cnames - cname1
	tempname mo1 mo2
	mat `mo1'=	`mo'[1,1..`n1']
	mat `mo2'=`mo'[1,`=`n1'+1'..`n']
	return mat mo = `mo'
	return mat mo1= `mo1'
	return mat mo2= `mo2'
    return local  cnames "`cnames'"
    return local  cname1 "`cname1'"
    return local  cname2 "`cname2'"
end



*!This program is a minor modification of getcinfo written by David M. Drukker
cap program drop getcinfoalist
program define getcinfoalist, rclass
    syntax, alist(varlist ts fv) 
    _rmcoll `alist', noconstant expand
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
/*
local n=colsof(`mat')
local yname `vlist'
forv j=1/`n'{
	
	local r=`mat'[1,`j']
	gettoken v yname: yname
	*disp "`v'"
	if `r'==0 {
		local cnames `cnames' `v'
	}
}
*/

mata: rmcollv("`vlist'","`mat'","cnames")

return local cnames `cnames'
end



/////////////////////////////////////////////////////////////////////////////////
*********************************************************************************
cap program drop xtdplfc_fast
program define xtdplfc_fast, eclass prop(xt)

	version 14


syntax varlist, GENerate(string)  Uvars(varlist) [Zvars(varlist) endox(varlist)  QUANtile TE ///
				   ENDOZflag(numlist integer >0) ivx(varlist)  lags(integer 1) lagyinz(numlist integer >0) ONLYivxz ///
				   Ivtype(numlist integer >=0 )  NKnots(numlist integer >=2) power(numlist integer >0)  ///
				   grid(string) pctile(integer 0) brep(integer 200) ivz(string) ///
				   level(integer 95)  WILD  predict(string) NODOTS ///
				   MINNKnots(numlist integer >=2) MAXNKnots(numlist integer >=2) TENfoldcv]	


		  
*************************************************************
if `lags'<1{
		display as error "No lag dependent variables included."
		display as error "Use xtplfc or ivxtplfc instead."
		exit	
}

**************************************************************
	
	local nz0: word count `zvars'
	local nu: word count `uvars'	
	local comv: list zvars & uvars
	if "`comv'"!=""{
		di as error "`comv' included in both covariates and functions"
		di as error "Check zvars() and uvars()"
		exit 198
	}


	gettoken depvar indepvars : varlist

****************************************************************
****Analyze the time structure	
 	allocy `depvar',lags(`lags') lagyinz(`lagyinz')
    local Dy1inx   `r(Dy1inx)'
    local laby1inx `r(laby1inx)'
    local Doyinx 	`r(Doyinx)'
    local oyinx `r(oyinx)'
    local laboyinx  `r(laboyinx)'
    local laby1inz `r(laby1inz)'
    local oyinz `r(oyinz)'
    local laboyinz `r(laboyinz)' 
    local yinz `r(yinz)'
    local yinx `r(yinx)'

	local nyinz: word count `yinz'
	local nz=`nz0'+`nyinz'
	
	if `nz'!=`nu'{
		di as error "# of variables in uvars() not equal to that specified by zvars() and lagyinz() "
		exit
	}
**********************************************************************
***Add lagyinz in the front of zvars
    
	local zvars `laby1inz' `laboyinz' `zvars'

*************************************************************************

if `"`ivz'"'!=""{
	checkivz `ivz'
	local ivzlist `r(ivlist)'
	local ivuflag `r(uflag)'
	local ivztypelist `r(typelist)'	
} 


foreach v of local ivuflag{
	if `v'>`nz'{
		di as error "Error in ivz()"
		di as error "# specified by uflag() > # of variables in uvars()"
		exit
	}
}

************************************************************************
if "`endozflag'"!=""{
	numlist "`endozflag'",sort
	local endozflag `r(numlist)'
	local endozflag: list uniq endozflag	
}


foreach v of local endozflag{
	if `v'>`nz'{
		di as error "Error in endozflag()/lagyinz()"
		di as error "# specified by endozflag()/lagyinz() > # of variables in uvars()"
		exit
	}
	if `v'<`nyinz'{
		di as error "The first `nyinz' variables in uvars() interact with the lags of depvar specified by lagyinz()"
		di as error "The # in endozflag() should be > `nyinz'"
		exit
	}

}


if "`lagyinz'"!=""{
	local endozflag2 1 `endozflag'	// adding the first lag of depvar
}
else{
	local endozflag2  `endozflag'
}


//local endozfalg2: list uniq endozflag2

*********************************************************************
**********************************************************************

//index the uvars
qui numlist "1/`nu'"
local allflags=r(numlist)

local exgozflag: list allflags - endozflag2

local nendoz: word count `endozflag2'
if "`ivtype'"==""{
	forv i=1/`nendoz'{
		local ivtype `ivtype' 1   // by default
	}
	
}
else{
	local nivtype: word count `ivtype'
	if `nivtype'!=`nendoz'{
		di as error "# of integer in ivtype()!= # of variables specified by endozflag() and lagyiz()"
		exit
	}
}

//mata: ivtype=tokens(st_local("ivtype"))


***********************************************************************

	if `"`endox'"'!="" & "`ivx'"==""{
			disp as error "Instruments should be provided for variables in endox()." 
			disp as error "Using ivx()."
			exit 198
	}


*************************************************************************


*************************************************************************

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



********************************************************************************	
* powermat nknotsmat maxnknotsmat minnknotsmat nknots4cv nsplinesmat nknots4cv


	tempname powermat nknotsmat maxnknotsmat minnknotsmat mxmi nsplinesmat nknots4cv
	if "`power'"==""{
		//disp "By default, the number of power (degree) is set to 3 for all functions."
		mat `powermat'=J(1,`nu',3)
		mata: powermat=st_matrix("`powermat'")
	}
	else{
		local npow: word count `power'
		if `npow'!=`nu' {
			di as error "# of integers in power() != # of variables in uvars()."
			exit			
		}
		mata: powermat=strtoreal(tokens(st_local("power")))
		mata: st_matrix("`powermat'",powermat)	
		
	}	

	
	if "`maxnknots'"=="" & "`nknots'"!="" {

		local nk0: word count `nknots'
		if `nk0'!=`nu'{
			di as error "# of integers in knots() != # of variables in uvars()."
			exit
		}

		mata: nknotsmat=strtoreal(tokens(st_local("nknots")))
		mata: st_matrix("`nknotsmat'",nknotsmat)

	}
	
	
	if "`maxnknots'"=="" & "`nknots'"=="" {
		mat `nknotsmat'=J(1,`nz',2)
		mata: nknotsmat=st_matrix("`nknotsmat'")
		//disp "By default, # of knots is set to 2 for all functions."
	}
	
	if "`maxnknots'"!=""{
		local nmaxnk: word count `maxnknots'
		if `nmaxnk'!=`nu'{
			di as error "# of integers in maxnknots != # of variables in uvars()."
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
			di as error "# of integers in minnknots != # of variables in uvars()."
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
*deal with collinearity
*********************************************************************
/*
	local comv: list indepvars & endox
	if "`comv'"!=""{
		di as error "`comv' included in both exogenous and endogenous variable lists"
		exit 198
	}
*/
	local comv: list indepvars | endox
	local exgox: list comv - endox

	local comv: list comv & zvars

	if !(`"`comv'"'==""){
		disp as red "`comv' included in both linear and nonlinear variable lists." 
		exit 198
	}


	
	if `"`endox'"'!=""&`"`exgox'"'!=""{
	 		getcinfo2list, alist(`endox') blist(`exgox')
	 		local endox `r(cname1)'
	 		local exgox `r(cname2)'
	 	}
	 else if `"`endox'"'==""&`"`exgox'"'!=""{
	 		getcinfoalist, alist(`exgox') 
	 		local exgox `r(cnames)'
	 	}
	 else if `"`endox'"'!=""&`"`exgox'"'==""{
	  		getcinfoalist, alist(`endox') 
	 		local endox `r(cnames)'	
	 	}

 


	if "`ivx'"!=""{
	 		getcinfoalist, alist(`ivx') 
	 		local ivx `r(cnames)'
		}


	local xvars `endox' `exgox'

/////////////////////////////////////////////////////////////////////
* construct IVs

	if "`onlyivxz'"==""{
		local ivtypecopy `ivtype'
		if "`laby1inz'"!=""{
			local ivzlist L2.`depvar' `ivzlist'
			local ivuflag  1 `ivuflag'
			gettoken ivtemp ivtypecopy: ivtypecopy
			local ivztypelist `ivtemp' `ivztypelist' 
		}
		else{
			local ivx L2.`depvar' `ivx'
		}

		foreach v of local endox{
			local ivx `ivx' L2.`v'
		}

		foreach v of local endozflag{
			local zv: word `v' of `zvars'
			local ivzlist `ivzlist' L2.`zv'
		}

		local ivuflag `ivuflag' `endozflag'
		local ivztypelist `ivztypelist' `ivtypecopy'

	}

	local ivlist `ivx'
	local ivzlistcopy `ivzlist'
	local ivztypelistcopy `ivztypelist'
	foreach uk of local ivuflag{
		gettoken zk ivzlistcopy: ivzlistcopy
		gettoken ivk ivztypelistcopy: ivztypelistcopy
		local uk0: word `uk' of `uvars'
		local ivname `zk'*L`ivk'.[Splines(`uk0')]
		local ivname=subinstr(`"`ivname'"',"L0.[","L.[",.)
		local ivlist `ivlist' `ivname'
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
		local gridname `gridname' `grid'_`k'
		label var `grid'_`k' "Grid points of `uk' over [`pctile'-`upctile'] percentiles"	
		local k=`k'+1	
	}
  }
  
  

 /*   
	local zvars2name
	
	forv uk=1/`nu'{
		
		local zk: word `uk' of `zvars'
		local uk0: word `uk' of `uvars'
		if `"`zk'"'==`"L.`depvar'"'{
			local zkname L1_`depvar'
		}
		else{
			local zkname=strtoname(`"`zk'"',1)	
		}
		local zvars2name `zvars2name' `zkname'
		confirm new variable `generate'_`zkname'
		confirm new variable `generate'_`zkname'_sd
		qui gen `generate'_`zkname'=.
		qui gen `generate'_`zkname'_sd=.	
		label var `generate'_`zkname' `"Estimated g(`uk0') at the `vlabel' points w.r.t `zk'*g(`uk0')"'	
		label var  `generate'_`zkname'_sd  `"`vcetype' S.E. of g(`uk0') w.r.t `zk'*g(`uk0')"'
	
	}	
*/
	forv uk=1/`nu'{
		
		local zk: word `uk' of `zvars'
		local uk0: word `uk' of `uvars'
		confirm new variable `generate'_`uk'
		confirm new variable `generate'_`uk'_sd
		qui gen `generate'_`uk'=.
		qui gen `generate'_`uk'_sd=.	
		label var `generate'_`uk' `"Estimated g(`uk0') at the `vlabel' points w.r.t `zk'*g(`uk0')"'	
		label var  `generate'_`uk'_sd  `"`vcetype' S.E. of g(`uk0') w.r.t `zk'*g(`uk0')"'
	
	}	
    
  
*********************************************************************	
/*
	_xt, trequired 
	local id=r(ivar)
	local time=r(tvar)
*/	

/*
//fill the gap
	tempvar dataflag
	qui gen `dataflag'=1
	qui tsfill,full

//continuous id and time
	tempvar cid ctime
	qui egen `cid'=group(`id')
	qui egen `ctime'=group(`time')
	qui tab `ctime',nofreq
	local ntime=r(r)
	local te=cond("`te'"=="",0,`ntime'-1)
*/

	 qui xtset
	 if r(balanced)!="strongly balanced" {
		noi display as error "{bf:xtplfc} needs strongly balanced panel, use {bf:tsfill, full} to rectangularize your data!."
		exit	 	
	 }

	local id=r(panelvar)    
    local time=r(timevar) 
	local tdelta=r(tdelta)
	//continuous id and time
	tempvar cid ctime
	qui egen `cid'=group(`id')
	qui egen `ctime'=group(`time')

/*	
	// Rectangularize the data
	// Note: For survey data, there might be year gaps which must not be filled in.
	// re-xtset before using tsfill
	qui xtset `cid' `ctime' 
	tempvar dataflag
	qui gen `dataflag'=1
	qui tsfill,full		
*/	
*********************************************************************
   sort `id' `time'

	qui count if `time'-`time'[_n-1]!=`tdelta' & `cid'==`cid'[_n-1]
	local gaps=r(N)
	if `gaps'!=0{
		noi di as error `"`time' is regularly spaced, but does not have intervals of 1"'
		exit
	}  
	
   qui count if `cid'==`cid'[1]
   local ntime=r(N) 
   local te=cond("`te'"=="",0,`ntime'-1)
   local labtematrix
   forv j=1/`te'{
		local temp=`time'[`j']
		local labtematrix `labtematrix' _I`time'_`temp'
   }
   
*********************************************************************
*********************************************************************
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
		
	
*********************************************************************
	**********************************

	local yX `depvar' `laby1inx' `laboyinx' `xvars'
	local exgox `laboyinx' `exgox'

	local didots=cond("`nodots'"=="",1,0)
	local qflag=cond("`quantile'"=="",1,0)						
	tempname b V


   // qui xtset
	di 
/*	
	mata: xtdplfc("`cid'","`ctime'","`yX'","`zvars'","`zvars2name'","`uvars'","`exgox'","`exgozflag'","`ivx'","`ivzlist'", ///
	              "`ivuflag'","`ivztypelist'","`nknotsmat'", "`maxnknotsmat'", "`minnknotsmat'","`powermat'",`qflag',`te',`brep',"`V'","`b'", ///
				  "`generate'", "`grid'",`didots',"`wild'", "`yinx'","`yinz'","`onlyivxz'")
*/
/*
	mata: xtdplfc("`cid'","`ctime'","`yX'","`zvars'","`uvars'","`exgox'","`exgozflag'","`ivx'","`ivzlist'", ///
	              "`ivuflag'","`ivztypelist'","`nknotsmat'", "`maxnknotsmat'", "`minnknotsmat'","`powermat'",`qflag',`te',`brep',"`V'","`b'", ///
				  "`generate'", "`grid'",`didots',"`wild'", "`yinx'","`yinz'","`onlyivxz'")
*/
    qui xtset `id' `time'
	mata: xtdplfc("`cid'","`ctime'","`yX'","`zvars'","`uvars'","`exgox'","`exgozflag'","`ivx'","`ivzlist'", ///
	              "`ivuflag'","`ivztypelist'","`nknotsmat'", "`maxnknotsmat'", "`minnknotsmat'","`powermat'",`qflag',`te',`brep',"`V'","`b'", ///
				  "`generate'", "`gridname'",`didots',"`wild'", "`yinx'","`yinz'","`onlyivxz'","`splitdummy'")

				  
				  
    qui replace `touse'=0 if missing(`touse')
	tempname nsplinesmat
	mat `nsplinesmat'=`nknotsmat'+`powermat'
	local nobs= `nobs'
	local dof=`DFR'	  
	
	local nx: word count `laby1inx' `laboyinx' `xvars' `labtematrix'
	if `nx'>0{
		tempname b1 V1
		mat `b1'=`b'[1,1..`nx']
		mat `V1'=`V'[1..`nx',1..`nx']
		mat colnames `b1'= `laby1inx' `laboyinx' `xvars' `labtematrix'
		mat colnames `V1'= `laby1inx' `laboyinx' `xvars' `labtematrix'
		mat rownames `V1'= `laby1inx' `laboyinx' `xvars' `labtematrix'
	
	}


/*
	local labz
	local kk=1
	foreach z of local zvars{
		local uk: word `kk' of `uvars'
		local nspline=`nsplinesmat'[1,`kk']
		forv j=1/`nspline'{
			local labz `labz' `z'*Sp_`uk'_`j'
		
		}
		local kk=`kk'+1
	
	}
*/

local labz
forv uk=1/`nu'{
		
		local zk: word `uk' of `zvars'
		local uk0: word `uk' of `uvars'
		local nspline=`nsplinesmat'[1,`uk']
		forv j=1/`nspline'{
			local labz `labz' `zk'*S(`uk0')_`j'
		}


	}	



	
	mat colnames `b'= `laby1inx' `laboyinx' `xvars' `labtematrix' `labz'
	mat colnames `V'= `laby1inx' `laboyinx' `xvars' `labtematrix' `labz'
	mat rownames `V'= `laby1inx' `laboyinx' `xvars' `labtematrix' `labz'
	//mat list `b'

	//qui keep if `dataflag'==1
	//qui xtset `id' `time'
	sort `id' `time'

	noi disp in green""
	noi disp in green "Partially linear functional-coefficient panel data model."
	noi disp in green "Fixed-effect series semiparametric estimation."
    //noi di in green "" 
if (`nx'>0){	
		ereturn local vcetype `vcetype'
		//mat list `b1'
		//mat list `V1'
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
		noi di in green "{col 48} Root MSE             =" in yellow %8.4f `RMSE'		
		ereturn display, level(`level') 		
	}
	else{
		ereturn post, dep(`depvar') obs(`nobs') esample(`touse')	buildfvinfo
	
	}

		*noi di in green "Instrumented: `endox' `lablHZ2'"
		*noi di in green "Instruments: `ivlist'"
		noi disp in green "Estimated functional coefficient(s) are saved in the data."
		noi disp in green "Instruments for differenced equation: `ivlist'"
		//ereturn clear
		//ereturn post  `b' `V',esample(`touse') dep(`depvar') obs(`nobs') dof(`dof')
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
		//ereturn mat    knotv  = `knotv'
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
		ereturn local cmd "xtdplfc" 
		ereturn local title "Partially linear functional-coefficient dynamic panel data model"
		ereturn local model "Fixed-effect Sieve 2SLS Estimation" 
		ereturn local depvar "`depvar'" 
		ereturn local estfun `estfun'
		ereturn local basis "B-spline"
		ereturn local vcetype `vcetype'
		ereturn local ivlist `ivlist'
		

	
	
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

////////////////////////////////////////////////////////////////////////////
****************************************************************************


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




