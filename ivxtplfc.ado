*! version 7.0
* 2019-8-26
* change cvtwo option to tenfoldcv
* bugs fix
*! version 6.0.5 
* 2019-01-25
*12Dec2018
*15Nov2018
*By Kerry Du
cap program drop ivxtplfc
program define ivxtplfc, eclass prop(xt)

	version 14

if replay()& "`e(cmd)'"=="ivxtplfc" {
	ereturn  display
	exit
}	

syntax varlist, GENerate(string)  Uvars(varlist) Zvars(varlist)[ endox(varlist) ///
				   ENDOZflag(numlist integer >0) ivx(varlist)  ivz(string)  TE  ///
				    NKnots(numlist integer >=2) power(numlist integer >0)  ///
				   grid(string) pctile(integer 0) brep(integer 200) QUANtile ///
				   level(integer 95)  WILD spgen(string) predict(string) ///
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
/*
	capture which lrcov
	local rc = _rc
	if (`rc') {
		display as error "{bf:xtplfc} requires {bf:lrcov}"
		exit 198
		}				   
*/
*************************************************************
if "`fast'"!=""{
	ivxtplfc_fast `varlist', generate(`generate')  uvars(`uvars') zvars(`zvars') ///
				   endox(`endox') endozflag(`endozflag') ivx(`ivx')  ivz(`ivz')  ///
				   nknots(`nknots') power(`power') `quantile' ///
				   grid(`grid') pctile(`pctile') brep(`brep')  `te' ///
				   level(`level')  `wild' spgen(`spgen') predict(`predict') ///
				   minnknots(`minnknots') maxnknots(`maxnknots') `nodots' `tenfoldcv'
		exit
}				   
	


*************************************************************
	
	local nz: word count `zvars'
	local nu: word count `uvars'	
	local comv: list zvars & uvars
	if "`comv'"!=""{
		di as error "`comv' included in both covariates and functions"
		di as error "Check zvars() and uvars()"
		exit 198
	}
	if `nz'!=`nu'{
		di as error "# of variables in uvars() not equal to that specified by zvars()"
		exit
	}	


	gettoken depvar indepvars : varlist

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

qui numlist "1/`nu'"
local allflags=r(numlist)

local exgozflag: list allflags - endozflag

if (`"`endox'"'!="" & "`endozflag'"!=""){
	di as error "No endogeous variables specified"
	di as error "Use {bf:xtplfc} instead"
	exit
	
}

	if `"`endox'"'!="" & "`ivx'"==""{
			disp as error "Instruments should be provided for variables in endox()." 
			disp as error "Using ivx()."
			exit 198
	}




*************************************************************************

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
	local nx: word count `xvars'

/////////////////////////////////////////////////////////////////////
* construct IVs

	local ivlist `ivx'
	local ivzlistcopy `ivzlist'
	local ivztypelistcopy `ivztypelist'
	foreach uk of local ivuflag{
		gettoken zk ivzlistcopy: ivzlistcopy
		gettoken ivk ivztypelistcopy: ivztypelistcopy
		local uk0: word `uk' of `uvars'
		local ivname `zk'*L`ivk'.[Splines(`uk0')]
		local ivlist `ivlist' `ivname'
	}		
*********************************************************************	
  local vlabel=cond("`grid'"=="","observed","grid")

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
  
  

	
	forv uk=1/`nu'{
		
		local zkname: word `uk' of `zvars'
		local uk0: word `uk' of `uvars'
		//local zvars2name `zvars2name' `zkname'
		confirm new variable `generate'_`uk'
		confirm new variable `generate'_`uk'_sd
		qui gen `generate'_`uk'=.
		qui gen `generate'_`uk'_sd=.	
		label var `generate'_`uk' `"Estimated g(`uk0') at the `vlabel' points w.r.t `zkname'*g(`uk0')"'	
		label var  `generate'_`uk'_sd  `"`vcetype' S.E. of g(`uk0') w.r.t `zkname'*g(`uk0')"'
	
	}	

  
 local qflag=cond("`quantile'"=="",1,0)	 
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
	
	//qui tab `cid', nofreq
	//local nid=r(r)
	//qui tab `ctime',nofreq
	//local ntime=r(r)
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
	
	qui count if `ctime'==`ctime'[1]
	local nid=r(N)
	qui count if `cid'==`cid'[1]
	local ntime=r(N)
	
	local te=cond("`te'"=="",0,`ntime'-1)
	
	sort `id' `time'
	qui count if `ctime'-`ctime'[_n-1]!=`tdelta' & `cid'==`cid'[_n-1]
	local gaps=r(N)
	if `gaps'!=0{
		noi di as error `"`time' is regularly spaced, but does not have intervals of 1"'
		exit
	}  		
	
	
	if `te'>0{
		mata: tematrix=gendtime(`nid',`te')
		mata: _sttempvar("tematrix",tematrix)
	}
	
   local labtematrix
   forv j=1/`te'{
		local temp=`time'[`j']
		local labtematrix `labtematrix' _I`time'_`temp'
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
	local nx2=`nx'+`te'
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
					foreach uk of local endozflag {
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



					  qui ivregress 2sls D.`depvar' (`Dendox' `DENDZ' = `ivx' `instz' )  ///
								   `DEXGZ'  `Dexgox2' if `cid'~=`jt', noconstant 


						qui predict `yhat' if `cid'==`jt',xb
						qui replace `ehat_sq'=(D.`depvar'-`yhat')^2 if `cid'==`jt'
						qui cap drop `yhat'
					}
				}
				else{
					qui xtset `id' `time'
				   forvalues jt=1/10{


					  //qui putmata tempiv=`ivx'  if `splitdummy'~=`jt', replace
					  //mata: tempiv[1::10,.]
					  qui ivregress 2sls D.`depvar' (`Dendox' `DENDZ' = `ivx' `instz' )  ///
								   `DEXGZ'  `Dexgox2' if `splitdummy'~=`jt', noconstant 


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
			//local udist=(`max_`ii''-`min_`ii'')/`nknot`ii''
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
				foreach uk of local endozflag {
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

//su  `EXGZ' 
//su `Dendox' `DENDZ' `DEXGZ'  `Dexgox'
//su `ivx' `instz'
				
	qui ivregress 2sls D.`depvar' (`Dendox' `DENDZ' = `ivx' `instz' )  ///
	 					   `DEXGZ'  `Dexgox2', noconstant 			


   	local Xalls `Dendox' `Dexgox2'
	//local nbxs: word count `Xalls'
	local allreg  `Dendox' `DENDZ' `DEXGZ' `Dexgox2'
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
	mata: XXX=st_data(.,"`endox' `ENDZ' `EXGZ' `exgox'")
	mata: b00=b0[1,1..(length(b0)-`te')]
	mata: yhat=XXX*b00'
	if `te'>0{
		mata: timeeffects=b0[1,(length(b0)-`te'+1)..length(b0)],-sum(b0[1,(length(b0)-`te'+1)..length(b0)])
		mata: yhat=yhat+J(`nid',1,timeeffects')
	}

	qui getmata `yhat'=yhat, replace
	qui gen `ehat'=`depvar'-`yhat'
	qui bys `cid' (`ctime'): egen `aui'=mean(`ehat') 
	qui replace `ehat'=`ehat'-`aui'
	qui replace `yhat'=`yhat'+`aui' 
	
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
	
	mata: bz=b0[1,zpos]
	
	local estfun
	mata: gfun=estgfun(bz',HH,bzindex)
	
  **************************************

   if "`grid'"!=""{	
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
	  mata: gfun=estgfun(bz',HH,bzindex)
   }
		
	local k=1
	foreach j of local zvars {
		
		mata: temvec=gfun[,`k']
		qui getmata `generate'_`k'=temvec,replace
		local estfun `estfun' `generate'_`k'
		local k=`k'+1
	}
	


   
   
	
*************************************************************	


*************************************************************

if `brep'>0 {
	di
	local regcmd (`Dendox' `DENDZ' = `ivx' `instz' ) `DEXGZ' `Dexgox2', noconstant 
	local bootcmd bs_xtresid2, id(`cid') xb(`Dxb') residual(`Dehat') regcmd(`regcmd') 
	if "`wild'"!=""{
		local bootcmd bs_xtwild2, xb(`yhat') residual(`ehat') regcmd(`regcmd') 
	}
	disp "Computing the bootstrap standard errors..."
	
	mata: bsmat=J(`brep',length(b0),.)


	forvalues z = 1/`nz'{
		mata: gf_`z'=J(`N',`brep',.)
	}
	
		
	mata: bsmat=J(`brep',length(b0),.)	
	
	forvalues z = 1/`nz'{
		mata: gf_`z'=J(`N',`brep',.)
	}
	
	sort `id' `time'
	//qui xtset
	tempname bb 
	forv k=1/`brep' {
			
		`nodots' displaydot, dot(`k')
		qui `bootcmd'
 
		mat `bb'=e(bs)
		mata: bb=st_matrix("`bb'")
		mata: bs2=bb[1,zpos]
		mata: bsmat[`k',.]=bb[1,xpos],bs2
			forvalues z=1/`nz'{
				mata: bz=bs2[1,(1+bzindex[1,`z'])..bzindex[1,`z'+1]]
				mata: gf_`z'[.,`k']=HH[.,(1+bzindex[1,`z'])..bzindex[1,`z'+1]]*bz'

			}
		

   }
	


	mata: V=quadvariance(bsmat)
	mata: st_matrix("`V'",V)
	if "`Xalls'"!=""{
		mat `V1'=`V'[1..`nx2',1..`nx2']
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
		mat `V1'=`V'[1..`nx2',1..`nx2']
	}	
}

/*
else{
	
		   tempname omega
		   qui sort `id' `time'
		   local WW `ivx' `instz' `DEXGZ' `Dexgox2'
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
			   
			   qui xtlrcov `WW', id(`cid') time(`ctime') wvar(`Dehat') `hacopt' 
			   mat `omega'=r(omega)
			   mata: omega=st_matrix("`omega'")	
		   }


	    mata: omega=omega/rows(W)
	

		if "`Xalls'"!=""{

			sort `id' `time'
			mata: X=st_data(.,"`Xalls'","`flag'")
			mata: gsd=estcov(W,(X,HZ),omega,"`V'",`nx2',bzindex,HH)
			mat `V1'=`V'[1..`nx2',1..`nx2']

			
		}
		else{
		
			mata: gsd=estcov(W,HZ,omega,"`V'",`nx2',bzindex,HH)
		
		}
	
		
	  //// qui getmata `generate'_sd*=gsd,replace //////
	   local k=1

	   foreach v of local zvars{
	   		mata: tempgsd=gsd[.,`k']
			qui getmata `generate'_`v'_sd=tempgsd,replace
			local k=`k'+1
	   
	   }
	   	   


}
*/



	//qui keep if `dataflag'==1
	//qui xtset `id' `time'
	sort `id' `time'
	
	mat colnames `B' =  `endox' `exgox' `labtematrix' `lablHZ' 
	mat colnames `V' =  `endox' `exgox' `labtematrix' `lablHZ'  
	mat rownames `V' =  `endox' `exgox' `labtematrix' `lablHZ'  	

	noi disp in green""
	noi disp in green "Partially linear functional-coefficient dynamic panel data model."
	noi disp in green "Fixed-effect sieve 2SLS estimation."

	if ("`Xalls'"!=""){
			mat colnames `b1' = `endox' `exgox' `labtematrix'
			mat colnames `V1' = `endox' `exgox' `labtematrix'
			mat rownames `V1' = `endox' `exgox'  `labtematrix'
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
		//ereturn clear
		//ereturn post  `B' `V',esample(`flag') dep(`depvar') obs(`nobs') dof(`dof')
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
		ereturn local cmd "ivxtplfc" 
		ereturn local title "Partially linear functional-coefficient panel data model"
		ereturn local model "Fixed-effect Sieve 2SLS Estimation" 
		ereturn local depvar "`depvar'" 
		ereturn local estfun `estfun'
		ereturn local basis "B-spline"
		ereturn local vcetype `vcetype'
		ereturn local ivlist `ivlist'		
		
		forv j=1/`nu'{
			ereturn local k`j' `knotlist_`j''
		}
				  		
	

	

end



////////////////////////////////////////////////////////////////////
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

cap program drop  bs_xtresid2
program define bs_xtresid2,eclass
        version 14
        syntax, id(varname numeric) xb(varname numeric) RESidual(varname numeric) regcmd(string) 

        tempvar idx ys
		* idx randomly selects the observations with replacement
        qui putmata id=`id',replace
        mata: mm_panels(id,idinfo=.)
       // mata: length(idinfo)
        mata: idx=mm_sample(length(idinfo),.,idinfo)
        qui getmata `idx'=idx, replace

        * the new dependent variable using resample residuals
        qui gen double `ys' = `xb' + `residual'[`idx']
		
		qui ivregress 2sls `ys' `regcmd'
		
        mat b=e(b)
        ereturn matrix bs=b
end  



cap program drop  bs_xtwild2
program define bs_xtwild2,eclass
        version 14
        syntax,  xb(varname numeric) RESidual(varname numeric) regcmd(string)

		tempvar randu ys flag Dys
		
		qui gen `randu'=uniform()
		qui gen `flag'=(1-sqrt(5))/2
		qui replace `flag'=(sqrt(5)+1)/2 if `randu'>((1+sqrt(5))/(2*sqrt(5)))
		
		qui gen `ys'=`xb'+`residual'*`flag'
		
		qui ivregress 2sls D.`ys' `regcmd'
		//list D.`ys'
		
        mat b=e(b)
        
        ereturn matrix bs=b
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

	rmcollv, vlist(`cnames') mat(`mo') 
	local cnames `r(cnames)'
	
	qui _rmcoll `alist', noconstant expand
	local cname1 `r(varlist)'
	*local n1: word count `cname1'
	*local n: word count `cnames'
	
	qui _rmcoll `blist', noconstant expand
	local cname2 `r(varlist)'

	local cname1: list cname1 & cnames
	local cname2: list cnames - cname1

    return local  cnames "`cnames'"
    return local  cname1 "`cname1'"
    return local  cname2 "`cname2'"
end




///////////////////////////////////////////////////////////


cap program drop ivxtplfc_fast
program define ivxtplfc_fast, eclass prop(xt)

	version 14

syntax varlist, GENerate(string)  Uvars(varlist) Zvars(varlist)[ endox(varlist) ///
				   ENDOZflag(numlist integer >0) ivx(varlist)  ivz(string) TE  ///
				   NKnots(numlist integer >=2) power(numlist integer >0)  ///
				   grid(string) pctile(integer 0) brep(integer 200) QUANtile ///
				   level(integer 95)  WILD spgen(string) predict(string) ///
				   MINNKnots(numlist integer >=2) MAXNKnots(numlist integer >=2) NODOTS TENfoldcv]	
				   
	
	local nz: word count `zvars'
	local nu: word count `uvars'	
	local comv: list zvars & uvars
	if "`comv'"!=""{
		di as error "`comv' included in both covariates and functions"
		di as error "Check zvars() and uvars()"
		exit 198
	}
	if `nz'!=`nu'{
		di as error "# of variables in uvars() not equal to that specified by zvars()"
		exit
	}	


	gettoken depvar indepvars : varlist

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

qui numlist "1/`nu'"
local allflags=r(numlist)

local exgozflag: list allflags - endozflag

if (`"`endox'"'!="" & "`endozflag'"!=""){
	di as error "No endogeous variables specified"
	di as error "Use {bf:xtplfc} instead"
	exit
	
}

	if `"`endox'"'!="" & "`ivx'"==""{
			disp as error "Instruments should be provided for variables in endox()." 
			disp as error "Using ivx()."
			exit 198
	}




*************************************************************************

	
	local vcetype=cond("`wild'"=="","Clustered Bootstrap","Wild Bootstrap")
	
/*	
	else{
		local vcetype Robust 
	}	


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

	local ivlist `ivx'
	local ivzlistcopy `ivzlist'
	local ivztypelistcopy `ivztypelist'
	foreach uk of local ivuflag{
		gettoken zk ivzlistcopy: ivzlistcopy
		gettoken ivk ivztypelistcopy: ivztypelistcopy
		local uk0: word `uk' of `uvars'
		local ivname `zk'*L`ivk'.[Splines(`uk0')]
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
		label var `grid'_`k' "Grid points of `uk' over [`pctile'-`upctile'] percentiles"	
		local gridname `gridname' `grid'_`k'		
		local k=`k'+1	
	}
  }
  
  

	
	forv uk=1/`nu'{
		
		local zkname: word `uk' of `zvars'
		local uk0: word `uk' of `uvars'
		//local zvars2name `zvars2name' `zkname'
		confirm new variable `generate'_`uk'
		confirm new variable `generate'_`uk'_sd
		qui gen `generate'_`uk'=.
		qui gen `generate'_`uk'_sd=.	
		label var `generate'_`uk' `"Estimated g(`uk0') at the `vlabel' points w.r.t `zkname'*g(`uk0')"'	
		label var  `generate'_`uk'_sd  `"`vcetype' S.E. of g(`uk0') w.r.t `zkname'*g(`uk0')"'
	
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
   


local yX `depvar' `endox' `exgox'	
tempname V b
local didots=cond("`nodots'"=="",1,0)
local qflag=cond("`quantile'"=="",1,0)	
					
di

/*

 mata: ivxtplfc("`cid'","`ctime'","`yX'","`zvars'","`uvars'","`exgox'","`exgozflag'", ///
		          "`ivx'","`ivzlist'", "`ivuflag'","`ivztypelist'","`nknotsmat'", "`maxnknotsmat'", ///
		          "`minnknotsmat'","`powermat'",`qflag',`te',`brep',"`V'","`b'", "`generate'", "`gridname'",`didots', ///
		         "`wild'")
*/
 //qui xtset
 mata: ivxtplfc("`cid'","`ctime'","`yX'","`zvars'","`uvars'","`exgox'","`exgozflag'", ///
		          "`ivx'","`ivzlist'", "`ivuflag'","`ivztypelist'","`nknotsmat'", "`maxnknotsmat'", ///
		          "`minnknotsmat'","`powermat'",`qflag',`te',`brep',"`V'","`b'", "`generate'", "`gridname'",`didots', ///
		         "`wild'","`splitdummy'")				 
				  
	qui replace `touse'=0 if missing(`touse')
	local nx: word count `xvars' `labtematrix'
	//local nx=`nx'+`te'
	mat `nsplinesmat'=`nknotsmat'+`powermat'
	local nobs= `nobs'
	local dof=`DFR'	   
	//local nx: word count `xvars'
	if `nx'>0{
		tempname b1 V1
		mat `b1'=`b'[1,1..`nx']
		mat `V1'=`V'[1..`nx',1..`nx']
		mat colnames `b1'= `xvars' `labtematrix'
		mat colnames `V1'= `xvars' `labtematrix'
		mat rownames `V1'= `xvars' `labtematrix'
	
	}
	local labz
	local kk=1
	foreach z of local zvars{
		local uk: word `kk' of `uvars'
		local nspline=`nsplinesmat'[1,`kk']
		forv j=1/`nspline'{
			local labz `labz' `z'*Sp(`uk')_`j'
		
		}
		local kk=`kk'+1
	
	}
	
	mat colnames `b'= `xvars' `labtematrix' `labz'
	mat colnames `V'= `xvars' `labtematrix' `labz'
	mat rownames `V'= `xvars' `labtematrix' `labz'
	//mat list `b'

	//qui keep if `dataflag'==1
	//qui xtset `id' `time'
	sort `id' `time'
	noi disp in green""
	noi disp in green "Partially linear functional-coefficient dynamic panel data model."
	noi disp in green "Fixed-effect sieve 2SLS estimation."

	if (`nx'>0){
			//mat colnames `b1' = `xvars' `labtematrix'
			//mat colnames `V1' = `xvars' `labtematrix'
			//mat rownames `V1' = `xvars' `labtematrix'
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
			ereturn post, dep(`depvar') obs(`nobs') esample(`touse')
		
		}		

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
		ereturn local cmd "ivxtplfc" 
		ereturn local title "Partially linear functional-coefficient panel data model"
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
