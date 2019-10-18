# xtplfc_Stata
Stata module for estimating partially linear functional-coefficient panel data models

- **xtplfc.ado** estimates partially linear functional-coefficient static panel data models.
- **ivxtplfc.ado** estimates partially linear functional-coefficient static panel data models with endogeneous variables.
- **xtdplfc.ado** estimates partially linear functional-coefficient dynamic panel data models.


## Install within Stata
* Method 1: **net install**
```
net install xtplfc, from("https://raw.githubusercontent.com/kerrydu/xtplfc_Stata/master/")
```
* Method 2: **github install**
```
net install github, from("https://haghish.github.io/github/") 
github install kerrydu/xtplfc_Stata
```



## Requirement 
- Stata version 14 or later
- Depends on the `moremata` and `bspline` packages
- Only be used if data are declared to be panel data through the xtset or tsset command



## Citation

If you use this module, please cite  the following papers:

```bibtex
@TechReport {xtplfc2019,
  Author = {Du, K., and Zhang, Y. and Zhou, Q.},
  Title = {Estimating partially linear functional-coefficient panel data models with Stata},
  Note = {Working Paper},
  Year = {2019},
}
```
> Du, K., Zhang, Y. and Zhou, Q. 2019. Estimating partially linear functional-coefficient panel data models with Stata.
> Working Paper.
> https://github.com/kerrydu/xtplfc_Stata/blob/master/manuscript.pdf

```bibtex
@TechReport {zhang2018,
  Author = {Zhang, Y., and Q. Zhou.},
  Title = {Partially linear functional-coefficient panel data models: Sieve Estimation and Specification testing},
  Note = {Working Paper},
  Year = {2018},
}
```
> Zhang, Y., and Q. Zhou. 2018. Partially linear functional-coefficient panel data models: Sieve Estimation and Specification testing. Working paper.



####  Things to be aware of:
- Sometimes  `lxtplfc.mlib` might not be indexed automatically. Consequently, the program would abort with error information that mata function could not be found. In this case, please prompt `mata mata mlib index` in the Stata and re-run the command. 
- If you find any issues during using this module, please report issues on [**Github**](https://github.com/kerrydu/xtplfc_Stata/issues), or email to [**Kerry Du**](https://kerrydu.github.io/) 
