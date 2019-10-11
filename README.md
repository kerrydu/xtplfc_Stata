# xtplfc_Stata
Stata module for estimating partially linear functional-coefficient panel data models

- **xtplfc.ado** estimates partially linear functional-coefficient static panel data models.
- **ivxtplfc.ado** estimates partially linear functional-coefficient static panel data models with endogeneous variables.
- **xtdplfc.ado** estimates partially linear functional-coefficient dynamic panel data models.


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

> Du, K., Zhang, Y. and Zhou, Q. 2019. "Estimating partially linear functional-coefficient panel data models with Stata"
> Working Paper.
> https://github.com/kerrydu/xtplfc_Stata/blob/master/manuscript.pdf



####  Things to be aware of:
- Sometimes  **lxtplf.mlib** might not be indexed automatically. Consequently, the program would abort with error information that mata function count be found. In this case, please prompt "mata mata mlib index" in the Stata and re-run the command. 
- If you find any issues during using this module, please email to [**Kerry Du**](https://kerrydu.github.io/) (kerrydu@xmu.edu.cn)
