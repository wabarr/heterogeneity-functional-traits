Code for paper submitted to Journal of Human Evolution

Run the lines of code below in the R terminal to run all analyses. This code assumes the operating system sorts the files numerically, which it does on Mac OS X. If your OS doesn't do this you may need to run the files manually

```
analysis_dir <- "~/Dropbox/heterogeneity_funcric_JHE/JHE_revise_resubmit/heterogeneity-functional-traits/"
setwd(analysis_dir)
lapply(
  list.files(".", pattern="R$", full=T),
  FUN=source
)
```
