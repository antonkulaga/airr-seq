if("BiocInstaller" %in% rownames(installed.packages())) remove.packages("BiocInstaller")
install.packages("BiocManager")

to_install <- c(
  "tidyverse", # The tidyverse is an opinionated collection of R packages designed for data science.
  "plotly", #charts
  "docopt", #for CLI
  "magick" # for image processing
)
install.packages(to_install, dependencies = TRUE)

BiocManager::install("ShortRead", update=FALSE, ask=FALSE)
