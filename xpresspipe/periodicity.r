#!/usr/bin/env Rscript

# Measure periodicity of footprint sample

# Install dependencies
install.packages("devtools")
library(devtools)
install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE,
				build_opts = c("--no-resave-data", "--no-manual"))
library(riboWaltz)
