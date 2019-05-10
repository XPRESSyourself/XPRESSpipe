# Information for User Testing

### Round 1: v0.1.2-beta
##### Features
- Currently I have only tested most aspects of the pipeline as a pipeline, running an individual sub-modules alone may still be a bit buggy in how the inputs are passed
- This works best on a supercomputer -- if you need access and are a member of the Rutter lab, fill out this [form](https://www.chpc.utah.edu/role/user/account_request.php) and have Jared submit it to CHPC @ helpdesk@chpc.utah.edu or whatever relevant form the application suggests (sorry, I can't view the application since I already have an account)
- Most features of XPRESSpipe are documented [here](https://github.com/XPRESSyourself/XPRESSpipe); however, if you notice something missing or a mistake, please submit an [issue request](https://github.com/XPRESSyourself/XPRESSpipe/issues)

##### Requests
- If you find an issue, it is easiest for me to address if you submit an issue on [GitHub](https://github.com/XPRESSyourself/XPRESSpipe/issues)
  - When filling out an issue, please describe the issue in detail, what you expect as output, and relevant output or log files (that were output in the output directory you specified). You should be able to attach these files when submitting an issue
- If you find a feature is missing you would like to see implemented, submit an [issue]((https://github.com/XPRESSyourself/XPRESSpipe/issues))
- After you have finished playing around with the software, please send me an email or slack describing your impressions, ease of use, etc

##### Installation
- While [Docker Images](https://cloud.docker.com/repository/docker/jordanberg/xpresspipe/general) exist for XPRESSpipe, it may be best to download the [latest version](https://github.com/XPRESSyourself/XPRESSpipe/releases) directly
- Assuming you have conda and git installed on your system, you can install using the ```install.sh``` script included in the [repo](https://github.com/XPRESSyourself/XPRESSpipe/)
  - Use the current one, not the one versioned
- You will also need to install the R dependencies required:
```
$ R

> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("Rsubread", version = "3.8")
> BiocManager::install("dupRadar", version = "3.8")
> BiocManager::install("DESeq2", version = "3.8")
```

##### Getting Started
- Please refer to [this](https://xpresspipe.readthedocs.io/en/latest/content/general-usage.html) page and [this](https://xpresspipe.readthedocs.io/en/latest/content/reference-building.html) page in the docs to help with getting started
