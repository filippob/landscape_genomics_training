# landscape_genomics_training

Crash course on basic landscape genomics (+ a touch on the genetic structure of populations)

## Overall structure of the training

1. theoretical background
2. practical session

## Install and using this repository

This repository is structured as a *R project* which uses
[renv](https://rstudio.github.io/renv/articles/renv.html) to manage package
dependencies. Sample scripts are placed in `scripts` folder and can be used
to follow the *landscape genomics course* held in Leon before the SMARTER
annual meeting.

### Cloning project

This project is managed using git. The simplest way to start with this project is 
to clone this repository:

```bash
git clone https://github.com/filippob/landscape_genomics_training.git
```

and then open `landscape_genomics_training.Rproj` using
[Rstudio](https://www.rstudio.com/).

You can also start by opening *Rstudio*, then on the *right-top* corner click on
the *Project* button, then on *New Project* -> *Version Control* -> *Git* and
then by filling the form using <https://github.com/filippob/landscape_genomics_training>
as the repository url.

### Setting up libraries

This project is managed through `renv`, in order to track packages required to
run the scripts inside this folder. All the packages installed here (except one)
are standard r packages which you can find in [CRAN](https://cran.r-project.org/)
and that you could install using the `install.packages` R command. 
The last package you require to access to the SMARTER data is the
[smarterapi](https://github.com/cnr-ibba/r-smarter-api) package, which is not yet
published in cran, is a source package which need to be installed using 
the `devtools::install_github` command from
[devtools](https://cran.r-project.org/web/packages/devtools/index.html) R package.
All this requirements could be managed by `renv` simply by calling:

```r
renv::restore()
```

in the same directory of this R project (where the file `renv.lock` is placed)

## Configure SMARTER backend credentials

SMARTER metadata are not yet public and need to be accessed using the API provided
by the [SMARTER-backend](https://webserver.ibba.cnr.it/smarter-api/docs/). This API
require to generate a *json web token* in order to access SMARTER data. To simplify
this process, this project makes use of the `smarterapi` R packages, which automates
some tasks like the token generation and requests. You require to set up your SMARTER
credentials (which are not stored in this project and **must not be tracked in this
project or any other public repository**). The simplest way to store your SMARTER
credentials is your main `.Rprofile` file in your `$HOME` directory. You can open
such file for editing by calling:

```r
source("scripts/0.set_up_credentials.R")
# or usethis::edit_r_environ()
```

This will open your `$HOME/.Rprofile` file in Rstudio or in your preferred terminal
editor. In this file you have to define your smarter credentials like this:

```text
SMARTER_API_USERNAME=<smarter username>
SMARTER_API_PASSWORD=<smarter password>
```

Then you have to save this file and *restart R* to see effects. After that you 
can use the SMARTER API backend without worring about tokens and authentication.

## Track additional libraries

This project is a github project and should managed using git. You could modify it
according your needs, however since packages are managed using `renv` you require
to track your new dependencies, in order that your scripts
could be shared within this project. The simplest way to ensure that your
environment is consistent with your libraries is by calling:

```r
renv::status()
```

this will tell you if your environment need to be updated or not. If you have
installed or removed packages which are not tracked with this repository, you have
to call:

```r
renv::snapshot()
```

to track new dependencies within the project. The `renv.lock` file need be tracked
like any other scripts in order to share your dependencies with others.
