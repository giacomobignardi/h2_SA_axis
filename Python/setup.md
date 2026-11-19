# Python setup (reticulate, RStudio)

`Python/03_FC_G.py` (individual functional gradients) runs in Python from RStudio through
**reticulate**, using a project-local virtual environment in `.venv/`. The environment itself
is not tracked by git; package versions are pinned in `Python/requirements.txt`.

## Before you start
- Recent R and RStudio.
- macOS: Apple's Command Line Tools (Terminal: `xcode-select --install`).
- Open `SA_axis.Rproj`, so the working directory is the repository root.
- Data folders are not part of the repository. The scripts expect this layout:

```
MPI_CBS_SA_SEM_part1/
├── 10_github/                      # this repository
├── 03_rawData/
│   ├── noa/                        # subj_inclusion.csv, HCP behavioural files, ...
│   └── private/rsfMRIavg/          # day1/, day2/, day_avg/  -> <subject>.csv (FC matrices)
└── 11_noa_output/                  # outputs read and written by the scripts
```
HCP data are available from the Human Connectome Project under its data use terms.

## First-time setup (R console)
1. Update reticulate, then restart R (Session → Restart R):
   ```r
   install.packages("reticulate")
   ```
2. Install Python 3.11 and create the environment in the project folder:
   ```r
   library(reticulate)
   py_bin <- install_python(version = "3.11")
   virtualenv_create("./.venv", python = py_bin)
   ```
3. Install the pinned packages:
   ```r
   system2(".venv/bin/python", c("-m", "pip", "install", "-r", "Python/requirements.txt"))
   ```
4. Restart R. The project `.Rprofile` points reticulate to `.venv` automatically. Check:
   ```r
   reticulate::py_config()                          # python: <repo>/.venv/bin/python
   reticulate::py_module_available("brainspace")    # TRUE
   ```

On Windows, I think replacing `.venv/bin/python` with `.venv/Scripts/python.exe`, should be sufficient. But I am not sure. 

## Running the script
```r
reticulate::source_python("Python/03_FC_G.py")
```
Or from a Terminal opened in the repository root: `.venv/bin/python Python/03_FC_G.py`.

## Updating packages
Install or upgrade, then re-pin and commit `Python/requirements.txt`:
```r
system2(".venv/bin/python", c("-m", "pip", "install", "--upgrade", "brainspace"))
system2(".venv/bin/python", c("-m", "pip", "freeze"), stdout = "Python/requirements.txt")
```