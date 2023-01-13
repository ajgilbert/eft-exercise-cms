# EFT parameterisations for STXS qqH

To generate the EFT parameterisations for STXS qqH using the cards provided in this directory, copy them to `EFT2Obs/cards/qqH-SMEFTsim3` and follow the instructions [here](../README.md). The full setup starting from scratch is described here:

## Process definition

Create a directory `qqH-SMEFTsim3` in `EFT2Obs/cards` and add a `proc_card.dat` specifying the model and defining the process:
```
import model SMEFTsim_topU3l_MwScheme_UFO-massless

generate p p > h j j QCD=0 NP<=1

output qqH-SMEFTsim3 -nojpeg
```
Next, run
```sh
./scripts/setup_process.sh eft-exercise/qqH-SMEFTsim3
```
This creates the process directory `MG5_aMC_v2_6_7/qqH-SMEFTsim3`, and adds two cards to `cards/qqH-SMEFTsim3`: `pythia8_card.dat` and `run_card.dat`. 

In `pythia8_card.dat`, the line `partonlevel:mpi = off` can be uncommented for now. This disables multiparton interactions, which makes the event generation faster.

In `run_card.dat`, two changes have to be made:
- `True = use_syst`: change to `False = use_syst`
- `systematics = systematics_program`: change to `none = systematics_program`

This is necessary because the extra weights that appear when `True = use_syst` is set are not handled correctly by EFT2Obs.


## Identify EFT operators

To automatically detect the EFT operators this process is sensitive to, and set up the reweighting, run
```sh
python scripts/auto_detect_operators.py -p qqH-SMEFTsim3
```
This creates two new files in the `cards/qqH-SMEFTsim3` directory: `config.json` and `reweight_card.dat`. 

The last file we need to add is `param_card.dat`. To create it, run
```sh
python scripts/make_param_card.py -p qqH-SMEFTsim3 -c cards/qqH-SMEFTsim3/config.json \
  -o cards/qqH-SMEFTsim3/param_card.dat
```
At this point, the contents of `cards/qqH-SMEFTsim3` should be exactly the same as the cards provided in this directory.


## Make gridpack and generate events

For this process, the gridpack can be created locally:
```sh
./scripts/make_gridpack.sh eft-exercise/qqH-SMEFTsim3 0 8
```
The file `gridpack_qqH-SMEFTsim3.tar.gz` will be copied to the main EFT2Obs directory. Now we can generate events in a set of slurm jobs:
```sh
python scripts/launch_jobs.py --gridpack gridpack_qqH-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p HiggsTemplateCrossSections -o qqH-SMEFTsim3 --task-name qqH-SMEFTsim3 --dir jobs --job-mode slurm \
  --env "HIGGSPRODMODE=VBF"
```
This runs through the full event generation with `MG5_aMC@NLO`, EFT model reweighting, showering with `Pythia8`, and processing with `Rivet`, using the Rivet routine `HiggsTemplateCrossSections.cc` in the `RivetPlugins` directory. The environment variable `HIGGSPRODMODE=VBF` is needed for the `HiggsTemplateCrossSections` Rivet routine.


## EFT parameterisation

The output of the previous command is 50 yoda files containing all the Rivet routine histograms with a copy for each weight. Merge the yoda files:
```sh
yodamerge -o qqH-SMEFTsim3/RivetTotal.yoda qqH-SMEFTsim3/Rivet_* --no-veto-empty
```
Then use the script `get_scaling.py` to calculate the EFT scaling parameters $A_{i}$ and $B_{ij}$:
```sh
python scripts/get_scaling.py -i qqH-SMEFTsim3/RivetTotal.yoda -o scaling_qqH-SMEFTsim3 \
  --hist "/HiggsTemplateCrossSections/HTXS_stage1_2_pTjet30" --bin-labels eft_exercise_bin_labels.json \
  -c cards/qqH-SMEFTsim3/config.json --rebin 18,19,20,21,22,23,24,25,26,27,28,29  
```
