# EFT parameterisations for Higgs STXS

To generate the EFT parameterisations for the ttbar measurement using the cards provided in this directory, copy each subdirectory to `EFT2Obs/cards` and follow the instructions [here](../README.md). The full setup starting from scratch is described in this README.


## Process definition

Create a subdirectory `ttbar-SMEFTsim3` in `EFT2Obs/cards` and add a `proc_card.dat` specifying the model and defining the process:
```
import model SMEFTsim_topU3l_MwScheme_UFO-massless

generate p p > t t~ NP<=1 @0
add process p p > t t~ j SMHLOOP<=1 NP<=1 @1

output ttbar-SMEFTsim3 -nojpeg
```


## Set up MadSpin for top quark decays

We decay the top quarks using `MadSpin`: one leptonically, the other hadronically. To `cards/ttbar-SMEFTsim3`, add a `madspin_card.dat` with the following content:
```
set ms_dir ./madspingrid
set Nevents_for_max_weigth 250
set max_weight_ps_point 400
set max_running_process 1

decay t > w+ b, w+ > l+ vl
decay t~ > w- b~, w- > j j

launch
```


## Set up cards

From the main EFT2Obs directory, run
```sh
./scripts/setup_process.sh ttbar-SMEFTsim3
```
This creates the process directory `MG5_aMC_v2_6_7/ttbar-SMEFTsim3`, and adds two cards to `cards/ttbar-SMEFTsim3`: `pythia8_card.dat` and `run_card.dat`. 

In `pythia8_card.dat`, the line `partonlevel:mpi = off` can be uncommented for now. This disables multiparton interactions, which makes the event generation faster.

In `run_card.dat`, two changes have to be made:
- `True = use_syst`: change to `False = use_syst`
- `systematics = systematics_program`: change to `none = systematics_program`

This is necessary because the extra weights that appear when `True = use_syst` is set are not handled correctly by EFT2Obs.


## Identify EFT operators

To automatically detect the EFT operators this process is sensitive to, and set up the reweighting, run
```sh
python scripts/auto_detect_operators.py -p ttbar-SMEFTsim3
```
This creates two new files in the `cards` subdirectory: `config.json` and `reweight_card.dat`. 

The last file we need to add is `param_card.dat`. To create it, run
```sh
python scripts/make_param_card.py -p ttbar-SMEFTsim3 -c cards/ttbar-SMEFTsim3/config.json \
  -o cards/ttbar-SMEFTsim3/param_card.dat
```


## Parton Matching/Merging

Since we added a process with extra jets, we need to use MLM matching/merging to avoid double counting of parton shower and matrix element events.

To enable MLM merging, in `run_card.dat`, set
- `1 = ickkw`
- `30.0 = xqcut`

For this process, MadGraph does not allow to use MLM merging. As a workaround we can do the following:
- Delete the process directory `MG5_aMC_v2_6_7/ttbar-SMEFTsim3`
- In `proc_card.dat`, change the process definitions:
```
generate p p > t t~ NP=0 @0
add process p p > t t~ j SMHLOOP<=1 NP=0 @1
```
- In `reweight_card.dat`, add these lines between the first and the second line:
```
change process p p > t t~ NP<=1
change process p p > t t~ j SMHLOOP<=1 NP<=1 --add
```
- Then, run again
```sh
./scripts/setup_process.sh ttbar-SMEFTsim3
```
The existing cards in `cards/ttbar-SMEFTsim3` will not be overwritten. At this point, the contents of `cards/ttbar-SMEFTsim3` should be exactly the same as the cards provided in this directory.


## Make gridpack

The gridpack can be created locally with
```sh
./scripts/make_gridpack.sh ttbar-SMEFTsim3 0 16
```
If this takes too long, use the script `launch_gridpack.py`:
```sh
python scripts/launch_gridpack.py ttbar-SMEFTsim3 -c 16 --job-mode slurm
```
The file `gridpack_ttbar-SMEFTsim3.tar.gz` will be copied to the main EFT2Obs directory.


## Generate events

Now we can generate events. This will run through the event generation with `MG5_aMC@NLO`, EFT reweighting, showering with `Pythia`, decaying the top quarks with `MadSpin`, and finally event selection with `Rivet`. Generate 1 million events in a set of slurm jobs:
```sh
python scripts/launch_jobs.py --gridpack gridpack_ttbar-SMEFTsim3.tar.gz -j 200 -s 1 -e 5000 \
  -p CMS_2018_I1663958 -o ttbar-SMEFTsim3 --task-name ttbar --dir jobs --job-mode slurm
```
When using condor, replace `--job-mode slurm` by `--job-mode condor` and add `--sub-opts '+MaxRuntime = 14400\nrequirements = (OpSysAndVer =?= "CentOS7")'`.


## EFT parameterisation

The output of the previous command is 200 YODA files containing all the Rivet routine histograms with a copy for each weight. Merge the yoda files using
```sh
yodamerge -o ttbar/RivetTotal.yoda ttbar/Rivet_* --no-veto-empty
```
Then use the script `get_scaling.py` to produce the JSON file with the EFT scaling parameters $A_{i}$ and $B_{ij}$ (first, copy this file to the main EFT2Obs directory: `eft_exercise_bin_labels.json`):
```sh
python scripts/get_scaling.py -i ttbar-SMEFTsim3/RivetTotal.yoda -o scaling_ttbar-SMEFTsim3 \
  --hist "/CMS_2018_I1663958/d01-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/ttbar-SMEFTsim3/config.json
```
