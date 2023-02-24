# EFT parameterisations for CMS-TOP-17-002

To generate the EFT parameterisations for the ttbar measurement using the cards provided in this directory, copy each subdirectory to `EFT2Obs/cards` and follow the instructions [here](../README.md). The full setup starting from scratch is described in this README.


## Process definition

Create two subdirectories `ttbar-tlep-SMEFTsim3` and `ttbar-tbarlep-SMEFTsim3` in `EFT2Obs/cards` and add a `proc_card.dat` to both, specifying the model and defining the process:
```
import model SMEFTsim_topU3l_MwScheme_UFO-massless

generate p p > t t~ NP<=1 @0
add process p p > t t~ j SMHLOOP<=1 NP<=1 @1

output ttbar-{tlep/tbarlep}-SMEFTsim3 -nojpeg
```


## Set up MadSpin for top quark decays

We decay the top quarks using `MadSpin`: one leptonically, the other hadronically. To `cards/ttbar-tlep-SMEFTsim3`, add a `madspin_card.dat` with the following content:
```
set ms_dir ./madspingrid
set Nevents_for_max_weigth 250
set max_weight_ps_point 400
set max_running_process 1

decay t > w+ b, w+ > l+ vl
decay t~ > w- b~, w- > j j

launch
```
We also need to consider the opposite decay: To `cards/ttbar-tbarlep-SMEFTsim3`, add a `madspin_card.dat` with
```
decay t > w+ b, w+ > j j
decay t~ > w- b~, w- > l- vl~
```


## Set up cards

From the main EFT2Obs directory, run
```sh
./scripts/setup_process.sh ttbar-tlep-SMEFTsim3
./scripts/setup_process.sh ttbar-tbarlep-SMEFTsim3
```
This creates the process directories `MG5_aMC_v2_6_7/ttbar-tlep-SMEFTsim3` and `MG5_aMC_v2_6_7/ttbar-tbarlep-SMEFTsim3`, and adds two new cards to `cards/ttbar-tlep-SMEFTsim3` and `cards/ttbar-tbarlep-SMEFTsim3`: `pythia8_card.dat` and `run_card.dat`. 

In `pythia8_card.dat`, the line `partonlevel:mpi = off` can be uncommented for now. This disables multiparton interactions, which makes the event generation faster.

In `run_card.dat`, two changes have to be made:
- `True = use_syst`: change to `False = use_syst`
- `systematics = systematics_program`: change to `none = systematics_program`

This is necessary because the extra weights that appear when `True = use_syst` is set are not handled correctly by EFT2Obs.


## Identify EFT operators

To automatically detect the EFT operators this process is sensitive to, and set up the reweighting, run
```sh
python scripts/auto_detect_operators.py -p ttbar-tlep-SMEFTsim3
python scripts/auto_detect_operators.py -p ttbar-tbarlep-SMEFTsim3
```
This creates two new files in the `cards` subdirectories: `config.json` and `reweight_card.dat`. 

The last file we need to add is `param_card.dat`. To create it, run
```sh
python scripts/make_param_card.py -p ttbar-tlep-SMEFTsim3 -c cards/ttbar-tlep-SMEFTsim3/config.json \
  -o cards/ttbar-tlep-SMEFTsim3/param_card.dat
python scripts/make_param_card.py -p ttbar-tbarlep-SMEFTsim3 -c cards/ttbar-tbarlep-SMEFTsim3/config.json \
  -o cards/ttbar-tbarlep-SMEFTsim3/param_card.dat
```


## Parton Matching/Merging

Since we added a process with extra jets, we need to use MLM matching/merging to avoid double counting of parton shower and matrix element events.

To enable MLM merging, in `run_card.dat`, set
- `1 = ickkw`
- `30.0 = xqcut`

For this process, MadGraph does not allow to use MLM merging. As a workaround we can do the following:
- Delete the process directory `MG5_aMC_v2_6_7/ttbar-tlep-SMEFTsim3` and `MG5_aMC_v2_6_7/ttbar-tbarlep-SMEFTsim3`
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
./scripts/setup_process.sh ttbar-tlep-SMEFTsim3
./scripts/setup_process.sh ttbar-tbarlep-SMEFTsim3
```
The existing cards in `cards/ttbar-tlep-SMEFTsim3` and `cards/ttbar-tbarlep-SMEFTsim3` will not be overwritten. At this point, the contents of the two `cards` subdirectories should be exactly the same as the cards provided in this directory.


## Make gridpack

The gridpacks can be created locally with
```sh
./scripts/make_gridpack.sh ttbar-tlep-SMEFTsim3 0 16
./scripts/make_gridpack.sh ttbar-tbarlep-SMEFTsim3 0 16
```
If this takes too long, use the script `launch_gridpack.py`:
```sh
python scripts/launch_gridpack.py ttbar-tlep-SMEFTsim3 -c 16 --job-mode slurm
python scripts/launch_gridpack.py ttbar-tbarlep-SMEFTsim3 -c 16 --job-mode slurm
```
The files `gridpack_ttbar-tlep-SMEFTsim3.tar.gz` and `gridpack_ttbar-tbarlep-SMEFTsim3.tar.gz` will be copied to the main EFT2Obs directory.


## Generate events

Now we can generate events. This will run through the event generation with `MG5_aMC@NLO`, EFT reweighting, showering with `Pythia`, decaying the top quarks with `MadSpin`, and finally event selection with `Rivet`. The Rivet routine for this analysis contains a large number of histograms, and there are many reweighting points, due to the large number of Wilson coefficients this process is sensitive too. To keep the size of the output YODA files manageable, you can edit the Rivet routine that is saved in `EFT2Obs/Rivet-3.0.1/analyses/pluginCMS/CMS_2018_I1663958.cc` and remove the histograms that aren't needed. 

Now you can generate 1 million events for each decay mode in a set of slurm jobs:
```sh
python scripts/launch_jobs.py --gridpack gridpack_ttbar-tlep-SMEFTsim3.tar.gz -j 200 -s 1 -e 5000 \
  -p CMS_2018_I1663958 -o ttbar-tlep-SMEFTsim3 --task-name ttbar --dir jobs --job-mode slurm
python scripts/launch_jobs.py --gridpack gridpack_ttbar-tbarlep-SMEFTsim3.tar.gz -j 200 -s 1 -e 5000 \
  -p CMS_2018_I1663958 -o ttbar-tbarlep-SMEFTsim3 --task-name ttbar --dir jobs --job-mode slurm
```
When using condor, replace `--job-mode slurm` by `--job-mode condor` and add `--sub-opts '+MaxRuntime = 14400\nrequirements = (OpSysAndVer =?= "CentOS7")'`.


## EFT parameterisation

The output of the previous command is 200 YODA files per decay mode, containing all the Rivet routine histograms with a copy for each weight. Merge the yoda files using
```sh
yodamerge -o ttbar-tlep-SMEFTsim3/RivetTotal.yoda ttbar-tlep-SMEFTsim3/Rivet_* --no-veto-empty
yodamerge -o ttbar-tbarlep-SMEFTsim3/RivetTotal.yoda ttbar-tbarlep-SMEFTsim3/Rivet_* --no-veto-empty
```
If you have used the unedited Rivet routine, each of these files will be around to 400 MB, so you might need to merge them in batches of a few files at a time. Next, add the histograms from the two ttbar decay modes using `yodamerge` with the option `--add`:
```sh
mkdir ttbar-merged
yodamerge -o ttbar-merged/RivetTotal.yoda ttbar-tlep-SMEFTsim3/RivetTotal.yoda \
  ttbar-tbarlep-SMEFTsim3/RivetTotal.yoda --no-veto-empty --add
```
Then use the script `get_scaling.py` to produce the JSON file with the EFT scaling parameters $A_{i}$ and $B_{ij}$ (first, copy this file to the main EFT2Obs directory: `eft_exercise_bin_labels.json`):
```sh
python scripts/get_scaling.py -i ttbar-merged/RivetTotal.yoda -o scaling_ttbar-SMEFTsim3 \
  --hist "/CMS_2018_I1663958/d05-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/ttbar-tlep-SMEFTsim3/config.json
```
