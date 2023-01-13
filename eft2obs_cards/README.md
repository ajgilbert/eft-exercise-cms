# Generating the EFT parameterisations

First set up the SMEFTsim3 model:
```sh
cd EFT2Obs
source env.sh
./scripts/setup_model_SMEFTsim3.sh
```
The cards for each process are in:
- W $\gamma$: `WG-SMEFTsim3`
- single top: `st_tch_4f-SMEFTsim3`
- WW: `WW-SMEFTsim3`
- Zjj: `Zjj-SMEFTsim3`
- H $\rightarrow \gamma\gamma$ (STXS): 
  - qqH: `qqH-SMEFTsim3`
  - ggH: `ggH-SMEFTsim3`
  - WH: `WH-SMEFTsim3`
  - ZH: `ZH-SMEFTsim3`
  - ttH: `ttH-SMEFTsim3`
  - tH: `tH-SMEFTsim3`

Copy each of the card directories in `eft2obs_cards` to `EFT2Obs/cards/` and go to the main EFT2Obs directory. Then set up each process:
```sh
./scripts/setup_process.sh qqH-SMEFTsim3
./scripts/setup_process.sh ggH-SMEFTsim3
./scripts/setup_process.sh WH-SMEFTsim3
./scripts/setup_process.sh ZH-SMEFTsim3
./scripts/setup_process.sh ttH-SMEFTsim3
./scripts/setup_process.sh tH-SMEFTsim3
./scripts/setup_process.sh st_tch_4f-SMEFTsim3
./scripts/setup_process.sh WG-SMEFTsim3
./scripts/setup_process.sh WW-SMEFTsim3
./scripts/setup_process.sh Zjj-SMEFTsim3
```
Create the gridpacks:
```sh
./scripts/make_gridpack.sh qqH-SMEFTsim3 0 8
./scripts/make_gridpack.sh ggH-SMEFTsim3 0 16
./scripts/make_gridpack.sh WH-SMEFTsim3 0 16
./scripts/make_gridpack.sh ZH-SMEFTsim3 0 32
./scripts/make_gridpack.sh ttH-SMEFTsim3 0 64
./scripts/make_gridpack.sh tH-SMEFTsim3 0 32
./scripts/make_gridpack.sh st_tch_4f-SMEFTsim3 0 32
./scripts/make_gridpack.sh WG-SMEFTsim3 0 32
./scripts/make_gridpack.sh WW-SMEFTsim3 0 32
./scripts/make_gridpack.sh Zjj-SMEFTsim3 0 32
```
In some cases, the gridpack generation might take too long to run locally. In that case, the script `launch_gridpack.py` can be used:
```sh
python scripts/launch_gridpack.py [X]-SMEFTsim3 -c 64 --job-mode slurm --task_name [X]_gridpack
```
Generate 1 million events (in processes that are sensitive to a large number of Wilson coefficients, the reweighting can take along time, so the number of events per job should be reduced to make the event generation faster):
```sh
python scripts/launch_jobs.py --gridpack gridpack_qqH-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p HiggsTemplateCrossSections -o qqH-SMEFTsim3 --task-name qqH --dir jobs --job-mode slurm \
  --env "HIGGSPRODMODE=VBF"
python scripts/launch_jobs.py --gridpack gridpack_ggH-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p HiggsTemplateCrossSections -o ggH-SMEFTsim3 --task-name ggH --dir jobs --job-mode slurm \
  --env "HIGGSPRODMODE=GGF"
python scripts/launch_jobs.py --gridpack gridpack_WH-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p HiggsTemplateCrossSections -o WH-SMEFTsim3 --task-name WH --dir jobs --job-mode slurm \
  --env "HIGGSPRODMODE=WH"
python scripts/launch_jobs.py --gridpack gridpack_ZH-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p HiggsTemplateCrossSections -o ZH-SMEFTsim3 --task-name ZH --dir jobs --job-mode slurm \
  --env "HIGGSPRODMODE=QQ2ZH"
python scripts/launch_jobs.py --gridpack gridpack_ttH-SMEFTsim3.tar.gz -j 200 -s 1 -e 5000 \
  -p HiggsTemplateCrossSections -o ttH-SMEFTsim3 --task-name ttH --dir jobs --job-mode slurm \
  --env "HIGGSPRODMODE=TTH"
python scripts/launch_jobs.py --gridpack gridpack_tH-SMEFTsim3.tar.gz -j 100 -s 1 -e 10000 \
  -p HiggsTemplateCrossSections -o tH-SMEFTsim3 --task-name tH --dir jobs --job-mode slurm \
  --env "HIGGSPRODMODE=TH"
python scripts/launch_jobs.py --gridpack gridpack_st_tch_4f-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p CMS_2019_I1744604 -o st_tch_4f-SMEFTsim3 --task-name st_tch_4f --dir jobs --job-mode slurm
python scripts/launch_jobs.py --gridpack gridpack_WG-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p CMS_2021_PAS_SMP_20_005 -o WG-SMEFTsim3 --task-name WG --dir jobs --job-mode slurm
python scripts/launch_jobs.py --gridpack gridpack_WW-SMEFTsim3.tar.gz -j 50 -s 1 -e 20000 \
  -p CMS_2021_PAS_SMP_20_005 -o WW-SMEFTsim3 --task-name WW --dir jobs --job-mode slurm
python scripts/launch_jobs.py --gridpack gridpack_Zjj-SMEFTsim3.tar.gz -j 100 -s 1 -e 10000 \
  -p ATLAS_2020_I1803608:TYPE=EW_ONLY -o Zjj-SMEFTsim3 --task-name Zjj --dir jobs --job-mode slurm
```
When using condor, replace `--job-mode slurm` by `--job-mode condor` and add `--sub-opts '+MaxRuntime = 14400\nrequirements = (OpSysAndVer =?= "CentOS7")'`.

Merge the output yoda files using
```sh
yodamerge -o RivetTotal.yoda Rivet_* --no-veto-empty
```
Produce the JSON files with the EFT scaling parameters $A_{i}$ and $B_{ij}$:
```sh
python scripts/get_scaling.py -i qqH-SMEFTsim3/RivetTotal.yoda -o scaling_qqH-SMEFTsim3 \
  --hist "/HiggsTemplateCrossSections/HTXS_stage1_2_pTjet30" --bin-labels eft_exercise_bin_labels.json \
  -c cards/qqH-SMEFTsim3/config.json --rebin 18,19,20,21,22,23,24,25,26,27,28,29
python scripts/get_scaling.py -i ggH-SMEFTsim3/RivetTotal.yoda -o scaling_ggH-SMEFTsim3 \
  --hist "/HiggsTemplateCrossSections/HTXS_stage1_2_pTjet30" --bin-labels bin_labels_ggh.json \
  -c cards/ggH-SMEFTsim3/config.json --rebin 2,3,4,6,7,8,9,10,11,12,13
python scripts/get_scaling.py -i WH-SMEFTsim3/RivetTotal.yoda -o scaling_WH-SMEFTsim3 \
  --hist "/HiggsTemplateCrossSections/HTXS_stage1_2_pTjet30" --bin-labels bin_labels_wh.json \
  -c cards/WH-SMEFTsim3/config.json --rebin 30,31,32,35
python scripts/get_scaling.py -i ZH-SMEFTsim3/RivetTotal.yoda -o scaling_ZH-SMEFTsim3 \
  --hist "/HiggsTemplateCrossSections/HTXS_stage1_2_pTjet30" --bin-labels bin_labels_zh.json \
  -c cards/ZH-SMEFTsim3/config.json --rebin 36,41
python scripts/get_scaling.py -i ttH-SMEFTsim3/RivetTotal.yoda -o scaling_ttH-SMEFTsim3 \
  --hist "/HiggsTemplateCrossSections/HTXS_stage1_2_pTjet30" --bin-labels bin_labels_tth.json \
  -c cards/ttH-SMEFTsim3/config.json --rebin 48,49,50,51,52,53
python scripts/get_scaling.py -i tH-SMEFTsim3/RivetTotal.yoda -o scaling_tH-SMEFTsim3 \
  --hist "/HiggsTemplateCrossSections/HTXS_stage1_2_pTjet30" --bin-labels bin_labels_th.json \
  -c cards/tH-SMEFTsim3/config.json --rebin 56,57
python scripts/get_scaling.py -i st_tch_4f-SMEFTsim3/RivetTotal.yoda -o scaling_st_tch_4f-SMEFTsim3 \
  --hist "/CMS_2019_I1744604/d13-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/st_tch_4f-SMEFTsim3/config.json
python scripts/get_scaling.py -i WG-SMEFTsim3/RivetTotal.yoda -o scaling_WG-SMEFTsim3-d54-x01-y01 \
  --hist "/CMS_2021_PAS_SMP_20_005/d54-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/WG-SMEFTsim3/config.json
python scripts/get_scaling.py -i WG-SMEFTsim3/RivetTotal.yoda -o scaling_WG-SMEFTsim3-d55-x01-y01 \
  --hist "/CMS_2021_PAS_SMP_20_005/d55-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/WG-SMEFTsim3/config.json
python scripts/get_scaling.py -i WG-SMEFTsim3/RivetTotal.yoda -o scaling_WG-SMEFTsim3-d56-x01-y01 \
  --hist "/CMS_2021_PAS_SMP_20_005/d56-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/WG-SMEFTsim3/config.json
python scripts/get_scaling.py -i WW-SMEFTsim3/RivetTotal.yoda -o scaling_WW-SMEFTsim3 \
  --hist "/ATLAS_2019_I1734263/d04-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/WW-SMEFTsim3/config.json
python scripts/get_scaling.py -i Zjj-SMEFTsim3/RivetTotal.yoda -o scaling_Zjj-SMEFTsim3 \
  --hist "/ATLAS_2020_I1803608:TYPE=EW_ONLY/d04-x01-y01" --bin-labels eft_exercise_bin_labels.json \
  -c cards/Zjj-SMEFTsim3/config.json
```

