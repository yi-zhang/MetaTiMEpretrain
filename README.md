
# Pipeline for calling Meta-Components from multiple scRNA-seq data.

This MetaTiMEpretrain repo, when ran on a large set of scRNA-seq samples, can generate gene programs corresponding to cell types, cell states, or signaling pathways in the provided data. Meta-components called from scRNA-seq represent independent transcriptional variations, reproducibly seen in data, in format as weighted gene contribution vectors. For now the low-dim method is independent component analysis. Check [MetaTiME paper](https://www.nature.com/articles/s41467-023-38333-8) or [MetaTiME annotator repo](https://github.com/yi-zhang/MetaTiME) for the scenario of Tumor Microenvironment.
[in progress :)]

## Dependency

- scanpy, pandas, multiprocessing, sklearn, seaborn, scikitnetwork >=0.28.2

## Test run

- `git clone https://github.com/yi-zhang/MetaTiME.git`
- `cd MetaTiMEpretrain`
- Collect your scRNA datasets in  <input dir> as h5ad format. Or, a few tumor scRNA datasets are provided in this [sample input link](https://www.dropbox.com/scl/fo/udl7ep9juxqn79bj64vlm/h?rlkey=jzg5m9zqxmnl7ec5iaj1l8cqu&dl=0). Download folder and cp to `MetaTiMEpretrain/test/`, or point to it as input `datadir` in `scpp.py`.
- `sh test.sh` # This is a script containing a few sequential steps, explained below.
- An output folder for this sample input data is available from [result dir for testrun](https://www.dropbox.com/scl/fo/udl7ep9juxqn79bj64vlm/h?rlkey=jzg5m9zqxmnl7ec5iaj1l8cqu&dl=0).

## Steps in pipeline

Each steps can be ran separately.

### Step 1: Parallel single dataset preprocessing

- `python metatimetrain/scpp.py --datadir ./test/testdata/ --listfile ./test/testdata/datalst1.tsv --outdir ./test/analysis/pp`
- `--datadir`: Input directory with scRNA files.
- `--outdir` : output directory with preprocessed scRNA files. Counts will be depth normalized and log transformed.
- Ran in parallel. Can take longer time (~hr) and big memory depends on data size.

### Step 2: Per dataset decomposition.

- `python metatimetrain/decompose.py -d ./test/analysis/pp/ -t 4 -o ./test/analysis/decompose/ -k 100`
- `-d` : Input directory with preprocessed scRNA files.
- `-t`: Number of threads.
- `-o`: Output directory of per-dataset decomposition table.
- `-k`: Number of low-dim components.
- Ran in parallel. Can take longer time (~hr) and big memory depends on data size.

### STEP-3: Align low-dim components

- `python metatimetrain/aligncomp.py -c ./test/analysis/decompose -o ./test/analysis/decompose_align/ -k 100`
- `-c` : Input directory with per-dataset decomposition table.
- `-o`: Output directory of aligned decomposition tables.
- `-k`: Number of low-dim components. Must be the same as in input decomposition tables.

### STEP-4: Pull low-dim components

- `python metatimetrain/pullcomp.py -c ./test/analysis/decompose_align/ -o ./test/analysis/decompose_pull/`
- `-c` : Input directory with aligned, per-dataset decomposition table.
- `-o`: Output directory to store a table gathering all low-dim components

### STEP-5: Call Meta-Components!

- `python metatimetrain/callmec.py -c ./test/analysis/decompose_pull/ -o ./test/analysis/MeC/ -u True -s 2`
- `-c` : Input directory with pulled components.
- `-u` : Unit test or not. `True`  for test datasets
- `-s` : Minimum number of components per meta-component cluster.

---

### [optional] Step-0: merge scRNA samples by cohort.

- `python [mergecohort.py](http://mergecohort.py/) --datadir ../mmmetatime202306/data --metafile ../mmmetatime202307cohort/breast_final_downloaded_metatable_gsm.csv --outdir ../mmmetatime202307cohort/gsedata -t 8`


