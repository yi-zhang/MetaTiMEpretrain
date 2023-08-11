
# 
# echo "[STEP-1 Optional ] Merge GSM dataset into GSE."
#python mergecohort.py --datadir ../mmmetatime202306/data --metafile ../mmmetatime202307cohort/breast_final_downloaded_metatable_gsm.csv --outdir ../mmmetatime202307cohort/gsedata -t 8
echo "[STEP-1] Per dataset preprocessing"
python metatimetrain/scpp.py --datadir ./test/testdata/ --listfile ./test/testdata/datalst1.tsv --outdir ./test/analysis/pp
echo "[STEP-2] Per dataset Decomposition"
python metatimetrain/decompose.py -d ./test/analysis/pp/ -t 4 -o ./test/analysis/decompose/ -k 100
echo "[STEP-3] Align Decomposed components"
python metatimetrain/aligncomp.py -c ./test/analysis/decompose -o ./test/analysis/decompose_align/ -k 100
echo "[STEP-4] Pull together Decomposed components"
python metatimetrain/pullcomp.py -c ./test/analysis/decompose_align/ -o ./test/analysis/decompose_pull/ 
echo "[STEP-5] Call Meta-Components"
python metatimetrain/callmec.py -c ./test/analysis/decompose_pull/ -o ./test/analysis/MeC/ -u True -s 2

