########
python sumstats.py csv --sumstats BRCA_GCST010098.txt.gz --out BRCA_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats BRCA_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > BRCA_EUR_qc_noMHC.tsv.gz

python sumstats.py csv --sumstats CRC_GCST90255675.txt.gz --out CRC_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats CRC_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > CRC_EUR_qc_noMHC.tsv.gz

python sumstats.py csv --sumstats EC_GCST006464.txt.gz --out EC_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats EC_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > EC_EUR_qc_noMHC.tsv.gz

python sumstats.py csv --sumstats LC_GCST004748.txt.gz --out LC_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats LC_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > LC_EUR_qc_noMHC.tsv.gz

python sumstats.py csv --sumstats OC_GCST004462.txt.gz --out OC_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats OC_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > OC_EUR_qc_noMHC.tsv.gz

python sumstats.py csv --sumstats PRCA_GCST90274714.txt.gz --out PRCA_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats PRCA_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > PRCA_EUR_qc_noMHC.tsv.gz

python sumstats.py csv --sumstats RCC_GCST90320057.txt.gz --out RCC_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats RCC_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > RCC_EUR_qc_noMHC.tsv.gz

python sumstats.py csv --sumstats THCA_GCST90399736.txt.gz --out THCA_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats THCA_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > THCA_EUR_qc_noMHC.tsv.gz


python sumstats.py csv --sumstats ESCA_GCST003739.txt.gz --out ESCA_EUR.csv --force --auto --head 5
python sumstats.py zscore --sumstats ESCA_EUR.csv \
| python sumstats.py qc \
--sumstats - \
--exclude-ranges 6:26000000-34000000 --max-or 1e37 \
| gzip > ESCA_EUR_qc_noMHC.tsv.gz

cd /mnt/d/linux
wget https://archives.boost.io/release/1.69.0/source/boost_1_69_0.tar.gz
tar -xzvf boost_1_69_0.tar.gz && cd boost_1_69_0
./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time
./b2 --clean && ./b2 --j12 -a


cd /mnt/d/linux
git clone --recurse-submodules -j8 https://github.com/precimed/mixer.git
mkdir mixer/src/build && cd mixer/src/build
cmake .. -DBOOST_ROOT=/mnt/d/linux/boost_1_69_0 && make bgmg -j16  

cd /mnt/d/linux/mixer
MIXER=/mnt/d/linux/mixer/precimed/mixer.py
LIBBGMG=/mnt/d/linux/mixer/src/build/lib/libbgmg.so
LDSR=/mnt/d/1kg.v3

for i in {1..22}; do
echo "Processing chr$i"
python $MIXER ld \
--lib $LIBBGMG \
--bfile $LDSR/UK10K_1KG_qc_rsid_$i \
--out $LDSR/UK10K_1KG_qc_rsid_${i}.run4.ld \
--r2min 0.05 --ldscore-r2min 0.05 --ld-window-kb 30000
done

for i in {1..20}; do
echo "Pruning chr$i"
python $MIXER snps \
--lib $LIBBGMG \
--bim-file $LDSR/UK10K_1KG_qc_rsid_@.bim \
--ld-file $LDSR/UK10K_1KG_qc_rsid_@.run4.ld \
--out $LDSR/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--maf 0.05 --subset 2000000 --r2 0.8 --seed $i
done

cd /mnt/d/linux/mixer

WKDIR=/mnt/d/linux/mixer
MIXER=/mnt/d/linux/mixer/precimed/mixer.py
LIBBGMG=/mnt/d/linux/mixer/src/build/lib/libbgmg.so
LDSR=/mnt/d/1kg.v3


ID=BRCA
mkdir -p ${WKDIR}/BRCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/BRCA/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=PRCA
mkdir -p ${WKDIR}/PRCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/PRCA/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  

ID=CRC
mkdir -p ${WKDIR}/CRC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/CRC/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  

ID=LC
mkdir -p ${WKDIR}/LC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/LC/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  

ID=EC
mkdir -p ${WKDIR}/EC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/EC/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  

ID=RCC
mkdir -p ${WKDIR}/RCC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/RCC/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  

ID=OC
mkdir -p ${WKDIR}/OC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/OC/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  

ID=THCA
mkdir -p ${WKDIR}/THCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/THCA/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  

ID=ESCA
mkdir -p ${WKDIR}/ESCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} fit1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--out ${WKDIR}/THCA/${ID}.fit.rep${i} \
--extract ${LDSR}/UK10K_1KG_qc_rsid_prune_maf0p05_rand2M_r2p8.rep${i}.snps \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done  


ID=BRCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/BRCA/${ID}.fit.rep${i}.json \
--out ${WKDIR}/BRCA/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=PRCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/PRCA/${ID}.fit.rep${i}.json \
--out ${WKDIR}/PRCA/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=CRC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/CRC/${ID}.fit.rep${i}.json \
--out ${WKDIR}/CRC/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=LC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/LC/${ID}.fit.rep${i}.json \
--out ${WKDIR}/LC/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=EC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/EC/${ID}.fit.rep${i}.json \
--out ${WKDIR}/EC/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=RCC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/RCC/${ID}.fit.rep${i}.json \
--out ${WKDIR}/RCC/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=OC
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/OC/${ID}.fit.rep${i}.json \
--out ${WKDIR}/OC/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=THCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/THCA/${ID}.fit.rep${i}.json \
--out ${WKDIR}/THCA/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID=ESCA
for i in {1..20}
do
echo "Processing step $i"
python ${MIXER} test1 \
--trait1-file ${WKDIR}/${ID}_EUR_qc_noMHC.tsv.gz \
--load-params-file ${WKDIR}/THCA/${ID}.fit.rep${i}.json \
--out ${WKDIR}/THCA/${ID}.test.rep${i} \
--bim-file ${LDSR}/UK10K_1KG_qc_rsid_@.bim \
--ld-file ${LDSR}/UK10K_1KG_qc_rsid_@.run4.ld \
--lib ${LIBBGMG}
done

ID1=THCA
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=LUNG
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=BRCA
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=PRCA
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=CRC
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=ENDO
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=RCC
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=OVCA
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std
ID1=ESCA
python ${MFIGURE} combine --json ${WKDIR}/${ID1}/${ID1}.fit.rep@.json  --out ${WKDIR}/${ID1}/${ID1}.fit
python ${MFIGURE} one --json ${WKDIR}/${ID1}/${ID1}.fit.json --out ${WKDIR}/${ID1}/${ID1} --statistic mean std

./run_fit2_parallel.sh
./run_test2_parallel.sh
ID1=PRCA
ID2=BRCA
ID=${ID1}_${ID2}
WKDIR=/mnt/g/Linux/mixer/${ID}
MFIGURE=/mnt/g/Linux/mixer/precimed/mixer_figures.py
python ${MFIGURE} combine --json ${WKDIR}/${ID}.test.rep@.json --out ${WKDIR}/${ID}.test
python ${MFIGURE} combine --json ${WKDIR}/${ID}.fit.rep@.json  --out ${WKDIR}/${ID}.fit
python ${MFIGURE} two --json-fit ${WKDIR}/${ID}.fit.json --json-test ${WKDIR}/${ID}.test.json --out ${WKDIR}/${ID} --statistic mean std
