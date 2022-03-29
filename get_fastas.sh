rm -rf local_db_solo comparatives_db

mkdir local_db_solo
mkdir comparatives_db
touch local_db_solo/multifasta.fasta
touch comparatives_db/COMPBASE.fasta

a=$(ls ./*/*TAREAN*.fasta | sed 's#/TAREAN.*##g' | sort | uniq)
b=$(ls ./*/COMP* | sed 's#/COMP.*##g')
c=$(comm -23 <(echo $a| sed 's# #\n#g'|sort) <(echo $b| sed 's# #\n#g'|sort))

for j in `echo $c`; do k=$(echo $j|sed 's#\./##'); cat $j/*TAREAN* | sed "s#>#>${k}_#g" >> local_db_solo/multifasta.fasta; done
for j in `echo $b`; do k=$(echo $j|sed 's#\./##'); cat $j/*TAREAN* | sed "s#>#>${k}_#g" >> comparatives_db/COMPBASE.fasta; done

cat 'local_db_solo/multifasta.fasta'| sed -Ez 's#\n##g'|sed -E 's#([ATCG]{2,})#\n\1\1\1\n#g' > local_db_solo/multifasta_x3.fasta
cd local_db_solo
makeblastdb -dbtype nucl -in multifasta.fasta
makeblastdb -dbtype nucl -in multifasta_x3.fasta
cd ../

cat 'comparatives_db/COMPBASE.fasta'| sed -Ez 's#\n##g'|sed -E 's#([ATCG]{2,})#\n\1\1\1\n#g' > comparatives_db/COMPBASE_x3.fasta
cd comparatives_db
makeblastdb -dbtype nucl -in COMPBASE.fasta
makeblastdb -dbtype nucl -in COMPBASE_x3.fasta
cd ../