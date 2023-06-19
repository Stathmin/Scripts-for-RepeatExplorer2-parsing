#for i in $(ls *.fasta); do ./muscle5.1.linux_intel64 -super5 $i -output $i.ada || echo $i; done

for i in $(ls *.fasta); do (mafft --localpair --maxiterate 1000 --adjustdirection $i > $i.maf ) || echo $i; done
