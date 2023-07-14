# Scripts-for-RepeatExplorer2-parsing
![image](https://github.com/Stathmin/Scripts-for-RepeatExplorer2-parsing/assets/55657873/c23d7f66-8c10-4dcd-b468-b2a72771fdea)

./screenshot_repex.py - script used to copy the most infromative lightweight files from RepeatExplorer2 output

./KA16, ./KA27, ./KA25KK6(1to2) - examples of RepeatExplorer2 output, copied with screenshot_repex.py
./cypher.csv - an example of species annotation used for table creation
./README.txt - demo analysis

./ncbi_repeats_db - the local variant of NCBI Nucleotide we have used with notes on how to adapt it for your experiment
./ncbi_naming.csv - an example of NCBI seqs annotation
./solo_annotation - an output directory for the most scripts

./get_fastas.sh - used to create blastdbs for blast.py, stores them in ./local_db_solo and ./comparatives_db
e.g.:
bash ./get_fastas.sh

./blast.py - required for represenation.py

./represenation.py - produces tables for the RepeatExplorer2 output
e.g.: 
python3 ./represenation.py "KA16" "local"
python3 ./represenation.py "KA25KK6(1to2)" "comp"

./repeat_sums_by_group.py - used to create tables with percentages of mobile elements by type
