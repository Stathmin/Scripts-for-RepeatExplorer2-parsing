#required folders:
#./comparatives_db, ./local_db_solo, ./ncbi_repeats_db,
#./solo_annotation

#./comparatives_db, ./local_db_solo and the current folder are the
#only folders accessed and modified by the given scripts
#ncbi_repeats_db stores note.txt file with an example of NCBI
#Nucleotide query we used to build this db

#required files:
#./ncbi_naming.csv, ./cypher.csv, filled as shown

#In our study we had to work with multiple RepeatExplorer2
#assemblies, done on multiple servers. To avoid transfering numerous
#heavy files, we have produced a script to selectively copy only the
#most important files of an assembly and later performed analysis
#with them
#e.g.:
#screenshot_repex.py foo/bar/KA16Assembly/RepExRES ./KA16
#~10-20Mb of report tables, tarean images and sequences are copied to
#the selected folder

#after all the needed folders for the experiment are created in ./,
#we run:

bash get_fastas.sh

#which creates dbs for further blast from TAREAN fastas

python3 ./represenation.py 'KA27' 'local'
#creates report in ./solo_annotation for single species

python3 ./represenation.py 'KA25KK6(1to2)' 'comp'
#creates report in ./solo_annotation for comparative assembly

python3 repeat_sums_by_group.py 'KA27'
#creates in ./solo_annotation the table and, if plotly library is
#installed, sunburst plot html and png