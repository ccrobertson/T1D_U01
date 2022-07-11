pipeline=$1
if [ -s logs/${pipeline}/snakemake.log ]; then
   
	python parse_snakemake_log.py logs/${pipeline}_snakemake.log

else 

	echo "logs/${pipeline}_snakemake.log is empty"

fi
