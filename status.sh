if [ -s logs/snakemake.log ]; then
   
	python parse_snakemake_log.py logs/snakemake.log

else 

	echo "logs/snakemake.log is empty"

fi
