snakemake all -n --jobname "islets.{jobid}" --jobs 100 \
		--keep-going \
		--rerun-incomplete \
		--snakefile workflow/src/multiome.smk \
		--configfile workflow/src/multiome.yaml \
		--use-conda \
		--printshellcmds \
		--cluster-config workflow/envs/cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus}" \
		> logs/snakemake.log 2>&1 &
