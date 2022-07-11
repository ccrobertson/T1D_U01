pipeline=$1
snakemake all -n --jobname "${pipeline}.{jobid}" --jobs 100 \
		--keep-going \
		--rerun-incomplete \
		--snakefile workflow/src/${pipeline}.smk \
		--configfile workflow/src/${pipeline}.yaml \
		--use-conda \
		--use-singularity \
		--printshellcmds \
		--cluster-config workflow/envs/cluster_${pipeline}.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus}" \
		> logs/${pipeline}_snakemake.log 2>&1 &
