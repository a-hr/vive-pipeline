launch_cluster:
	sbatch launch_cluster.sh

clean:
	rm -rf .nextflow.log*
	rm -rf work
	rm -f logs/*
	clear