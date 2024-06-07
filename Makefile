launch_cluster:
	sbatch launch_cluster.sh

clean:
	@rm -rf .nextflow.log*
	@rm -rf work
	@rm -f logs/*

pull:
	echo "Pulling containers ..."
	@mkdir -p containers
	@for container in `grep -oP "(?<=container = ').*(?=')" confs/slurm.config`; do \
		echo "Checking if $$container is already downloaded ..."; \
		containerName=`echo $$container | sed 's/\//-/g' | sed 's/:/-/g'`; \
		if [ -f containers/$$containerName.img ]; then \
			echo "$$containerName.img already exists. Skipping download."; \
			echo ""; \
		else \
			echo "Pulling $$container ..."; \
			singularity -s pull --name $$containerName.img --dir containers/ docker://$$container; \
			echo ""; \
		fi; \
	done
	@echo "Done!"