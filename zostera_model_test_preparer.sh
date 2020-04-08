#!/bin/bash
echo "Zostera-marina-IBM model preparer"
if [ ! -d inputs/ ]; then 
	mkdir inputs
fi
if [[ ! -f inputs/sample_environment_2000.csv && ! -f inputs/sample_founding_rhizomes_2000.csv && ! -f inputs/sample_world.py ]]; then
	mv sample_environment_2000.csv sample_founding_rhizomes_2000.csv sample_world.py inputs
	cd inputs
	python sample_world.py
fi
