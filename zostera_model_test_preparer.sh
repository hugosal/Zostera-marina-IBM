#!/bin/bash
echo "Zostera-marina-IBM model preparer"
echo ""
mkdir inputs
mv sample_environment2000.csv sample_founding_rhizomes_2000.csv sample_world.py inputs
cd inputs
python sample_world.py
echo "Succes"
exit 0