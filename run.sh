#!/bin/bash
mkdir data
python punta_banda_world.py
mv founding_rhizomes_2000.csv founding_rhizomes_2018.csv environment_2000.csv environment_2018.csv cannal_200m_broad_4m_prof.dat data
python zostera_model.py founding_rhizomes_2000.csv environment_2000.csv cannal_200m_broad_4m_prof.dat False
