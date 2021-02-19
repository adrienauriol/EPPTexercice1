#!/bin/bash

for number in {1..5}
do
python analysis.py case${number} momentumMode
python analysis.py case${number}_PbBefore momentumMode
python analysis.py case${number}_PbAfter momentumMode
done

python analysis.py case1_protons calorimetryMode
python analysis.py case1_positrons calorimetryMode


exit 0