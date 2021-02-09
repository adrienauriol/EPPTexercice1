#!/bin/bash

for number in {1..5}
do 
python analysis.py case${number}
done

exit 0