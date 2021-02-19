#!/bin/bash

for number in {7..8}
do
./../build/exampleB5 mymacro${number}.mac | tee mymacro${number}.out
done

exit 0
