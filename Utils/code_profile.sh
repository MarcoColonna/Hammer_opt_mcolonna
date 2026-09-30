#!/bin/bash

rm -R *.trace

for i in `find . -type f -perm +ugo+x -maxdepth 1 | sort -t '\0' -n`
do
    xctrace record \
        --template 'Time Profiler' \
        --output $i.trace \
        --launch -- ./$i

    xctrace export --input $i.trace \
        --xpath '/trace-toc/run[@number="1"]/data/table[@schema="time-profile"]' \
        --output $i.xml

 done

