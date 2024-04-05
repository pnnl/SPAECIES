#!/bin/sh

make -j
for i in {0..11}
do
    ./rainshaft i
done