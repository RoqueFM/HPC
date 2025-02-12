#! /bin/bash

for i in run_??_*sh; do
	echo -- "$i"
	./$i || exit 1;
done
