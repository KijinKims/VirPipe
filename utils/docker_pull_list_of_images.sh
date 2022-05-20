#!/bin/bash
input=$1

grep -v '^#' < $input | { while read line; do docker pull $line; done; }
