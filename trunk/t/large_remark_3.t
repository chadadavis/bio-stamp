#!/bin/bash

# This file has many REMARK  3 lines
# STAMP previously crashed before reading them all and did not output 'refinement'
$srcdir/pdbc -d 2jkv | grep -q refinement
