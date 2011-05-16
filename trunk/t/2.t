#!/bin/bash

${bindir}alignfit -f $test.blc -d $test.dom -out $tmp$test.alignfit &> $tmp$test.alignfit.out
diff -I '^%' --brief $test.alignfit $tmp$test.alignfit
