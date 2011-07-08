#!/bin/bash

rm $tmp$test.alignfit.out 2>/dev/null
$srcdir/alignfit -f $test.blc -d $test.dom -out $tmp$test.alignfit &> $tmp$test.alignfit.out
diff -I '^%' $test.alignfit $tmp$test.alignfit && rm $tmp$test.*
