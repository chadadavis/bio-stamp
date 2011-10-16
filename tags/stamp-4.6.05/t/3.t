#!/bin/bash

rm $tmp$test.stamp.out 2>/dev/null
$srcdir/stamp -f $test.dom -rough -prefix $tmp$test.stamp &> $tmp$test.stamp.out
diff $test.stamp.out $tmp$test.stamp.out && rm $tmp$test.*
