#!/bin/bash

${bindir}stamp -f $test.dom -rough -prefix $tmp$test.stamp &> $tmp$test.stamp.out
diff --brief $test.stamp.out $tmp$test.stamp.out
