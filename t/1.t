#!/bin/bash

rm $tmp$test.aa 2>/dev/null
$srcdir/pdbseq -f $test.dom -tl 30 > $tmp$test.aa
diff --brief $test.aa $tmp$test.aa && rm $tmp$test.*
