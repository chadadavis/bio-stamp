#!/bin/bash

${bindir}pdbseq -f $test.dom -tl 30 > $tmp$test.aa
diff --brief $test.aa $tmp$test.aa
