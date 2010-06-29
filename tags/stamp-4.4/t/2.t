#!/bin/bash

../bin/alignfit -f $test.blc -d $test.dom -out $tmp$test.alignfit
diff -I '^%' --brief $test.alignfit $tmp$test.alignfit
