#!/bin/sh
################################################################################
##
## FILE: extract_all
##
## DESCRIPTION: extracts the body from each of the SGML files listed, in turn
##
## SYNOPSIS: extract_all
##
################################################################################

FILES="reut2-001.sgm reut2-002.sgm reut2-003.sgm reut2-004.sgm reut2-005.sgm reut2-006.sgm reut2-007.sgm reut2-008.sgm reut2-009.sgm reut2-010.sgm reut2-011.sgm reut2-012.sgm reut2-013.sgm reut2-014.sgm reut2-015.sgm reut2-016.sgm reut2-017.sgm reut2-018.sgm reut2-019.sgm reut2-020.sgm reut2-021.sgm"

for file in ${FILES}
do
    echo creating ${file:6:3}.body.txt
    ./get_reuters_body.sh < $file > ${file:6:3}.body.txt
done

