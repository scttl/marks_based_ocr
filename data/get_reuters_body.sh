#!/bin/sh
################################################################################
##
## FILE: get_reuters_body.sh 
##
## DESCRIPTION: strips out the main body text from reuters SGML articles
##
## SYNOPSIS:  get_reuters_body.sh < infile > outfile
##
################################################################################

START_TAG="<BODY>"
END_TAG="</BODY>"
IGNORE_LINE="[Rr][Ee][Uu][Tt][Ee][Rr]"   #ignore lines that just say Reuter
IGNORE_LINE_LEN=6

START_FOUND=1
while read LINE
do
    if [ `echo "$LINE" | grep "$START_TAG" | wc -l` -eq 1 ]
    then
        #start tag found, take everything on the line after the tag
        START_FOUND=0
        #@@@currently hard coding start tag due to quote madness!
        DATA=`echo "$LINE" | sed -e 's/\\(.*\\)\\(<BODY>\\)\\(.*\\)/\\3/g'`
        echo "$DATA"
    elif [ `echo "$LINE" | grep "$END_TAG" | wc -l` -eq 1 ]
    then
        START_FOUND=1
        #only junk appears before the tag so ignore this line
        #DATA=`echo "$LINE" | sed -e \'s/\\(^.*\\)\\(${END_TAG}\\)/\\1/g\'`
        #echo "$DATA"
    else
        if [ $START_FOUND -eq 0 ]
        then
            if [ `expr "$LINE" : "$IGNORE_LINE"` -ne $IGNORE_LINE_LEN ]
            then
                echo "$LINE"
            fi
        fi
    fi
done
