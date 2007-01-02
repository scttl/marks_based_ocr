#!/bin/sed -f
###############################################################################
##
## FILE: srilm_char_ngram.sed
##
## DESCRIPTION: Quick sed shell script to process input text files for use as a
##              character N-gram model for SRILM.  Essentially converts spaces
##              to the '_' character, and adds a space after each character
##
## SYNOPSIS: ./srilm_char_ngram.sed  < infile > outname
##
###############################################################################

#change spaces to underscores
y/ /_/

#add a space character to each character
s/\(.\)/\1 /g
