#!/bin/bash

INFILE=$1
OUTFILE="${INFILE%.*}_embedded.pdf"

gs -dNOPAUSE -dBATCH -dEmbedAllFonts=true -sDEVICE=pdfwrite -o ${OUTFILE} -c ".setpdfwrite <</NeverEmbed [ ]>> setdistillerparams" -f $INFILE
