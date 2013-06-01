#!/bin/bash

INFILE=$1
OUTFILE="${INFILE%.*}_embedded.pdf"

gs -dNOPAUSE -dBATCH -dAutoRotatePages=/None -dEmbedAllFonts=true -sDEVICE=pdfwrite -o ${OUTFILE} -c ".setpdfwrite <</NeverEmbed [ ]>> setdistillerparams" -f $INFILE
