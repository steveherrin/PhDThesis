#!/bin/bash

INFILE=$1
OUTFILE="${INFILE%.*}_embedded.pdf"

AUTHOR="Steven Herrin"
TITLE="Double beta decay in xenon-136"

gs -dNOPAUSE -dBATCH -dAutoRotatePages=/None -dEmbedAllFonts=true -sDEVICE=pdfwrite -o ${OUTFILE} -c ".setpdfwrite <</NeverEmbed [ ]>> setdistillerparams" -f $INFILE
exiftool -overwrite_original -Title="${TITLE}" -Author="${AUTHOR}" -Creator="" ${OUTFILE}

echo "Author set to: ${AUTHOR}"
echo "Title  set to: ${TITLE}"
