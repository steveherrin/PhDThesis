#!/bin/bash

INFILE=$1
OUTFILE="${INFILE%.*}_final.pdf"

AUTHOR="Steven Herrin"
TITLE="Double beta decay in xenon-136"

# Embed all the fonts, even the universal ones (with the .setpdfwrite part)
gs -dNOPAUSE -dBATCH -dAutoRotatePages=/None -dEmbedAllFonts=true -sDEVICE=pdfwrite -o ${OUTFILE} -c ".setpdfwrite <</NeverEmbed [ ]>> setdistillerparams" -f $INFILE
# Add author and title information, and strip out creator information
exiftool -overwrite_original -Title="${TITLE}" -Author="${AUTHOR}" -Creator="" ${OUTFILE}

echo "Author set to: ${AUTHOR}"
echo "Title  set to: ${TITLE}"
echo "Final PDF written to: ${OUTFILE}"
