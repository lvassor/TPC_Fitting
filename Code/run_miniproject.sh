#!/bin/bash
# Author: Luke Vassor ljv3618@ic.ac.uk
# Script: miniproject.sh
# Desc: Runs Luke Vassor's CMEE miniproject and compiles LaTeX report.
# Date: Mar 2019

echo Welcome to Luke Vassor miniproject, submitted for the M.Sc. CMEE course at Imperial College London.

echo Shell script will now run project.
SECONDS=0
python3 mini_project_wrangling_fitting.py
Rscript plotting.R
cd ../Write-Up/
# bibtex write_up.aux
# pdflatex write_up.tex
: | bibtex write_up.aux | grep '^!.*' -A200 --color=always                      # silences BibTex output
: | pdflatex -halt-on-error write_up.tex | grep '^!.*' -A200 --color=always     # silences pdflatex output

echo COMPLETE. Project has been run and report has been compiled. Now opening report.

if (( $SECONDS > 60 )) ; then
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo "Project ran in $minutes minute(s) and $seconds second(s)"
else
    echo "Project ran in $SECONDS seconds"
fi

xdg-open write_up.pdf