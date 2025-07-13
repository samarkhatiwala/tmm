#!/bin/tcsh

printf "%s$1.petsc," `tr '\n' ' ' < MOBI_tracer_names.txt` | sed 's/,*$//g'

