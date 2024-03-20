#!/bin/bash

pdftk doc_jury.pdf cat 2-3 output PV.pdf

pdftk doc_jury.pdf cat 4-6 output rapport.pdf

pdftk doc_jury.pdf cat 1 7-10 output guide.pdf
