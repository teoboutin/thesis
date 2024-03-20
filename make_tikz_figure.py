#!/usr/bin/python3



import sys

import os
import shutil

if __name__=="__main__":
    
    tikzfile = sys.argv[1]
    openpdf = sys.argv[2] if len(sys.argv)>2 else ""
    
    cwd = os.getcwd()
    os_home = os.environ.get('HOME')
    
    
    
    pdfoutput_filename=cwd + "/" + tikzfile.replace(".tex", ".pdf")
    
    
    print(os_home)
    print(tikzfile)
    
    
    tempfiledir=fr"{os_home}/.tikz"
    try:
        os.mkdir(tempfiledir)
    except FileExistsError:
        print("exist")
        
    tempfilename=fr"{tempfiledir}/figure"
    
    
    with open(tikzfile, "r") as content:
        print(content)

        with open(fr"{tempfilename}.tex", "w") as tempfile:
            tempfile.write("\\documentclass[a4paper]{article}\n")
            
            tempfile.write("\\usepackage{tikz}\n")
            tempfile.write("\\usepackage{pgfplots}\n")
            tempfile.write("\\usepackage{bm}\n")
            tempfile.write("\\pgfplotsset{compat=1.16}\n")
            
            line=content.readline()
            while line.startswith("%"):
                print(line.rstrip("\n"))
                if line.startswith("%RequestPackage") or line.startswith("%RequirePackage"):
                    request, package=line.rstrip("\n").split(" ")
                    print(fr"requested package '{package}'")
                    tempfile.write("\\usepackage{%s}\n" % (package))
                line=content.readline()
            
            
            tempfile.write("\\begin{document}\n")
            
            tempfile.write("\\hoffset=-1in\n")
            tempfile.write("\\voffset=-1in\n")
            tempfile.write("\\setbox0\\hbox{\n")
            
            
            
            tempfile.write(line)
            for l in content.readlines():
                tempfile.write(l)
                
            tempfile.write("}\n")
            
            tempfile.write("\\pdfpageheight=\\dimexpr\\ht0+\\dp0\\relax\n")
            tempfile.write("\\pdfpagewidth=\\wd0\n")
            tempfile.write("\\shipout\\box0\n")
            
            
            
            #  tempfile.write("\\end{document}")
            tempfile.write("\\stop")
    
    os.system(fr"cd {tempfiledir} && pdflatex -halt-on-error {tempfilename}.tex")
    #  stream = os.popen('pdflatex temp')
    #  output = stream.read()
    #  print(output)
    
    print("copy ", fr"{tempfilename}.pdf", "to", pdfoutput_filename)
    shutil.copyfile(fr"{tempfilename}.pdf", pdfoutput_filename)
    #  stream = os.popen('evince temp.pdf')
    if openpdf=="open":
        os.system(fr"evince {tempfilename}.pdf")
    
