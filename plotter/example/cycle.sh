#Plotting Routine, *.plt
input='unstable_check'

#Run the plotting routine.
#Requires a local installation of gnuplot
gnuplot $input.plt

#The *.plt file creates a *.tex file,
#Which we will use latex (or pdflatex)
#to make a portable digital file.
for files in $input*.tex
do
    jj=$(basename $files .tex)
    echo $jj
    pdflatex $jj.tex
    #latex $jj.tex
    #dvips $jj.dvi
    #ps2pdf $jj.ps

    #If conversion to *.png is necesssary
    #Requires ImageMagick
    
    #convert -density 300 $jj.pdf -trim $jj.png
    #eog $jj.png
done

#You don't need those extra files.
rm -f *.log
rm -f *.eps
rm -f *.aux
rm -f *.ps
rm -f *.dvi
rm -f *.tex

