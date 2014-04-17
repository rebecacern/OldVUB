texfiles=$(ls 2DBtagEff_XS*.tex)

for i in $texfiles;do 

  echo "Compiling $i"

   pdflatex --file-line-error --shell-escape --synctex=1 $i

done

 #rm -rfv AnalysisResults

 mkdir AnalysisResults

 mv -v 2DBtagEff_XS*.pdf AnalysisResults

 rm -v 2DBtagEff_XS*
