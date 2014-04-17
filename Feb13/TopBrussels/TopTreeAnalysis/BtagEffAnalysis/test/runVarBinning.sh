export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../..:../../TMVA/lib/
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:../..:../../TMVA/lib/

g++ -m64 -L../../ -L../../TMVA/lib/ -lTopTreeAnaContent -lTopTreeAna -I `root-config --incdir` `root-config --libs` VarBinning.C -o VarBinning 

./VarBinning
