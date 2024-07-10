NGRAINS=30
MAXERROR=0.01
for NXNYNZ in 32
do
  FNAME="$NXNYNZ"Cube"$NGRAINS"Grains
  eval "neper -T -n "$NGRAINS" -id 5 -periodicity 1 -morpho 'gg,aspratio(2,1,1)' -oridescriptor euler-bunge -morphooptistop eps="$MAXERROR" -format tesr,vtk -o "$FNAME" -tesrsize "$NXNYNZ":"$NXNYNZ":"$NXNYNZ" -tesrformat ascii"
  eval "python convert_neper_to_LApx.py "$FNAME""
done
