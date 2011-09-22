
root -b -q -l 'makeExclusionHist.C+("bino","375","met100","1jet","1")'
root -b -q -l 'makeExclusionHist.C+("bino","375","met100","nojet","1")'
root -b -q -l 'makeExclusionHist.C+("bino","375","met100","1jet","0")'
root -b -q -l 'makeExclusionHist.C+("bino","375","met50","1jet","1")'
root -b -q -l 'makeExclusionHist.C+("wino","375","met100","1jet","1")'
root -b -q -l 'makeExclusionHist.C+("bino","150","met100","1jet","1")'
root -b -q -l 'makeExclusionHist.C+("wino","150","met100","1jet","1")'

root -b -q -l 'drawContour.C("bino","375","met100","1jet","1")'
root -b -q -l 'drawContour.C("bino","375","met100","nojet","1")'
root -b -q -l 'drawContour.C("bino","375","met100","1jet","0")'
root -b -q -l 'drawContour.C("bino","375","met50","1jet","1")'
root -b -q -l 'drawContour.C("wino","375","met100","1jet","1")'
root -b -q -l 'drawContour.C("bino","150","met100","1jet","1")'
root -b -q -l 'drawContour.C("wino","150","met100","1jet","1")'

