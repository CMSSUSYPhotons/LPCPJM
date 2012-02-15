
dest_dir=plots.20120215

root -q -l 'drawContour.C("bino","nojet",true)'
root -q -l 'drawContour.C("bino","1jet",true)'
root -q -l 'drawContour.C("wino","nojet",true)'
root -q -l 'drawContour.C("wino","1jet",true)'
root -q -l 'drawContour.C("bino_mNScan","nojet",true)'
root -q -l 'drawContour.C("bino_mNScan","1jet",true)'

mkdir -p $dest_dir
mv *.gif *.pdf $dest_dir
