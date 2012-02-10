
dest_dir="table_20120209/multiChannel"
mkdir -p $dest_dir

rm -f *.table

./makeTable.sh bino 1jet multiChannel
./makeTable.sh bino nojet multiChannel
./makeTable.sh wino 1jet multiChannel
./makeTable.sh wino nojet multiChannel
./makeTable.sh bino_mNScan 1jet multiChannel
./makeTable.sh bino_mNScan nojet multiChannel
mv *.table $dest_dir


#./makeTable.sh bino 1jet bin0
#mv *.table ${dest_dir}/bin0

#./makeTable.sh bino 1jet bin1
#mv *.table ${dest_dir}/bin1

#./makeTable.sh bino 1jet bin2
#mv *.table ${dest_dir}/bin2

#./makeTable.sh bino 1jet bin3
#mv *.table ${dest_dir}/bin3

#./makeTable.sh bino 1jet bin4
#mv *.table ${dest_dir}/bin4

#./makeTable.sh bino nojet bin0
#mv *.table ${dest_dir}/bin0

#./makeTable.sh bino nojet bin1
#mv *.table ${dest_dir}/bin1

#./makeTable.sh bino nojet bin2
#mv *.table ${dest_dir}/bin2

#./makeTable.sh bino nojet bin3
#mv *.table ${dest_dir}/bin3

#./makeTable.sh bino nojet bin4
#mv *.table ${dest_dir}/bin4
