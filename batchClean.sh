

for file in ./*; do
	sed '/^#/d' <$file > zz1.txt
	sed '/^a/d' <zz1.txt > zz2.txt
	sed '/^5/d' <zz2.txt > $file
	echo "done with file";
done

