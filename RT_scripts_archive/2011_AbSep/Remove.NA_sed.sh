for file in EE SS silenced induced
do
	grep NM $file | sed 's/\"//g' | cut -f2-7 > $file.list
done
