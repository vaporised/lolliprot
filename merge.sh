n=100

# Zip and index vcfs
for FILE in $1/*vcf; do
  bgzip < $FILE > $FILE.temp.gz;
  bcftools index $FILE.temp.gz > $FILE.temp.gz.tbi;
done 

# Merge in groups of 100
files=($1/*.vcf.temp.gz)
for ((i=0; i < ${#files[@]}; i+=n)); do
 bcftools merge -Oz "${files[@]:i:n}" >  $1/merged_$i.vcf.temp.gz
 bcftools index $1/merged_$i.vcf.temp.gz > $1/merged_$i.vcf.temp.gz.tbi
done

# Merge the merged vcfs
bcftools merge $1/merged_*.vcf.temp.gz > merged.vcf

# Clean up
# rm *vcf.temp.gz
#rm *vcf.temp.gz.tbi
rm $1/*temp.gz
rm $1/*temp.gz.tbi
rm $1/*temp.gz.csi

echo "Done"
