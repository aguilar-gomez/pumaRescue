#Save all the values not only discrete timepoints
for file in cfp_tx*out
do
grep "^15" $file |grep  "0."|uniq> $file.fitness
done


paste -d " " cfp_tx*_GR_*fitness|cut -d " " -f3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 > allfitnessGR
paste -d " " cfp_tx*_noGR_*fitness|cut -d " " -f3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 > allFitnessnoGR


cut -f1 -d " " cfp_tx*noGR_3.out.fitness > times
cut -f1 -d " " cfp_tx_recentsplit_noGR_3.out.fitness > times_recentsplit


#Current format
for file in cfp_tx*out
do
#grep "^43900\|^43911\|^43916\|^43921\|^43961\|^44011\|^44410" $file > $file.lines
grep "^153900\|^153911\|^153916\|^153921\|^153961\|^154011\|^154411\|^154911\|^155910" $file > $file.lines
done

paste -d " " cfp_tx_GR_*lines|cut -d " " -f2,5,8,11,14,17,20,23,26,29 > GR_TX_10runs
paste -d " " cfp_tx_noGR_*lines|cut -d " " -f2,5,8,11,14,17,20,23,26,29 > noGR_TX_10runs

paste -d " " cfp_tx_GR_*lines|cut -d " " -f3,6,9,12,15,18,21,24,27,30 > GR_FL_10runs
paste -d " " cfp_tx_noGR_*lines|cut -d " " -f3,6,9,12,15,18,21,24,27,30 > noGR_FL_10runs

#Fitnes in Texas
paste -d " " cfp_tx*_GR_*lines|cut -d " " -f2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59 > GR_TX_20runs
paste -d " " cfp_tx*_noGR_*lines|cut -d " " -f2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59  > noGR_TX_20runs
#Fitnes in Florida
paste -d " " cfp_tx*_GR_*lines|cut -d " " -f3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 > GR_FL_20runs
paste -d " " cfp_tx*_noGR_*lines|cut -d " " -f3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 > noGR_FL_20runs

awk '
{
    sum = 0
    for (i = 1; i <= NF; i++) {
        sum += $i
    }
    if (NF > 0) {
        average = sum / NF
        printf " %.3f\n", average
    }
}
'  GR_TX_20runs > averageGRTX
GR_FL_20runs > averageGRFL
noGR_FL_20runs > averagenoGRFL
noGR_TX_20runs > averagenoGRTX


################################################################################################################################

#Previous format
for file in c*_noGR_*out; do echo $file; tail -11850 $file>${file%%.out}_fitness.out ; done

# Define the list of line numbers to search for
line_numbers=("43900" "43911" "43916" "43921" "43961" "44011" "44410")

for file in *fitness*; do
  for number in "${line_numbers[@]}"; do
    linenumber=$(grep -n "$number" "$file" | cut -d: -f1)
    if [ -n "$linenumber" ]; then
      endline=$((linenumber + 2))
      # Extract lines and filter only floating-point numbers
      sed -n "${linenumber},${endline}p" "$file" | grep -o '[0-9]*\.[0-9]*' >> "fitness_extracted_$file"
    fi
  done
done


paste fitness_extracted_cfp_tx_noGR_* > combined_fitness_values_noGR

# Calculate the average per row
awk '{
  sum = 0;
  for (i = 1; i <= NF; i++) {
    sum += $i;
  }
  average = sum / NF;
  print average;
}' combined_fitness_values_noGR > row_averages_noGR
