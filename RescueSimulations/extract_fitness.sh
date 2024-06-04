
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

