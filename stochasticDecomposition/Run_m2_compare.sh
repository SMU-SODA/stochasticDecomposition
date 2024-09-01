#!/bin/bash

# Define the output file for the results
OUTPUT_FILE="1.csv"

# Create the header for the CSV file
echo "Problem Name,Iteration,Loose Tolerance (t = l),Nominal Tolerance (t = n),Tight Tolerance (t = t)" > $OUTPUT_FILE


# Define the problems and -t values
problems=("pgp2") # for testing purposes
#problems=("pgp2" "baa99" "baa99-20" "cep" "lands" "lands2" "lands3" "retail" "ssn" "storm" "4node" "20")
t_values=("l" "n" "t")

# Function to run sultan command and extract lower bound estimates
run_sultan() {
    local problem=$1
    local iterations=$2
    local loose_bounds=()
    local nominal_bounds=()
    local tight_bounds=()

    for t in "${t_values[@]}"; do
        local result=$(sultan -p "$problem" -i /Users/sultan/S2/spAlgorithms/spInput/ -o /Users/sultan/S2/spAlgorithms/spOutput -e 0 -d 1 -t "$t" -m "$iterations" -c 0 2>&1)

        # Check if the sultan command failed
        if [[ $? -ne 0 ]]; then
            echo "Error running sultan command for problem $problem with -t $t"
            echo "$result"
            continue
        fi

        # Loop through each replication and extract the lower bound estimate
        for ((i = 1; i <= iterations; i++)); do
            lower_bound=$(echo "$result" | awk -v rep="Replication-$i" '$0 ~ rep {getline; while ($0 !~ /Lower bound estimate/ && getline) {}; if ($0 ~ /Lower bound estimate/) {print $NF}}')
            if [[ "$t" == "l" ]]; then
                loose_bounds+=("$lower_bound")
            elif [[ "$t" == "n" ]]; then
                nominal_bounds+=("$lower_bound")
            elif [[ "$t" == "t" ]]; then
                tight_bounds+=("$lower_bound")
            fi
        done
    done

    # Write the results to the CSV file
    for ((i = 0; i < iterations; i++)); do
        if (( i == 0 )); then
            echo "$problem,$((i + 1))st,${loose_bounds[i]},${nominal_bounds[i]},${tight_bounds[i]}" >> $OUTPUT_FILE
        else
            echo ",$((i + 1))nd,${loose_bounds[i]},${nominal_bounds[i]},${tight_bounds[i]}" >> $OUTPUT_FILE
        fi
    done
}

# Run the sultan command for each problem
for problem in "${problems[@]}"; do
    run_sultan "$problem" 2
done

echo "Results have been saved to $OUTPUT_FILE"

