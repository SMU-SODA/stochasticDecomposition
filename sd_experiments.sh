# Command line options
#         -p string  -> problem name.
#         -i string  -> input directory where the problem SMPS files are saved.
#         -o string  -> output directory where the result files will be written.
#         -e {0,1}   -> evaluation flag which determines if the final solution will be evaluated through out of sample evalauation.
#         -d {0,1}   -> use the dual stability test.
#         -t {l,n,t} -> tolerance level to be employed.
#                        Suggested tolerance(EPSILON, SCAN_LEN) = 'l'oose (0.01, 128), 'n'omimal (0.001, 256) and 't'ight (0.0001, 512)
#         -m {0,1}   -> use multiple replication.
#         -c {0,1}   -> build and solve compromise problem.
#         -s {0,1}   -> sample increment size.

# Default values
outputDir=~/Documents/experiments/internalSampling/
inputDir=~/Documents/workspace/spAlgorithms/spInput/
evalFlag=1
dualTest=1
replicate=30
compromise=0

# prob_set=("pgp2" "cep" "4node" "baa99" "lands3" "stgorm" "fleet1" "fleet2" "ssn" "20" "baa99-20")
prob_set=("pgp2")
incSize=("1" "5" "10" "25" "50" "100")
toleranceLevel=("l" "n" "t")

# Problem set
for probName in "${prob_set[@]}"; do
    for sampleInc in "${incSize[@]}"; do
	for tol in "${toleranceLevel[@]}"; do
	    ./Debug/twoSD -i $inputDir -o $outputDir -p $probName -e $evalFlag -d $dualTest -t $tol -m $replicate -c $compromise -s $sampleInc	   
	    mv ${outputDir}twoSD/${probName} ${outputDir}twoSD/${probName}_${sampleInc}_${tol}
	done
    done
done
