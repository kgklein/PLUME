#Name of input file:
input='example_map_par'
#data subdirectory; see &params name list
#from input file.
dataName='example'
cd ../../
#Always make sure the data directory is present
mkdir -p data
#Always make sure the requested sub-directory
#indicated with the 'dataName' variable
#is present
mkdir -p data/$dataName
#Execute PLUME wiht the selected input file
./plume.e inputs/example/$input'.in'
#Return to original directory
cd inputs/example/
