
make -f MakefileSerial
L=16
folder="../data/"$L"x1/"
#folder="../data/"$L"x"$L"/"
mkdir $folder
mkdir $folder"DIST"
mkdir $folder"ERR"
mkdir $folder"PLOTS"
mkdir $folder"ANIM"
mkdir $folder"OUT"
mkdir $folder"CONFIG"
mkdir $folder"SDE"

# {14..18}}
for i in {0..1}    #field: nsteps goes in here           
#for i in 0
do 
    for b in 0   # processor 
    do
	#sqsub -q serial -n 1 --mpp=1G -r 7d -o $folder"OUT/"out"h"${i}"b"${b}"beta"${beta}.log ./betadoub$L ${i} ${b} 
	./hybrid$L ${i} ${b} &
    done 
done





