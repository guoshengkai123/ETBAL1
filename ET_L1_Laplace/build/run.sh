cmake ..
make
cd bin

cd output_temp
rm -r ./*
cd ..

cd simulate_data
rm -r ./*
cd ..

cd simulate_data_noOne
rm -r ./*
cd ..


rm -rf error_Smoothing.txt
rm -rf memory.txt
rm -rf error_begin.txt

gcc run_testmemory.c -o run_testmemory

for((i=0;i<100;i++));do

    echo "--------------------------------- "$[i]" -------------------------------------------"

    ./run_testmemory $[i]

done

cd ..


