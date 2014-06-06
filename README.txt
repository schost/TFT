
* Run "make" at the top-level to generate all object files. They 
  are in the build directory. 

* Run "make all" in src to generate all test files. All executables
  take an optional argument. If set to 1, runs a check. Otherwise,
  runs timings. Example: "./test_mul 1" checks multiplication, and
  "./test_mul" outputs timings. 

* After running 
      ./test_mul > ctft.dat
  and 
      ./test_middle_product > ctft_middle.dat
  running 
      gnuplot ctft.plt 
  outputs an eps file that compares ntl's FFT, this TFT and its middle 
  product.
