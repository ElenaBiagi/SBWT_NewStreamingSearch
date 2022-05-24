# SBWT

This is the code for the paper [Succinct k-mer Set Representations Using Subset Rank Queries on the Spectral Burrows-Wheeler Transform (SBWT)](https://www.biorxiv.org/content/10.1101/2022.05.19.492613v1). The repository includes implementations of the various SBWT variants described in the paper. Note that contrary to many other k-mer membership data structures, our code is not aware of DNA reverse complements. That is, it considers a k-mer and its reverse complement as separate k-mers.

This construction algorithm is based on the lightning-fast [k-mer counter KMC](https://github.com/refresh-bio/KMC). Our code links directly to the KMC binaries. We have made slight changes to the KMC codebase to make this possible. Our modified version of KMC is included in this repository.

We are currently actively working on the code. Top items on the to-do list are the following:

* Support for other file formats than FASTA.
* Support from multiple input files.
* Saving space by indexing only canonical k-mers.

# Compiling

```
git submodule init
git submodule update

# Build the KMC components
cd KMC
make -j4
cd ..

# Build the SBWT code
cd build
cmake .. -DMAX_KMER_LENGTH=32
make -j4
```

Change the parameter `-DMAX_KMER_LENGTH=32` to increase the maximum allowed k-mer length, up to 255.

**Troubleshooting**: If you run into problems involving the `<filesystem>` header, you probably need to update your compiler. The compiler `g++-8` should be sufficient. Install a new compiler and direct CMake to use it with the `-DCMAKE_CXX_COMPILER` option. For example, to set the compiler to `g++-8`, run CMake with the option `-DCMAKE_CXX_COMPILER=g++-8`. 

Note: the Elias-Fano variants make use of the `_pext_u64` instruction in the BMI2 instruction set. Older CPUs might not support this instruction. In that case, we fall back to a simple software implementation, which will ruin the performance of the Elias-Fano variants (those whose variant name starts with "mef").

# Index construction

Below is the command to build the SBWT for input data `example_data/coli3.fna` provided in this repository, with k = 30. The index is written to the file `index.sbwt`.

```
./build/bin/sbwt build -i example_data/coli3.fna -o index.sbwt -k 30
```

This builds the default variant, which is the plain matrix SBWT. Other variant can be specified with the `--variant` option.
The list of all command line options and parameters is below:

```
Construct an SBWT variant.
Usage:
  build [OPTION...]

  -i, --in-file arg           The input sequences as a FASTA file.
  -o, --out-file arg          Output file for the constructed index.
  -k, --kmer-length arg       The k-mer length.
      --variant arg           The SBWT variant to build. Available 
                              variants: plain-matrix rrr-matrix mef-matrix 
                              plain-split rrr-split mef-split plain-concat 
                              mef-concat plain-subsetwt rrr-subsetwt 
                              (default: plain-matrix)
      --no-streaming-support  Save space by not building the streaming 
                              query support bit vector. This leads to 
                              slower queries.
  -t, --n-threads arg         Number of parallel threads. (default: 1)
  -m, --ram-gigas arg         RAM budget in gigabytes (not strictly 
                              enforced). Must be at least 2. (default: 2)
      --temp-dir arg          Location for temporary files. (default: .)
  -h, --help                  Print usage

```

# Running queries

To query for existence of all k-mers in an index for all sequences in a fasta-file, run the following command:

```
./build/bin/sbwt search -i index.sbwt -q example_data/queries.fna -o out.txt
```

This prints for each query of length n in the input a line containing n-k+1 space-separated integers, which are the ranks of the columns representing the k-mer in the index. If the k-mer is not found, -1 is printed. If the index was built with `--streaming-support`, the faster streaming query algorithm is automatically used. The full options are:

```
Query all k-mers of all input reads.
Usage:
  ./build/bin/sbwt search [OPTION...]

  -o, --out-file arg    Output filename.
  -i, --index-file arg  Index input file.
  -q, --query-file arg  The query in FASTA format.
  -h, --help            Print usage
```

# For developers: building and running the tests 

```
git submodule init
git submodule update

# Build googletest
cd googletest
mkdir build
cd build
cmake ..
make

# Build SBWT
cd ../../build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DMAX_KMER_LENGTH=32
make
```

This will build the executable `./build/bin/sbwt_tests`. Make sure to run the test executable from the root of the repository, or otherwise it will not find the example data in ./example_data.
