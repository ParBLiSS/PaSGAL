PaSGAL
========================================================================

PaSGAL (**Pa**rallel **S**equence to **G**raph **Al**igner) is designed to accelerate the local sequence alignment of sequences to directed acyclic sequence graphs (e.g., variation graphs, splicing graphs). The underlying algorithm is a parallelization of dynamic programming procedure for sequence to DAG alignment. This routine is also commonly referred as *partial order alignment*. Because computing exact alignments is compute intensive, PaSGAL uses Advanced Vector Extensions (AVX) SIMD instructions and OpenMP to achieve high alignment performance on CPUs with wide SIMD width and multiple cores. Given a set of query sequences (e.g., long PacBio/ONT or short Illumina reads) and a reference DAG, PaSGAL produces an highest scoring optimal local alignment for each query sequence along a path in the graph. The algorithm uses inter-task parallelization, i.e., different execution units operate on different queries at a time. Further details about the algorithm and evaluation are available in our [preprint](https://www.biorxiv.org).


## Dependencies

- `cmake` version >= 3.1
- A C++ compiler with c++11 support, e.g., GNU `g++` (version 4.9+) or Intel `icpc` (version 15+)
- [Google Protobuf](https://github.com/protocolbuffers/protobuf) library

## Download and compile

```sh
git clone <GITHUB_URL>
```

Compile using the cmake tool.

```sh
mkdir build_directory && cd build_directory
cmake <options, see the NOTE below> ../PaSGAL
make -j4
```

NOTE: 1) `-DPROTOBUF_DIR=<path>` should point to installation directory of the protobuf library. 2) In case avx512 feature is not available on the CPU being used, `-DSIMD_SUPPORT=<avx512/avx2/none>` should be specified accordingly. 3) Cmake will automatically look for default C/C++ compilers. To modify the default, user can set variables `-DCXX_COMPILER=<path to C++ compiler>` and `-DC_COMPILER=<path to C compiler>` if needed. After the compilation completes, expect an executable `PaSGAL` in your build\_directory. 

## Usage

* Produce help page
```sh
PaSGAL -h
```

* Align set of query sequences against a reference DAG (in .vg format):
```sh
PaSGAL -m vg -r graph.vg -q reads.fq -o outputfile -t 24
```

* Map set of query sequences against a reference DAG (in .txt format):
```sh
PaSGAL -m txt -r graph.txt -q reads.fq -o outputfile -t 24
```

**Output file format:** The output is tab-delimited with each line consisting of query id, query length, 0-based start offset, end offset, strand, reference graph start, reference graph end, alignment score and cigar string. The reference offsets are indicated as tuples of the corresponding vertex id and character offset in it.

## Graph input format
PaSGAL currently accepts a DAG in two input formats: `.vg` and `.txt`. `.vg` is a protobuf graph format, defined by [VG](https://github.com/vgteam/vg) tool developers [here](https://github.com/vgteam/vg/wiki/File-Formats). `.txt` is a simple human readable format. The first line indicates the count of total vertices (say *n*). Each subsequent line contains information of vertex *i*, 0 <= *i* < *n*. The information conveys out-neighbor vertex ids followed by its DNA sequence. For example, the following graph is a chain of four vertices:

```sh
4
1 AC
2 GT
3 GCCTG
CT
```

## An example run

Sample input data is available [data](data) folder to do a quick test run. 

```sh
$ ./PaSGAL -m vg -r BRCA1.vg -q short_reads.fa -t 16 -o output.txt
```

Expect output log in the following format during execution:

```sh
$ PaSGAL -r data/BRCA1.vg -m "vg" -q data/reads.fa -t 16 -o output.txt
--------
Assert() checks     ON
AVX SIMD support    ON (AVX512)
VTUNE profiling     OFF
--------

INFO, psgl::parseandSave, reference file = data/BRCA1.vg (in vg format)
INFO, psgl::parseandSave, query file = data/reads.fa
INFO, psgl::parseandSave, output file = output.txt
INFO, psgl::parseandSave, thread count = 16
INFO, psgl::parseandSave, scoring scheme = [ match:1 mismatch:1 ins:1 del:1 ]
....
....
INFO, psgl::main, run finished
```

## TODOs

* Support semi-global alignment mode
* Support affine gap penalty
* Support .gfa input format for graphs
