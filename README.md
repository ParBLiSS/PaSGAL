PaSGAL
========================================================================
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)

PaSGAL (**Pa**rallel **S**equence to **G**raph **Al**igner) is designed to accelerate local sequence alignment of sequences to directed acyclic sequence graphs (DAGs), e.g., variation graphs, splicing graphs. The underlying algorithm is a parallelization of dynamic programming procedure for sequence to DAG alignment. With computing exact alignments being compute intensive, PaSGAL uses Advanced Vector Extensions (AVX) SIMD instructions and OpenMP to achieve high alignment performance on CPUs equipped with multiple cores and wide SIMD width. Given a set of query sequences (e.g., long PacBio/ONT or short Illumina reads) and a reference DAG, PaSGAL produces an highest scoring optimal local alignment for each query sequence along a path in the graph. Details about the algorithm and performance are available in our paper [below](#publication).

<p align="center">
<img src="https://1aaaa1f6-a-62cb3a1a-s-sites.googlegroups.com/site/chirgjain/readme-pasgal.visual.jpg?attachauth=ANoY7coRxb7Yqn_icioWTOKg1Pllk0_DEOnwE7DSeWeMC3N4-k2b595j4sPlOv_hP_gp8TSSEVFNrcV_8g5PPdqKPJYD6dQn3kP4u7L6PGFkG2bRZe2tge2Kl86asaYGyrwpmgBtH-PX-VbAdWf3yWbNYxRguJtRzb3do2BXvanYZyCdqeLH9KQst5k8lakcI2Bh3Vth-USKR-9g2q9Feo9VR7sBBEQQHPrrXwhjdkvOgrduRVweraw%3D&attredirects=0" height="350"/>
</p>

## Dependencies

- [cmake](https://cmake.org) version >= 3.1
- A C++ compiler with c++11 support, e.g., GNU `g++` (version 4.9+) or Intel `icpc` (version 15+)
- [Google Protobuf](https://github.com/protocolbuffers/protobuf) library

## Download and compile

The repository and external submodules can be downloaded using the recursive clone.

```sh
git clone --recursive <GITHUB_URL>
```

Next, compile the code using cmake utility:

```sh
mkdir build_directory && cd build_directory
cmake <OPTIONS> ../PaSGAL
make -j4
```

OPTIONS: 
1. `-DPROTOBUF_DIR=<path>` should provide path to installation directory of google protobuf library. 
2. If avx512 feature is not available on the CPU being used, `-DSIMD_SUPPORT=<avx512/avx2/none>` should be specified accordingly. 
3. Cmake will automatically look for default C/C++ compilers. To modify the default selection, users can set the two variables `-DCMAKE_CXX_COMPILER=<path to C++ compiler>` and `-DCMAKE_C_COMPILER=<path to C compiler>` if needed. 

After the compilation completes, expect an executable `PaSGAL` in your build\_directory. 

## Usage

* Produce help page
```sh
PaSGAL -h
```

* Align a set of query sequences against a reference DAG (in .vg format):
```sh
PaSGAL -m vg -r graph.vg -q reads.fq -o outputfile -t 24
```

* Align a set of query sequences against a reference DAG (in .txt format):
```sh
PaSGAL -m txt -r graph.txt -q reads.fq -o outputfile -t 24
```

**Output file format:** The output is tab-delimited with each line consisting of query id, query length, 0-based start offset, end offset, strand, reference graph start, reference graph end, alignment score and cigar string. The reference offsets are indicated as tuples of the corresponding vertex id and character offset in it.

## Graph input format
PaSGAL currently accepts a DAG in two input formats: `.vg` and `.txt`. `.vg` is a protobuf serialized graph format, defined by VG tool developers [here](https://github.com/vgteam/vg/wiki/File-Formats). `.txt` is a simple human readable format. The first line indicates the count of total vertices (say *n*). Each subsequent line contains information of vertex *i*, 0 <= *i* < *n*. The information in a single line conveys its zero or more out-neighbor vertex ids, followed by its non-empty DNA sequence (either space or tab separated). For example, the following graph is a directed chain of four vertices: `AC (id:0) -> GT (id:1) -> GCCGT (id:2) -> CT (id:3)`

```sh
4
1 AC
2 GT
3 GCCTG
CT
```

The first line above specifies the count of vertices as 4. The second line specifies that vertex 0 has an outgoing edge to vertex 1, and the label of vertex 0 is "AC". Similarly, the third line specifies that vertex 1 has an outgoing edge to vertex 2, and its label is "GT". 

## An example run

Sample input is available in [data](data) folder to do a quick test run. Larger data sets which we used for benchmarking in the [paper](#publication) can be accessed [here](https://alurulab.cc.gatech.edu/PaSGAL). Expect output log in the following format during execution:

```sh
$ PaSGAL -r data/BRCA1.vg -m "vg" -q data/reads.fa -t 36 -o output.txt
--------
Assert() checks     ON
AVX SIMD support    ON (AVX512)
VTUNE profiling     OFF
--------

INFO, psgl::parseandSave, reference file = data/BRCA1.vg (in vg format)
INFO, psgl::parseandSave, query file = data/reads.fa
INFO, psgl::parseandSave, output file = output.txt
INFO, psgl::parseandSave, thread count = 36
INFO, psgl::parseandSave, scoring scheme = [ match:1 mismatch:1 ins:1 del:1 ]
....
....
INFO, psgl::main, run finished
```

## TODOs

* Support semi-global alignment mode
* Support affine gap penalty
* Support .gfa input format for graphs
* Support intra-task parallelization

## <a name=“publication”></a>Publication

- **Chirag Jain, Sanchit Misra, Haowen Zhang, Alexander Dilthey and Srinivas Aluru**. "Accelerating Sequence Alignment to Graphs". *IEEE International Parallel and Distributed Processing Symposium (IPDPS) 2019*.
