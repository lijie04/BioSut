<<<<<<< HEAD
# biosut
## Introduction
### History
I started to write [`biosut`](https://github.com/jlli6t/BioSut) since the early 2015. At first, I named this 
package `basictools`, which included in many little python codes of mine.

However, `basictools` was not well organized. Therefore, 
[`biosut`](https://github.com/jlli6t/BioSut) came out. There was several 
updates of `biosut`. Start from version `2`, biosut contains six modules: 
`gt_file`, `gt_path`, `io_seq`, `alter_seq`, `io_bam`, `alter_bam`.

Later, I realized that there are some functionality redundancy, therefore,
start from version `3`, there is a big re-formatting of the whole package 
structure. So does the coding style, more concise and efficient.

Start from version `3`, module includes: `mapper`, `aligner`, `biosys`, 
`bioseq`, and `bamer`. Actually this kind of structure is quite similar to
the versions before `2`, like version `0.*` and `1.*`, if you are interested in
the old versions, please go the 
[release page](https://github.com/jlli6t/BioSut/releases).

#### Module 1 [biosys](./1.biosys.md)
This module is a collection of quite a few functions related to linux-system
operations, such as check programs are executable or not, executing linux 
command, check if files/directories are existed or not, remove file suffix,
find files with specified suffix, parse json format files and so on.

#### Module 2 [bioseq](./2.bioseq.md)
This module is a collection of quite a few functions related to sequence files
(FASTA/Q) operations, such as generating an iterator for the sequence file,
convert fastq to fasta format, count gc of input string, convert sequence file
into dictionary, assess genomes traits, select sequences according to their
length from sequence file and so on.

#### Module 3 [mapper](./3.mapper.md)
This module is collection of a few mappers such as bowtie2, bbmap, bwa, star,
salmon and so on. This module will return you the sorted bam file (alignments
were sorted by leftmost coordinates, by read name when -n is used).

#### Module 4 [aligner](./4.aligner.md)
This module is collection of a few aligners such diamond, blast and so on. This
module will return you the filter alignments table with header line.

#### Module 5 [bamer](./5.bamer.md)
This module is a collection of quite a few functions related to sam/bam file 
operations, such as covert sam to bam file, sort bam file by leftmost
coordinate, sort bam, last but not least, recover reads from bam file according 
reference sequences.

### Team
MM Team was founded Jun 2020.

### Issues and suggestions
Any issues and suggestions or enquires please use the 
[issue](https://github.com/jlli6t/BioSut/issues) portal. I will make response
as fast as I can.
=======
## Installation

Install biosut through PyPi:
```
pip3 install biosut
```

## usage
### Module biosys
```python3
from biosut.biosys import files,path
```

#### files related operations
For example, check the existance of a file/files and check if it is empty.
```python3
f = 'file.txt'
files.check_exist(f, check_empty=True)
```
Get prefix of a file.

```python3
files.get_prefix(f, include_path=True)
```

#### path related operations
For example, make sure a directory exists.
```python3
o = './test/a'
path.sure_exist(o)
```

or
```python3
new.o = path.sure_exist(o)
```


>>>>>>> cd8f5485cb8239645bb760ffe8dc9e2c8be1d37c
