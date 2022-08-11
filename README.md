<<<<<<< HEAD
# biosut
[![PyPI](https://shields.io/pypi/v/biosut.svg)](https://pypi.org/project/biosut)

## Introduction
[**biosut**](https://github.com/jlli6t/BioSut) is a python package that
integrated some biology-related bioinformatics operating modules on sequences,
files, directories, mappings, alignments and so on. The biosut is open source
and released under the
[GNU General Public License (Version 3)](https://pypi.org/project/biosut/).

## Installation
Install the package through PyPi:
```
pip3 install biosut
```

## Documentation
Detailed usage refer to [document](docs/documentation.md)

There are 5 modules:

[biosys](docs/1.biosys.md): Linux-system operations

[bioseq](docs/2.bioseq.md): FASTA/Q sequences file operations

[mapper](docs/3.mapper.md): A collection for mappings.

[aligner](docs/4.aligner.md): A collection for alignments.

[bamer](docs/5.bamer.md): A collection for sam/bam operations.

## Dependencies
There are dependencies in need for some modules, such as bammer:

[samtools](http://www.htslib.org/)

## Bugs and questions
For any bugs or questions please use
[Issue](https://github.com/jlli6t/BioSut/issues) portal.

## Copyright
Copyright 2017-2021 Jie Li. See [LICENSE](./LICENSE) for further details.
=======
**biosut ("biology suite tool") is a packge containing some biology-related bioinformatics operations on sequences, file**

## Installation
Install biosut through PyPi:
```
pip3 install biosut
```

## Usage
### Module biosys
```
from biosut.biosys import files,path
```
#### files related operations
For example, check the existance of a file/files and check if it's empty.
```
f = 'file.txt'
files.check_exist(f, check_empty=True)
```
Get prefix of a file.
```
files.get_prefix(f, include_path=True)
```

#### path related operations
For example, make sure a directory exists.
```
o = './test/a'
path.sure_exist(o)
```
or
```
new.o = path.sure_exist(o)
```

# Summary


>>>>>>> cd8f5485cb8239645bb760ffe8dc9e2c8be1d37c
